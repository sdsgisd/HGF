//
//  HGF.h
//  
//
//  Created by Fang Da on 10/27/14.
//
//  Eddited by Sadashige Ishida 2017.

#ifndef __MultiTracker__HGF__
#define __MultiTracker__HGF__

#include <iostream>
#include <memory>

#include <Eigen/SparseCore>

#include "eigenheaders.h"
#include "surftrack.h"

class Sim;
class Scenes;

class HGF : public LosTopos::SurfTrack::SolidVerticesCallback, public LosTopos::T1Transition::VelocityFieldCallback, public LosTopos::SurfTrack::MeshEventCallback
{
    friend class Sim;
    friend class Scenes;

public:
    //////////////////////////
    //Simulation parameters.
    //////////////////////////
    
    bool accel=true;//The flow is parabolic if it is false.
    
    double damp=1.0;
    double surface_tension=1.0;
    
    double smooth=0.1;
    bool smooth_from_neighbors=false;
    //The velocity of each vertex is smoothed using its neighborhoods,
    //if it is true.
    
    bool save_mesh=false;
    
    bool sparse=false;
    
    bool implicit_scheme=false;
    
    //Used for logging information.
    bool logging_geometry=true;
    bool write_geometry=true;
    
    bool logging_time=false;
    bool logging_detailed_time=false;
    double time_dAdx,time_volume_correction,time_los_topos;
    double ave_time_dAdx,ave_time_volume_correction,ave_time_los_topos;

    //Used for energy preservation.
    //Recommended to be all false,
    //since none of them does not work well
    //at this moment.
    bool total_energy_preservation=false;
    bool local_energy_preservation=false;
    bool energy_preserve_per_vertex=false;
    bool energy_preserve_with_max_velocity=false;

    //////////////////////////
    //Simulation parameters
    //dynamically used for a time step.
    //////////////////////////
    bool with_gravity=false;
    bool add_velocity=false;
    
    bool blowing_bubble0=false;
    bool absorbing_bubble0=false;
    
    bool move_right=false;
    bool move_left=false;
    bool do_stepScenes=false;
    
    bool bursting=false;

public:
    int num_bubbles;
    int num_region;
    Eigen::VectorXd initialVolumes;
    
    Eigen::MatrixXd Constrained_V;//Constrained vertices for rendering purposes.
    Eigen::MatrixXi Constrained_F;//Constrained triangles for rendering purposes.

    std::vector<std::unique_ptr<Eigen::Vector3d>> velocity_per_bubble;
    Eigen::Vector3d gravity_vec;

    HGF(const std::vector<LosTopos::Vec3d> & vs, const std::vector<LosTopos::Vec3st> & fs, const std::vector<LosTopos::Vec2i> & ls, const std::vector<size_t> & constrained_vertices= std::vector<size_t>(),  const std::vector<Vec3d> & constrained_positions = std::vector<Vec3d>(),const int num_bubbles=-1,const bool detect_boundary=false);
    
    ~HGF();
    
protected:
    LosTopos::SurfTrack * m_st;
 
    bool delta_dn_tq_juctions=true;
    //Special care for t/q-junction
    //when determining the amount of correction
    //for volume preservation.
    
    int region_offset=1;//Offset for open regions
    int AIR;//Maximum index of open regions.
    //Basically, AIR=region_offet-1;
    
    // Velocities of vertices
    LosTopos::NonDestructiveTriMesh::VertexData<Vec3d> * m_v;
    
    // constrained vertices
    std::vector<size_t> m_constrained_vertices;
    std::vector<Vec3d> m_constrained_positions;

public:
    //////////////////////////
    //Main functions
    //////////////////////////
    double step(double dt);

    void easy_orientation();
    //Set each triangle label (i,j) to be i>j.
    
    void set_parameters();
    
    void detect_boundary_vertices();
    //Currently not used.
    //Incompatible with a mesh structure with ghost vertices and open regions.

protected:
    //////////////////////////
    //For updating the state
    //////////////////////////
    double actual_dt;
    
    void stepHGF(double dt);
    
    void set_intermediate_symplectic_euler(double dt,Eigen::MatrixXd &U,Eigen::MatrixXi& F,Eigen::MatrixXd &dUdt,Eigen::MatrixXd &NewU);
    
    void set_intermediate_implicit(double dt,Eigen::MatrixXd &U,Eigen::MatrixXi& F,Eigen::MatrixXd &dUdt,Eigen::MatrixXd &NewU);
    
    void smoothVelocity_neighborhood();
    
    //Correct volume of each closed region.
    //This is an extension of [Muller 2009]"Fast and Robust Tracking of Fluid Surfaces" for multiple regions.
    void correct_volume(Eigen::MatrixXd &targetU,Eigen::MatrixXi& F);
    void computeDelta_d_n(const Eigen::MatrixXd &targetU, const Eigen::MatrixXi& F,Eigen::MatrixXd &Delta_d_n);

public:
    ////////////////////
    //Utilities
    ////////////////////
    void writeObj_FaceLabel_constrainedVertices(const bool with_imaginary_vertices=false);
    void geometric_information(double time);
    void write_mesh();
    void write_film_mesh(std::string=".");
    void write_constrained_mesh(std::string=".");
    void face_edge_vertex_style(Eigen::MatrixXi &FE,Eigen::MatrixXi &EV,Eigen::MatrixXd &V);
    //Make data structure where each row i of FE indicates the edges of the i-th triangle, and each row j of EV indicates the vertices of the j-th edge.

    void volumes_and_areas_sparse(const Eigen::MatrixXd &targetU,Eigen::VectorXd &volumes,Eigen::SparseMatrix<double>& area_matrix );
    void volumes_and_areas(const Eigen::MatrixXd &targetU,Eigen::VectorXd &volumes, Eigen::MatrixXd& area_matrix );
    void volumes_and_areas(Eigen::VectorXd &volumes, Eigen::MatrixXd& area_matrix );

    void total_area(double time,bool write=true);
    void tq_junc_acurracy(double time, bool write=true);
    
    void vertex_areas(std::vector<double> &v_areas);
    double total_energy();

    void test_update_mesh_via_LosTopos();

public:
    const LosTopos::SurfTrack * surfTrack() const { return m_st; }
    LosTopos::SurfTrack * surfTrack()       { return m_st; }
    const LosTopos::NonDestructiveTriMesh & mesh() const { return m_st->m_mesh; }
    LosTopos::NonDestructiveTriMesh & mesh()       { return m_st->m_mesh; }
    
    Vec3d pos(size_t v) const { return vc(m_st->pm_positions[v]); }
    const Vec3d & vel(size_t vi) const { return (*m_v)[vi]; }
    Vec3d & vel(size_t vi)       { return (*m_v)[vi]; }
    std::vector<Vec3d> vel() const {
        int nv=mesh().nv();
        std::vector<Vec3d> vels(nv);
        for (size_t i = 0; i <  nv; i++)
        { vels[i] = vel(i);}
        return vels;
    }
    VecXd    velv() const {
        int nv=mesh().nv();
        VecXd vels = VecXd::Zero(nv* 3);
        for (size_t i = 0; i < nv; i++)
        { vels.segment<3>(i * 3) = vel(i); }
        return vels;
    }

    const std::vector<size_t> & constrainedVertices() const { return m_constrained_vertices; }
    std::vector<size_t> & constrainedVertices()       { return m_constrained_vertices; }
    const std::vector<Vec3d> & constrainedPositions() const { return m_constrained_positions; }
    std::vector<Vec3d> & constrainedPositions()       { return m_constrained_positions; }

protected:
    //////////////////////////
    //For acquiring information on mesh connectivity.
    //////////////////////////
    
    // convenient mesh topology query
    size_t edge_other_vertex(size_t e, size_t v) const { LosTopos::Vec2st edge = mesh().m_edges[e]; return edge[0] == v ? edge[1] : edge[0]; }
    // convenient mesh geometry query
    double face_area(size_t f)    const { return m_st->get_triangle_area(f); }
    double vert_area(size_t v)    const { double a = 0; for (size_t i = 0; i < mesh().m_vertex_to_triangle_map[v].size(); i++) a += face_area(mesh().m_vertex_to_triangle_map[v][i]); return a / 3; }
    
    Vec3d  face_outward_normal(size_t f)  const { LosTopos::Vec3st t = mesh().m_tris[f]; Vec3d n = (pos(t[1]) - pos(t[0])).cross(pos(t[2]) - pos(t[0])).normalized(); LosTopos::Vec2i l = mesh().get_triangle_label(f); if (l[0] < l[1]) return -n; else return n; }

    //////////////////////////
    //For mesh update via LosTopos
    //////////////////////////
protected:
    // SurfTrack::SolidVerticesCallback method
    bool            generate_collapsed_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos);
    bool            generate_split_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos);
    LosTopos::Vec3c generate_collapsed_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1);
    LosTopos::Vec3c generate_split_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1);
    bool            generate_edge_popped_positions(LosTopos::SurfTrack & st, size_t oldv, const LosTopos::Vec2i & cut, LosTopos::Vec3d & pos_upper, LosTopos::Vec3d & pos_lower);
    bool            generate_vertex_popped_positions(LosTopos::SurfTrack & st, size_t oldv, int A, int B, LosTopos::Vec3d & pos_a, LosTopos::Vec3d & pos_b);
    bool            solid_edge_is_feature(const LosTopos::SurfTrack & st, size_t e);
    
    // T1Transition::VelocityFieldCallback methods
    LosTopos::Vec3d sampleVelocity(LosTopos::Vec3d & pos);
    bool sampleDirectionalDivergence(const LosTopos::Vec3d & pos, const LosTopos::Vec3d & dir, double & output);
    
    // SurfTrack::MeshEventCallback
    void pre_collapse(const LosTopos::SurfTrack & st, size_t e, void ** data);
    void post_collapse(const LosTopos::SurfTrack & st, size_t e, size_t merged_vertex, void * data);
    
    void pre_split(const LosTopos::SurfTrack & st, size_t e, void ** data);
    void post_split(const LosTopos::SurfTrack & st, size_t e, size_t new_vertex, void * data);
    
    void pre_flip(const LosTopos::SurfTrack & st, size_t e, void ** data);
    void post_flip(const LosTopos::SurfTrack & st, size_t e, void * data);
    
    void pre_t1(const LosTopos::SurfTrack & st, size_t v, void ** data);
    void post_t1(const LosTopos::SurfTrack & st, size_t v, size_t a, size_t b, void * data);
    
    void pre_facesplit(const LosTopos::SurfTrack & st, size_t f, void ** data);
    void post_facesplit(const LosTopos::SurfTrack & st, size_t f, size_t new_vertex, void * data);
    
    void pre_snap(const LosTopos::SurfTrack & st, size_t v0, size_t v1, void ** data);
    void post_snap(const LosTopos::SurfTrack & st, size_t v_kept, size_t v_deleted, void * data);
    
    void pre_smoothing(const LosTopos::SurfTrack & st, void ** data);
    void post_smoothing(const LosTopos::SurfTrack & st, void * data);
    
    std::ostream & log() { static std::stringstream ss; return ss; }

};

#endif /* defined(__MultiTracker__HGF__) */
