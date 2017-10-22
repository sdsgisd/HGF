//
//  HGF.cpp
//
//
//  Created by Fang Da on 10/27/14.
//
//  Modified by Sadashige Ishida in 2017.

#include "vec.h"
#include <iostream>
#include <fstream>
#include <chrono>

#include "HGF.h"
#include "Options.h"
#include "LosTopos/LosTopos3D/subdivisionscheme.h"

void HGF::set_parameters(){
    
    delta_dn_tq_juctions=1;

    accel=Options::boolValue("accel");
    damp= Options::doubleValue("damping-coef");//1.0;0.999 ;0.995;1-dt*damp_scale.
    surface_tension=Options::doubleValue("surface_tension");;

    smooth=Options::doubleValue("smoothing-coef");
    const double eps=0.0001;
    smooth_from_neighbors=smooth<eps?false:true;
    
    save_mesh=Options::boolValue("save_mesh");
    
    sparse=Options::boolValue("sparse");
    
    implicit_scheme=Options::boolValue("implicit_scheme");

    //Used for logging information.
    logging_geometry=Options::boolValue("logging_geometry");
    write_geometry=Options::boolValue("write_geometry");;
    logging_time=Options::boolValue("logging_time");;
    logging_detailed_time=Options::boolValue("logging_detailed_time");;
    
    if(logging_detailed_time){
        time_dAdx=0.0;
        time_volume_correction=0.0;
        time_los_topos=0.0;
    }

    with_gravity=Options::boolValue("with_gravity");
    add_velocity=Options::boolValue("add_velocity");
    double gravity_scale= 0.02;//0.007;
    gravity_vec<<0,0,-9.8*gravity_scale;

    //Used for energy preservation.
    //Recommended to be all false,
    //since none of them does not work well
    //at this moment.
    total_energy_preservation=0;
    energy_preserve_per_vertex=0;
    energy_preserve_with_max_velocity=false;
    
    local_energy_preservation=0;

    //////////////////////////
    //Simulation parameters
    //dynamically used for a time step.
    //////////////////////////
    move_right=false;
    move_left=false;
    blowing_bubble0=false;
    absorbing_bubble0=false;
    
    bursting=false;
    
    do_stepScenes=false;

}

HGF::HGF(const std::vector<LosTopos::Vec3d> & vs, const std::vector<LosTopos::Vec3st> & fs, const std::vector<LosTopos::Vec2i> & ls, const std::vector<size_t> & constrained_vertices,  const std::vector<Vec3d> & constrained_positions, const int _num_bubbles,bool detect_boundary):num_bubbles(_num_bubbles)
{

    // construct the surface tracker
    double mean_edge_len = Options::doubleValue("remeshing-resolution");
    if (mean_edge_len == 0)
    {
        for (size_t i = 0; i < fs.size(); i++)
        {
            mean_edge_len += mag(vs[fs[i][0]] - vs[fs[i][1]]);
            mean_edge_len += mag(vs[fs[i][1]] - vs[fs[i][2]]);
            mean_edge_len += mag(vs[fs[i][2]] - vs[fs[i][0]]);
        }
        mean_edge_len /= (fs.size() * 3);
    }
    double min_edge_len = mean_edge_len * 0.5;
    double max_edge_len = mean_edge_len * 1.5;
    
    std::cout << "mean edge length = " << mean_edge_len << " min edge length = " << min_edge_len << " max edge length = " << max_edge_len << std::endl;
    
    LosTopos::SurfTrackInitializationParameters params;
    params.m_proximity_epsilon = Options::doubleValue("lostopos-collision-epsilon-fraction") * mean_edge_len;
    params.m_merge_proximity_epsilon = Options::doubleValue("lostopos-merge-proximity-epsilon-fraction") * mean_edge_len;
    params.m_allow_vertex_movement_during_collapse = true;
    params.m_perform_smoothing = Options::boolValue("lostopos-perform-smoothing");
    params.m_min_edge_length = min_edge_len;
    params.m_max_edge_length = max_edge_len;
    params.m_max_volume_change = Options::doubleValue("lostopos-max-volume-change-fraction") * pow(mean_edge_len, 3);
    params.m_min_triangle_angle = Options::doubleValue("lostopos-min-triangle-angle");
    params.m_max_triangle_angle = Options::doubleValue("lostopos-max-triangle-angle");
    params.m_large_triangle_angle_to_split = Options::doubleValue("lostopos-large-triangle-angle-to-split");
    params.m_min_triangle_area = Options::doubleValue("lostopos-min-triangle-area-fraction") * pow(mean_edge_len, 2);
    params.m_verbose = false;
    params.m_allow_non_manifold = Options::boolValue("lostopos-allow-non-manifold");
    params.m_allow_topology_changes = Options::boolValue("lostopos-allow-topology-changes");
    params.m_collision_safety = true;
    params.m_remesh_boundaries = true;
    params.m_t1_transition_enabled = Options::boolValue("lostopos-t1-transition-enabled");
    params.m_pull_apart_distance = Options::doubleValue("lostopos-t1-pull-apart-distance-fraction") * mean_edge_len;
    
    params.m_velocity_field_callback = NULL;
    
    if (Options::boolValue("lostopos-smooth-subdivision"))
        params.m_subdivision_scheme = new LosTopos::ModifiedButterflyScheme();
    else
        params.m_subdivision_scheme = new LosTopos::MidpointScheme();
    
    params.m_use_curvature_when_collapsing = false;
    params.m_use_curvature_when_splitting = false;

    if(!detect_boundary){
        m_constrained_vertices = constrained_vertices;
        m_constrained_positions = constrained_positions;
        
        std::vector<LosTopos::Vec3d> masses(vs.size(), LosTopos::Vec3d(1, 1, 1));
        for (size_t i = 0; i < m_constrained_vertices.size(); i++)
            masses[m_constrained_vertices[i]] *= std::numeric_limits<double>::infinity();
        m_st = new LosTopos::SurfTrack(vs, fs, ls, masses, params);
        
    }else{
        std::vector<LosTopos::Vec3d> masses(vs.size(), LosTopos::Vec3d(1, 1, 1));
        
        m_st = new LosTopos::SurfTrack(vs, fs, ls, masses, params);
        
        m_constrained_vertices.clear();
        m_constrained_positions.clear();
        detect_boundary_vertices();

        for (size_t i = 0; i < m_constrained_vertices.size(); i++)
            m_st->m_masses[m_constrained_vertices[i]] *= std::numeric_limits<double>::infinity();
    }
    m_st->m_solid_vertices_callback = this;
    m_st->m_mesheventcallback = this;

    // find out the number of regions
    num_region = 0;
    for (size_t i = 0; i < mesh().nt(); i++)
    {
        LosTopos::Vec2i l = mesh().get_triangle_label(i);
        num_region = std::max(num_region, l[0] + 1);
        num_region = std::max(num_region, l[1] + 1);
    }

    // initialize the dynamicc quantities i.e. velocity.
    m_v = new LosTopos::NonDestructiveTriMesh::VertexData<Vec3d>(&(m_st->m_mesh));
    for (size_t i = 0; i < mesh().nv(); i++){
        (*m_v)[i] =Vec3d(0,0,0);
    }

    //compute initial volumes
    Eigen::VectorXd volumes;
    Eigen::MatrixXd area_matrix;

    if(num_bubbles==-1){
        num_bubbles=num_region-1;
        AIR=0;
    }
    
    region_offset=num_region-num_bubbles;
    
    AIR=region_offset-1;
    set_parameters();
    
    //Necessary for volume correction.
    easy_orientation();

    if(num_bubbles>0){
        volumes_and_areas( volumes, area_matrix);
    }

    initialVolumes=volumes;
    if(logging_geometry){
        using std::cout;using std::endl;
        cout<<"initial volumes:";
        for(int bi=0;bi<num_bubbles;++bi){
            cout<<initialVolumes[bi]<<" ";
        }
        cout<<endl;
        
        cout<<"initial areas:";
        for(int bi=0;bi<num_bubbles;++bi){
            cout<<area_matrix(bi,bi)<<" ";
        }
        cout<<endl;
    }
}

HGF::~HGF()
{
    if (m_st)
        delete m_st;
}

void HGF::smoothVelocity_neighborhood(){
    
    std::vector<Vec3d> newVelocity = vel();
    for (size_t i = 0; i < mesh().nv(); i++)
    {

        Vec3d neighborhood_mean_velocity(0,0,0);
        
        int neighborhood_counter = 0;
        for (size_t l = 0; l < mesh().m_vertex_to_edge_map[i].size(); l++)
        {
            LosTopos::Vec2st e = mesh().m_edges[mesh().m_vertex_to_edge_map[i][l]];
            
            size_t vother = (e[0] == i ? e[1] : e[0]);

            neighborhood_mean_velocity += vel(vother);
            neighborhood_counter++;
            
        }
        if (neighborhood_counter != 0){
            neighborhood_mean_velocity/=neighborhood_counter;
            
        }
        newVelocity[i]=vel(i)+(neighborhood_mean_velocity-vel(i))* smooth;

    }
    
    for (size_t i = 0; i < mesh().nv(); i++)
    {
        vel(i)=newVelocity[i];
    }
    
}

double HGF::step(double dt){
    
    actual_dt=dt;
    
    static int count=0;
    ++count;
    static std::chrono::system_clock::time_point  start, mid1, mid2,end;
    
    auto nv=mesh().nv();

    {
        
        double old_energy;
        if(total_energy_preservation){
            
            old_energy=total_energy();
        }

//        nv=mesh().nv();
//        for (size_t i = 0; i < nv; i++){
//            m_st->pm_velocities[i] = vc(vel(i));
//        }
        
        auto los_topos_part=[this](){
            // Mesh improvement
            for(int i = 0; i < Options::intValue("remeshing-iterations"); i++)
            {
                
                m_st->topology_changes();
                m_st->improve_mesh();
                
                //Necessary for volume correction.
                easy_orientation();
                
            }
        };

        if(logging_detailed_time and count>1){
            start = std::chrono::system_clock::now();
            
            los_topos_part();
            
            mid1 = std::chrono::system_clock::now();

        }else{
            los_topos_part();
        }
        
        if(total_energy_preservation){

            double new_energy=total_energy();
            
            double correction_scale= std::pow( old_energy/new_energy,0.5);
            
            nv=mesh().nv();
            for(int vi=0;vi<nv;++vi){
                vel(vi)*=correction_scale;
            }

        }

        // defrag the mesh in the end, to ensure the next step starts with a clean mesh
        m_st->defrag_mesh_from_scratch(m_constrained_vertices);

    }

    // damping by smoothing
    if( smooth_from_neighbors)
    {
        smoothVelocity_neighborhood();
    }
    
    //Enforce constrained vertices to be fixed.
    if(1)
    {
        // before enforcing constraints, first scan through the mesh to find any solid vertices not registered as constraints. they can appear due to remeshing (splitting an all-solid edge)
        std::vector<int> constrained_vertices_map(mesh().nv(), -1);
        for (size_t i = 0, n_consv=m_constrained_vertices.size(); i < n_consv; i++)
            constrained_vertices_map[m_constrained_vertices[i]] = i;
        
        for (size_t i = 0,nv=mesh().nv(); i <nv; i++)
        {
            if (m_st->vertex_is_any_solid(i) && constrained_vertices_map[i] < 0)
            {
                m_constrained_vertices.push_back(i);
                m_constrained_positions.push_back(pos(i));
            }
        }

        //    // enforce the constraints exactly.
        for (size_t ii = 0,n_consv=m_constrained_vertices.size(); ii <n_consv ; ii++)
        {
            size_t i = m_constrained_vertices[ii];
            m_st->pm_newpositions[i] = vc(m_constrained_positions[ii]);
            vel(i)<<0,0,0;
        }
        
    }

    stepHGF(dt);

    // move the mesh
    auto integrate_part=[this,dt](){
        double actual_dt;
        m_st->integrate(dt, actual_dt);
        if (actual_dt != dt)
            std::cout << "Warning: SurfTrack::integrate() failed to step the full length of the time step!" << std::endl;

    };
    
    if(logging_detailed_time and count>1){
        
        mid2=std::chrono::system_clock::now();
        integrate_part();
        end = std::chrono::system_clock::now();

        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>((end-mid2)+(mid1-start)).count(); //Convert time to ms.

        time_los_topos+=elapsed;
        
    }else{
        integrate_part();
    }

    return actual_dt;
    
}

//
//Below are the functions for LosTopos, including the interpolation of the velocity for a newly genereted vertex.
//

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Callbacks
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool HGF::generate_collapsed_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos)
{
    if (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1))
    {
        return false;
    } else if (st.vertex_is_any_solid(v0))
    {
        pos = st.pm_positions[v0];
        return true;
    } else if (st.vertex_is_any_solid(v1))
    {
        pos = st.pm_positions[v1];
        return true;
    } else
    {
        pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
        return true;
    }
}

bool HGF::generate_split_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos)
{
    std::cout << "solid callback: generate split position: " << v0 << " " << v1 << " " << (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1)) << std::endl;
    pos = (st.pm_positions[v0] + st.pm_positions[v1]) / 2;
    if (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1))
        return false;
    else
        return true;
}

LosTopos::Vec3c HGF::generate_collapsed_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1)
{
    return LosTopos::Vec3c(1, 1, 1);
}

LosTopos::Vec3c HGF::generate_split_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1)
{
    if (st.vertex_is_any_solid(v0) && st.vertex_is_any_solid(v1))
        return LosTopos::Vec3c(1, 1, 1);
    else if (st.vertex_is_any_solid(v0))
        return LosTopos::Vec3c(0, 0, 0);
    else if (st.vertex_is_any_solid(v1))
        return LosTopos::Vec3c(0, 0, 0);
    else
        return LosTopos::Vec3c(0, 0, 0);
}

bool HGF::generate_edge_popped_positions(LosTopos::SurfTrack & st, size_t oldv, const LosTopos::Vec2i & cut, LosTopos::Vec3d & pos_upper, LosTopos::Vec3d & pos_lower)
{
    return false;
}

bool HGF::generate_vertex_popped_positions(LosTopos::SurfTrack & st, size_t oldv, int A, int B, LosTopos::Vec3d & pos_a, LosTopos::Vec3d & pos_b)
{
    return false;
}

bool HGF::solid_edge_is_feature(const LosTopos::SurfTrack & st, size_t e)
{
    return false;
}

LosTopos::Vec3d HGF::sampleVelocity(LosTopos::Vec3d & pos)
{
    return LosTopos::Vec3d(0, 0, 0);
}

bool HGF::sampleDirectionalDivergence(const LosTopos::Vec3d & pos, const LosTopos::Vec3d & dir, double & output)
{
    return false;
}

struct CollapseTempData
{
    size_t v0;
    size_t v1;
    
    Vec3d old_x0;
    Vec3d old_x1;
    
    Vec3d old_u0;
    Vec3d old_u1;
    
    //For local energy preservation
    double old_local_energy;
    double local_area;
};

void HGF::pre_collapse(const LosTopos::SurfTrack & st, size_t e, void ** data)
{
    CollapseTempData * td = new CollapseTempData;
    td->v0 = st.m_mesh.m_edges[e][0];
    td->v1 = st.m_mesh.m_edges[e][1];
    
    td->old_x0 = vc(st.pm_positions[td->v0]);
    td->old_x1 = vc(st.pm_positions[td->v1]);
    
    td->old_u0 = vel(td->v0);
    td->old_u1 = vel(td->v1);

    *data = (void *)td;
    std::cout << "pre collapse: " << e << ": " << td->v0 << " " << td->v1 << std::endl;

    if(local_energy_preservation){
        
        auto target_faces0=st.m_mesh.m_vertex_to_triangle_map [td->v0];
        auto target_faces1=st.m_mesh.m_vertex_to_triangle_map [td->v1];
        
        std::vector<size_t> target_faces=target_faces0;
        
        for(auto face1:target_faces1){
            bool already_in=std::find(target_faces.begin(),target_faces.end(),face1)!=target_faces.end();
            
            if(!already_in){
                target_faces.push_back(face1);
            }
        }

        double local_energy=0;
        
        double local_area=0;
        
        for(auto face:target_faces){
            LosTopos::Vec3st tri=st.m_mesh.m_tris[face];
            
            double tri_area=face_area(face);
            local_area+=tri_area;
            
            for(size_t vi=0;vi<3;++vi){
                size_t v=tri[vi];
                Vec3d velocity=vel(v);
                local_energy+=1.0/3.0*tri_area*velocity.squaredNorm();
            }
        }
        
        td->old_local_energy=local_energy;
        
        td->local_area=local_area;

    }
}

void HGF::post_collapse(const LosTopos::SurfTrack & st, size_t e, size_t merged_vertex, void * data)
{
    CollapseTempData * td = (CollapseTempData *)data;
    std::cout << "post collapse: " << e << ": " << td->v0 << " " << td->v1 << " => " << merged_vertex << std::endl;
    assert((st.m_mesh.vertex_is_deleted(td->v0) && merged_vertex == td->v1) || (st.m_mesh.vertex_is_deleted(td->v1) && merged_vertex == td->v0));
    
    Vec3d merged_x = vc(st.pm_positions[merged_vertex]);
    double s = (merged_x - td->old_x0).dot(td->old_x1 - td->old_x0) / (td->old_x1 - td->old_x0).squaredNorm();

    if(true){
        if (s > 1) s = 1;
        if (s < 0) s = 0;
        
        Vec3d new_u= td->old_u0 * (1 - s) + td->old_u1 * s;

        if(energy_preserve_per_vertex){
            double new_norm=(1-s)*td->old_u0.norm()+s*td->old_u1.norm();
            
            if(energy_preserve_with_max_velocity){
                new_norm=std::max(td->old_u0.norm(),td->old_u1.norm());
            }
            new_u=new_norm*new_u.normalized();
            
        }
        
        vel(merged_vertex) = new_u;
    }
    
    if(local_energy_preservation){
        auto target_faces=st.m_mesh.m_vertex_to_triangle_map [merged_vertex];

        double new_local_energy=0;
        double vertex_area=0;
        
        double local_area=0;
        
        for(auto face:target_faces){
            LosTopos::Vec3st tri=st.m_mesh.m_tris[face];
            
            double tri_area=face_area(face);
            vertex_area+=1.0/3.0*tri_area;
            
            local_area+=tri_area;
            
            for(size_t vi=0;vi<3;++vi){
                size_t v=tri[vi];
                
                if(v==merged_vertex){continue;}
                Vec3d velocity=vel(v);
                new_local_energy+=1.0/3.0*tri_area*velocity.squaredNorm();
            }
        }
        
        double new_squared_norm=(td->old_local_energy-new_local_energy)/vertex_area;

        vel(merged_vertex)=vel(merged_vertex).normalized()*sqrt(new_squared_norm);

    }

    //Smmooth the velocity of the new vertex.
    if(false){
        
        Vec3d newVelocity(0,0,0);
        
        newVelocity=Vec3d(0,0,0);
        
        Vec3d neighborhood_mean_velocity(0,0,0);
        
        int neighborhood_counter = 0;
        for (size_t k = 0; k < st.m_mesh.m_vertex_to_edge_map[merged_vertex].size(); k++)
        {
            LosTopos::Vec2st e = st.m_mesh.m_edges[st.m_mesh.m_vertex_to_edge_map[merged_vertex][k]];
            size_t vother = (e[0] == merged_vertex ? e[1] : e[0]);

            neighborhood_mean_velocity += vel(vother);
            neighborhood_counter++;
            
        }
        if (neighborhood_counter != 0){
            neighborhood_mean_velocity/=neighborhood_counter;
            
        }
        
        double smoothing_coef=0.5;//0 ~ 1.
        
        newVelocity=vel(merged_vertex)+(neighborhood_mean_velocity-vel(merged_vertex))* smoothing_coef;

        vel(merged_vertex) = newVelocity;
        
    }
    
}

struct SplitTempData
{
    size_t v0;
    size_t v1;
    
    Vec3d old_x0;
    Vec3d old_x1;
    
    Vec3d old_u0;
    Vec3d old_u1;

    //For local energy preservation
    double old_local_energy;
    double local_area;
};

void HGF::pre_split(const LosTopos::SurfTrack & st, size_t e, void ** data)
{
    SplitTempData * td = new SplitTempData;
    td->v0 = st.m_mesh.m_edges[e][0];
    td->v1 = st.m_mesh.m_edges[e][1];
    
    td->old_x0 = vc(st.pm_positions[td->v0]);
    td->old_x1 = vc(st.pm_positions[td->v1]);
    
    td->old_u0 = vel(td->v0);
    td->old_u1 = vel(td->v1);

    *data = (void *)td;
    std::cout << "pre split: " << e << ": " << td->v0 << " " << td->v1 << std::endl;
    
    if(local_energy_preservation){
        auto edges=st.m_mesh.m_vertex_to_edge_map[td->v0];
        size_t common_edge=-1;
        for(auto edge:edges){
            size_t vother=edge_other_vertex(edge, td->v0);
            if(vother==td->v1){
                common_edge=edge;
                break;
            }
        }
        assert(common_edge!=-1);
        
        auto target_faces=st.m_mesh.m_edge_to_triangle_map[common_edge];
        
        assert(target_faces.size()==2 );
        
        double local_energy=0;
        double local_area=0;
        
        for(auto face:target_faces){
            LosTopos::Vec3st tri=st.m_mesh.m_tris[face];
            
            double tri_area=face_area(face);
            local_area+=tri_area;
            
            for(size_t vi=0;vi<3;++vi){
                size_t v=tri[vi];
                Vec3d velocity=vel(v);
                local_energy+=1.0/3.0*tri_area*velocity.squaredNorm();
            }
        }
        
        td->old_local_energy=local_energy;
        
        td->local_area=local_area;
        
    }

}

void HGF::post_split(const LosTopos::SurfTrack & st, size_t e, size_t new_vertex, void * data)
{
    SplitTempData * td = (SplitTempData *)data;
    std::cout << "post split: " << e << ": " << td->v0 << " " << td->v1 << " => " << new_vertex << std::endl;
    
    Vec3d midpoint_x = vc(st.pm_positions[new_vertex]);
    double s = (midpoint_x - td->old_x0).dot(td->old_x1 - td->old_x0) / (td->old_x1 - td->old_x0).squaredNorm();

    if(true){
        if (s > 1) s = 1;
        if (s < 0) s = 0;
        
        Vec3d new_u= td->old_u0 * (1 - s) + td->old_u1 * s;

        if(energy_preserve_per_vertex){
            double new_norm=(1-s)*td->old_u0.norm()+s*td->old_u1.norm();
            
            if(energy_preserve_with_max_velocity){
                new_norm=std::max(td->old_u0.norm(),td->old_u1.norm());
            }
            new_u=new_norm*new_u.normalized();

        }
        
        vel(new_vertex) = new_u;
    }

    if(local_energy_preservation){
        auto target_faces=st.m_mesh.m_vertex_to_triangle_map [new_vertex];

        assert(target_faces.size()==4);
        double new_local_energy=0;
        double vertex_area=0;
        double local_area=0;
        
        for(auto face:target_faces){
            LosTopos::Vec3st tri=st.m_mesh.m_tris[face];
            
            double tri_area=face_area(face);
            vertex_area+=1.0/3.0*tri_area;
            
            local_area+=tri_area;
            
            for(size_t vi=0;vi<3;++vi){
                size_t v=tri[vi];
                
                if(v==new_vertex){continue;}
                Vec3d velocity=vel(v);
                new_local_energy+=1.0/3.0*tri_area*velocity.squaredNorm();
            }
        }
        
        auto v_area=vert_area(new_vertex);
        double new_squared_norm=(td->old_local_energy-new_local_energy)/vertex_area;

        vel(new_vertex)=vel(new_vertex).normalized()*sqrt(new_squared_norm);

    }
    
    //Smmooth the velocity o the new vertex.
    if(false){
        
        Vec3d newVelocity(0,0,0);
        
        newVelocity=Vec3d(0,0,0);
        
        Vec3d neighborhood_mean_velocity(0,0,0);
        
        int neighborhood_counter = 0;
        for (size_t k = 0; k < st.m_mesh.m_vertex_to_edge_map[new_vertex].size(); k++)
        {
            LosTopos::Vec2st e = st.m_mesh.m_edges[st.m_mesh.m_vertex_to_edge_map[new_vertex][k]];
            size_t vother = (e[0] == new_vertex ? e[1] : e[0]);

            neighborhood_mean_velocity += vel(vother);
            neighborhood_counter++;
            
        }
        if (neighborhood_counter != 0){
            neighborhood_mean_velocity/=neighborhood_counter;
            
        }
        
        double smoothing_coef=0.5;//0 ~ 1.
        
        newVelocity=vel(new_vertex)+(neighborhood_mean_velocity-vel(new_vertex))* smoothing_coef;

        vel(new_vertex) = newVelocity;
        
    }
    
}

void HGF::pre_flip(const LosTopos::SurfTrack & st, size_t e, void ** data)
{
    
}

void HGF::post_flip(const LosTopos::SurfTrack & st, size_t e, void * data)
{
    
}

struct T1TempData
{
    
};

void HGF::pre_t1(const LosTopos::SurfTrack & st, size_t v, void ** data)
{
    T1TempData * td = new T1TempData;
    
    *data = (void *)td;
}

void HGF::post_t1(const LosTopos::SurfTrack & st, size_t v, size_t a, size_t b, void * data)
{
    std::cout << "v = " << v << " -> " << a << " " << b << std::endl;
    T1TempData * td = (T1TempData *)data;

    if(true){
        vel(a) = vel(v);
        vel(b) = vel(v);
    }
}

struct FaceSplitTempData
{
    size_t v0;
    size_t v1;
    size_t v2;
    
    Vec3d old_x0;
    Vec3d old_x1;
    Vec3d old_x2;
    
    Vec3d old_u0;
    Vec3d old_u1;
    Vec3d old_u2;
};

void HGF::pre_facesplit(const LosTopos::SurfTrack & st, size_t f, void ** data)
{
    FaceSplitTempData * td = new FaceSplitTempData;
    
    td->v0 = st.m_mesh.get_triangle(f)[0];
    td->v1 = st.m_mesh.get_triangle(f)[1];
    td->v2 = st.m_mesh.get_triangle(f)[2];

    td->old_x0 = vc(st.pm_positions[td->v0]);
    td->old_x1 = vc(st.pm_positions[td->v1]);
    td->old_x2 = vc(st.pm_positions[td->v2]);
    
    td->old_u0 = vel(td->v0);
    td->old_u1 = vel(td->v1);
    td->old_u2 = vel(td->v2);
    
    *data = (void *)td;
}

void HGF::post_facesplit(const LosTopos::SurfTrack & st, size_t f, size_t new_vertex, void * data)
{
    FaceSplitTempData * td = (FaceSplitTempData *)data;

    if(true){
        Vec3d new_x = vc(st.pm_positions[new_vertex]);
        Vec3d c = Vec3d::Zero();
        Vec3d n = (td->old_x1 - td->old_x0).cross(td->old_x2 - td->old_x0);
        double nsq = n.squaredNorm();
        c[0] = 1 - (new_x - td->old_x0).dot(n.cross(td->old_x1 - td->old_x2)) / nsq;
        c[1] = 1 - (new_x - td->old_x1).dot(n.cross(td->old_x2 - td->old_x0)) / nsq;
        if (c[0] > 1)        c[0] = 1;
        if (c[0] < 0)        c[0] = 0;
        if (c[1] > 1 - c[0]) c[1] = 1 - c[0];
        if (c[1] < 0)        c[1] = 0;
        c[2] = 1 - c[0] - c[1];
        
        vel(new_vertex) = td->old_u0 * c[0] + td->old_u1 * c[1] + td->old_u2 * c[2];
        
        if(energy_preserve_per_vertex){
            double new_norm=c[0]*td->old_u0.norm()+c[1]*td->old_u1.norm()*c[2]*td->old_u2.norm();
            
            if(energy_preserve_with_max_velocity){
                std::vector<double> norms={td->old_u0.norm(),td->old_u1.norm(),td->old_u2.norm()};
                new_norm=*std::max_element(norms.begin(),norms.end());
            }
            
            vel(new_vertex)=new_norm*vel(new_vertex).normalized();
        }
        
    }
}

struct SnapTempData
{
    size_t v0;
    size_t v1;
    
    Vec3d old_x0;
    Vec3d old_x1;

    Vec3d old_u0;
    Vec3d old_u1;
};

void HGF::pre_snap(const LosTopos::SurfTrack & st, size_t v0, size_t v1, void ** data)
{
    SnapTempData * td = new SnapTempData;
    td->v0 = v0;
    td->v1 = v1;
    
    td->old_x0 = vc(st.pm_positions[v0]);
    td->old_x1 = vc(st.pm_positions[v1]);

    td->old_u0 = vel(td->v0);
    td->old_u1 = vel(td->v1);

    *data = (void *)td;
    std::cout << "pre snap: " << v0 << " " << v1 << std::endl;
}

void HGF::post_snap(const LosTopos::SurfTrack & st, size_t v_kept, size_t v_deleted, void * data)
{
    SnapTempData * td = (SnapTempData *)data;
    std::cout << "post snap: " << td->v0 << " " << td->v1 << " => " << v_kept << std::endl;
    assert((td->v0 == v_kept && td->v1 == v_deleted) || (td->v1 == v_kept && td->v0 == v_deleted));
    assert(v_kept != v_deleted);
    assert(st.m_mesh.vertex_is_deleted(v_deleted));
    assert(!st.m_mesh.vertex_is_deleted(v_kept));
    
    Vec3d merged_x = vc(st.pm_positions[v_kept]);
    double s = (merged_x - td->old_x0).dot(td->old_x1 - td->old_x0) / (td->old_x1 - td->old_x0).squaredNorm();

    if(true){
        
        if (s > 1) s = 1;
        if (s < 0) s = 0;
        Vec3d new_u = td->old_u0 * (1 - s) + td->old_u1 * s;

        if(energy_preserve_per_vertex){
            double new_norm=(1-s)*td->old_u0.norm()+s*td->old_u1.norm();

            if(energy_preserve_with_max_velocity){
                new_norm=std::max(td->old_u0.norm(),td->old_u1.norm());
                
            }
            new_u=new_norm*new_u.normalized();
            
        }
        
        vel(v_kept) = new_u;
    }
    
    //Smmooth the velocity of the new vertex.
    if(false){
        Vec3d newVelocity(0,0,0);
        newVelocity=Vec3d(0,0,0);
        
        Vec3d neighborhood_mean_velocity(0,0,0);
        
        int neighborhood_counter = 0;
        for (size_t k = 0; k < st.m_mesh.m_vertex_to_edge_map[v_kept].size(); k++)
        {
            LosTopos::Vec2st e = st.m_mesh.m_edges[st.m_mesh.m_vertex_to_edge_map[v_kept][k]];
            size_t vother = (e[0] == v_kept ? e[1] : e[0]);

            neighborhood_mean_velocity += vel(vother);
            neighborhood_counter++;
            
        }
        if (neighborhood_counter != 0){
            neighborhood_mean_velocity/=neighborhood_counter;
            
        }
        
        double smoothing_coef=0.5;//0 ~ 1.
        
        newVelocity=vel(v_kept)+(neighborhood_mean_velocity-vel(v_kept))* smoothing_coef;

        vel(v_kept) = newVelocity;
        
    }
    
}

void HGF::pre_smoothing(const LosTopos::SurfTrack & st, void ** data)
{
    
}

void HGF::post_smoothing(const LosTopos::SurfTrack & st, void * data)
{
    
}
