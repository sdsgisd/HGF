//
//  HGF_Utility.cpp
//
//  Created by Sadashige Ishida on 2/4/17.
//
//

#include <stdio.h>
#include <map>
#include <numeric>

#include <Eigen/Sparse>
#include <igl/per_vertex_normals.h>
#include <igl/writeOBJ.h>
#include <igl/read_triangle_mesh.h>

#include "HGF.h"
#include "MeshIO.h"

void HGF::writeObj_FaceLabel_constrainedVertices(const bool with_imaginary_vertices){
    
    const bool with_normal=false;
    const bool change_y_z=true;
    MeshIO::saveOBJ(*this, "mesh.obj",with_imaginary_vertices,with_normal,change_y_z);

    std::string labelfile="flabel.txt";
    std::ofstream ofs_flabel(labelfile);
    
    int nt=(int)mesh().nt();
    ofs_flabel<<num_region<<std::endl;
    ofs_flabel<<nt<<std::endl;
    for(int i=0;i<nt;++i){
        
        LosTopos::Vec2i label = mesh().get_triangle_label(i);
        ofs_flabel<<label[0]<<" "<<label[1]<<std::endl;
    }

    //Write constrained vertices
    int nv_const=m_constrained_vertices.size();
    
    std::string constvertexfile="constvertex.txt";
    std::ofstream ofs_constvertex(constvertexfile);
    ofs_constvertex<<nv_const<<std::endl;
    for(int i=0;i<nv_const;++i){
        
        ofs_constvertex<<m_constrained_vertices[i]<<std::endl;
    }

}

void HGF::write_constrained_mesh(std::string outputfile){
    
    if(Constrained_V.size()>0){
        
        Eigen::MatrixXd y_z_changed_Constrained_V(Constrained_V.rows(),3);
        y_z_changed_Constrained_V<<Constrained_V.col(0),Constrained_V.col(2),Constrained_V.col(1);
        
        igl::writeOBJ(outputfile,y_z_changed_Constrained_V,Constrained_F);
    }
}

void HGF::write_film_mesh(std::string outputfile){
    
    MeshIO::saveOBJ(*this, outputfile,false,false,true);
}

void HGF::write_mesh(){
    static int count=0;

    std::string output_directory="outputmesh/";
    std::stringstream output_meshfile;
    output_meshfile<<output_directory<<"output_mesh"<<std::setfill('0') << std::setw(6)<<count<<".obj";

    std::stringstream constrained_meshfile;
    constrained_meshfile<<output_directory<<"output_frame"<<std::setfill('0') << std::setw(6)<<count<<".obj";
    
    write_film_mesh(output_meshfile.str());
    write_constrained_mesh(constrained_meshfile.str());
    ++count;
}

void HGF::geometric_information(double time){
    
    // compute the volumes of Fang Da program.
    std::vector<double> volumes(num_region, 0);
    Vec3d xref(0, 0, 0);
    for (size_t i = 0; i < mesh().nt(); i++)
    {
        LosTopos::Vec3st t = mesh().get_triangle(i);
        Vec3d x0 = pos(t[0]);
        Vec3d x1 = pos(t[1]);
        Vec3d x2 = pos(t[2]);
        double v = (x0 - xref).cross(x1 - xref).dot(x2 - xref);
        
        LosTopos::Vec2i l = mesh().get_triangle_label(i);
        volumes[l[0]] += v;
        volumes[l[1]] -= v;
    }
    
    std::cout<<"region volumes:";
    for (size_t i = 0; i < volumes.size(); i++)
        std::cout << " " << volumes[i]/6.0;
    std::cout << std::endl;

    {//Volumes of only closed regions i.e. bubbles.
        Eigen::VectorXd volume_vector;
        Eigen::MatrixXd area_matrix;
        volumes_and_areas( volume_vector, area_matrix);
        
        std::cout<<"bubble volumes:";
        for (size_t i = 0; i < volume_vector.size(); i++)
            std::cerr << " " << volume_vector[i];
        std::cerr << std::endl;
        //Compute volumes till here.
    }

    static int max_nv=mesh().nv();
    static int max_nt=mesh().nt();
    
    int nv=mesh().nv();
    max_nv=nv>max_nv?nv:max_nv;
    
    int nt=mesh().nt();
    max_nt=nt>max_nt?nt:max_nt;
    
    std::cout<<"nv:"<<nv<<",max_nv:"<<max_nv<<",nt:"<<nt<<",max_nt:"<<max_nt<<std::endl;

    tq_junc_acurracy(time,write_geometry);
    total_area(time,write_geometry);

}

void HGF::tq_junc_acurracy(double time,bool write){
    
    //Collect t-junctions.
    std::vector<size_t> t_junc_edge;
    int ne=mesh().ne();
    for(size_t ei=0;ei<ne;++ei){
        
        bool has3faces=(mesh().m_edge_to_triangle_map[ei].size()==3);
        bool notConstrained= !m_st->edge_is_all_solid(ei);
        
        if(has3faces and notConstrained){
            t_junc_edge.push_back(ei);
        }
    }
    
    //Compute dihedral angles.
    int num_te=(int)t_junc_edge.size();
    std::vector<double> t_junc_degrees(num_te*3);
    
    for(size_t tei=0;tei<num_te;++tei){
        size_t te=t_junc_edge[tei];
        std::vector<size_t> incident_tris= mesh().m_edge_to_triangle_map[te];
        
        Vec3d normal0=face_outward_normal(incident_tris[0]);
        Vec3d normal1=face_outward_normal(incident_tris[1]);
        Vec3d normal2=face_outward_normal(incident_tris[2]);
        
        //Angles are computed to be not approx 60 deg but approx 120 deg.
        double norm01=-std::abs(normal0.dot(normal1));
        double norm12=-std::abs(normal1.dot(normal2));
        double norm20=-std::abs(normal2.dot(normal0));
        
        double deg01_rag= std::acos(norm01);
        double deg12_rag= std::acos(norm12);
        double deg20_rag= std::acos(norm20);
        
        double deg01=deg01_rag*180/M_PI;
        double deg12=deg12_rag*180/M_PI;
        double deg20=deg20_rag*180/M_PI;
        
        const double degree120=120;
        
        t_junc_degrees[3*tei]=deg01;
        t_junc_degrees[3*tei+1]=deg12;
        t_junc_degrees[3*tei+2]=deg20;
    }

    //Output deviations of angles from 120 deg.
    double average=std::accumulate(t_junc_degrees.begin(), t_junc_degrees.end(), 0.0)/t_junc_degrees.size();

    class sum_error_t {
    public:
        const double degree120=120;
        double operator() (const double before, const double value) {
            
            return before + std::abs(value-degree120);
            
        }
    };
    
    class square_error_t {
    public:
        const double degree120=120;
        double operator() (const double before, const double value) {
            
            double abs_error=std::abs(value-degree120);
            return before + abs_error*abs_error;
        }
    };
    
    class maximam_error_t {
    public:
        const double degree120=120;
        double operator() (const double before, const double value) {
            
            return std::max(before, std::abs(value-degree120));
        }
    };

    double average_error_t=std::accumulate(t_junc_degrees.begin(), t_junc_degrees.end(), 0.0,sum_error_t())/t_junc_degrees.size();
    
    double rms_error_t=std::sqrt(std::accumulate(t_junc_degrees.begin(), t_junc_degrees.end(), 0.0,square_error_t())/t_junc_degrees.size());
    
    double max_error_t=std::accumulate(t_junc_degrees.begin(), t_junc_degrees.end(), 0.0,maximam_error_t());
    std::cout<<"average_error_t:"<<average_error_t<<",rms:"<<rms_error_t<<",max_error_T:"<<max_error_t<<std::endl;

    //Collect quad junctions.
    int nv=mesh().nv();
    std::vector<bool> is_tq_junc(nv);
    for(int vi=0;vi<nv;++vi){
        is_tq_junc[vi]=mesh().is_vertex_nonmanifold(vi);
    }

    std::vector<size_t>q_junctions;
    for(size_t v=0;v<nv;++v){
        if(!is_tq_junc[v]){
            continue;
        }
        auto incident_edges=mesh().m_vertex_to_edge_map[v];
        int num_incident_edges=incident_edges.size();
        
        int t_junc_counter=0;
        for(size_t ei=0;ei<num_incident_edges;++ei){
            size_t edge=incident_edges[ei];
            if(mesh().is_edge_nonmanifold(edge)){
                ++t_junc_counter;
            }
        }

        if(t_junc_counter==4 and !m_st->vertex_is_all_solid(v)){
            q_junctions.push_back(v);

        }

    }
    
    //Compute angles of t-junction edges at q-junctions.
    int num_qj=q_junctions.size();
    
    std::vector<double> q_junc_degrees(num_qj*6);
    
    for(size_t qi=0;qi<num_qj;++qi){
        size_t q=q_junctions[qi];
        
        auto incident_edges=mesh().m_vertex_to_edge_map[q];
        int num_incident_edges=incident_edges.size();
        
        std::vector<size_t>four_vertices;
        
        int t_junc_counter=0;
        for(size_t ei=0;ei<num_incident_edges;++ei){
            //Check if edge is a t-junc edge.
            size_t edge=incident_edges[ei];
            if(mesh().is_edge_nonmanifold(edge) ){
                size_t vother=edge_other_vertex(edge, q);
                four_vertices.push_back(vother);
            }
        }
        
        int size_four=four_vertices.size();
        
        //Compute angles.
        Vec3d edge0=vc(m_st->pm_positions[four_vertices[0]]-m_st->pm_positions[q]).normalized();
        Vec3d edge1=vc(m_st->pm_positions[four_vertices[1]]-m_st->pm_positions[q]).normalized();
        Vec3d edge2=vc(m_st->pm_positions[four_vertices[2]]-m_st->pm_positions[q]).normalized();
        Vec3d edge3=vc(m_st->pm_positions[four_vertices[3]]-m_st->pm_positions[q]).normalized();
        //The combination to choose two edges amoung four has 4C2=6 patterns.
        
        double norm01=edge0.dot(edge1);
        double norm02=edge0.dot(edge2);
        double norm03=edge0.dot(edge3);
        double norm12=edge1.dot(edge2);
        double norm13=edge1.dot(edge3);
        double norm23=edge2.dot(edge3);
        
        double deg01_rag= std::acos(norm01);
        double deg02_rag= std::acos(norm02);
        double deg03_rag= std::acos(norm03);
        double deg12_rag= std::acos(norm12);
        double deg13_rag= std::acos(norm13);
        double deg23_rag= std::acos(norm23);
        
        double deg01=deg01_rag*180/M_PI;
        double deg02=deg02_rag*180/M_PI;
        double deg03=deg03_rag*180/M_PI;
        double deg12=deg12_rag*180/M_PI;
        double deg13=deg13_rag*180/M_PI;
        double deg23=deg23_rag*180/M_PI;
        
        const double degree109=std::acos(-1.0/3)*180/M_PI;

        q_junc_degrees[6*qi]=deg01;
        q_junc_degrees[6*qi+1]=deg02;
        q_junc_degrees[6*qi+2]=deg03;
        q_junc_degrees[6*qi+3]=deg12;
        q_junc_degrees[6*qi+4]=deg13;
        q_junc_degrees[6*qi+5]=deg23;

    }
    
    //Output deviations of angles from arc(-1/3).
    
    class sum_error_q {
    public:
        const double degree109=std::acos(-1.0/3)*180/M_PI;
        double operator() (const double before, const double value) {
            
            return before + std::abs(value-degree109);
        }
    };
    
    class square_error_q {
    public:
        const double degree109=std::acos(-1.0/3)*180/M_PI;
        double operator() (const double before, const double value) {

            double abs_error=std::abs(value-degree109);
            return before + abs_error*abs_error;
        }
    };
    
    class maximam_error_q{
    public:
        const double degree109=std::acos(-1.0/3)*180/M_PI;
        double operator() (const double before, const double value) {
            
            return std::max(before, std::abs(value-degree109));
        }
    };

    double average_error_q=std::accumulate(q_junc_degrees.begin(), q_junc_degrees.end(), 0.0,sum_error_q())/q_junc_degrees.size();
    
    double max_error_q=std::accumulate(q_junc_degrees.begin(), q_junc_degrees.end(), 0.0,maximam_error_q());
    
    double rms_error_q=std::sqrt(std::accumulate(q_junc_degrees.begin(), q_junc_degrees.end(), 0.0,square_error_q())/q_junc_degrees.size());
    
    std::cout<<"average_error_q:"<<average_error_q<<",rms:"<<rms_error_q<<",max_error_q:"<<max_error_q<<std::endl;

    if(write){
        static std::string directory="./";
        static std::string t_filename="t_junc.csv";
        static std::ofstream t_output(directory+t_filename);
        
        t_output<<time<<","<<rms_error_t<<std::endl;
        
        static std::string q_filename="q_junc.csv";
        static std::ofstream q_output(directory+q_filename);
        
        q_output<<time<<","<<rms_error_q<<std::endl;
        
        static int count=0;
        int max_count=-1;
        if(count++==max_count){
            t_output.close();
            q_output.close();
            
            exit(0);
            std::cout<<"exit in write_total_area";
        }
        
    }

}

void HGF::total_area(double time,bool write){
    
    using namespace Eigen;
    double area=0;
    
    int nt=(int)mesh().nt();
    for (size_t t = 0; t < nt; t++)
    {
        
        if(m_st->triangle_is_all_solid(t)){
            continue;
        }
        auto f=mesh().get_triangle(t);

        LosTopos::Vec3d p0 = m_st->pm_positions[f.v[0]];
        LosTopos::Vec3d  p1 = m_st->pm_positions[f.v[1]];
        LosTopos::Vec3d  p2 = m_st->pm_positions[f.v[2]];
        
        LosTopos::Vec3d n=LosTopos::cross(p1-p0, p2-p0);

        double local_area=LosTopos::mag(n);
        area+=local_area;

    }
    
    area/=2.0;
    
    std::cout<<"total_area:"<<area<<std::endl;

    static std::string directory="./";
    static std::string filename="area_transition.csv";
    
    if(write){
        
        static std::ofstream output(directory+filename);
        static int count=0;
        
        output<<time<<","<<area<<std::endl;
        
        int max_count=-1;
        if(count++==max_count){
            output.close();
            exit(0);
            std::cout<<"exit in write_total_area";
        }
        
    }
    
}

void HGF::vertex_areas(std::vector<double> & v_areas){
    
    int nv=mesh().nv();
    v_areas.resize(nv);
    
    int nt=(int)mesh().nt();
    for (size_t t = 0; t < nt; t++)
    {
        
        auto tri=mesh().get_triangle(t);
        
        LosTopos::Vec3d p0 = m_st->pm_positions[tri[0]];
        LosTopos::Vec3d  p1 = m_st->pm_positions[tri[1]];
        LosTopos::Vec3d  p2 = m_st->pm_positions[tri[2]];
        
        LosTopos::Vec3d n=LosTopos::cross(p1-p0, p2-p0);
        
        double tri_area=LosTopos::mag(n);
        
        for(size_t v=0;v<3;++v){
            v_areas[tri[v]]+=1/3.0*tri_area;
        }
        
    }

}

double HGF::total_energy(){
    double energy=0;
    
    int nv=mesh().nv();
    
    std::vector<double>v_areas;
    vertex_areas(v_areas);
    
    for(int vi=0;vi<nv;++vi){
        double local_voronoi_area=v_areas[vi];
        
        energy+=local_voronoi_area* vel(vi).squaredNorm();
    }
    
    return energy;
}

void HGF::easy_orientation(){
    //Set each triangle label (i,j) to be i>j.
    
    int nt=mesh().nt();
    for(int ti=0;ti<nt;++ti){
        int &label0=mesh().m_triangle_labels[ti][0];
        int &label1=mesh().m_triangle_labels[ti][1];
        
        if(label0<label1){
            std::swap( mesh().m_tris[ti][0],mesh().m_tris[ti][1]);
            std::swap( label0,label1);
        }
    }
}

void HGF::detect_boundary_vertices(){

    int ne=mesh().ne();
    for(int ei=0;ei<ne;++ei){
        
    }
    
    int nv=mesh().nv();
    for(int vi=0;vi<nv;++vi){
        
        auto adjacent_adges=mesh().m_vertex_to_edge_map[vi];
        for(size_t e:adjacent_adges){

            if(mesh().m_edge_to_triangle_map[e].size()==1){
                m_constrained_vertices.push_back(vi);
                m_constrained_positions.push_back(vc(m_st->pm_positions[vi]));
                continue;
            }
        }
    }

}

void HGF::face_edge_vertex_style(Eigen::MatrixXi &FE,Eigen::MatrixXi &EV,Eigen::MatrixXd &V){
    //Make data structure where each row i of FE indicates the edges of the i-th triangle, and each row j of EV indicates the vertices of the j-th edge.
    
    using namespace Eigen;
    
    size_t nv=m_st->m_mesh.nv();
    size_t nt=m_st->m_mesh.nt();
    V=MatrixXd::Zero(nv, 3);
    Eigen::MatrixXi F=MatrixXi::Zero(nt, 3);

    //Convert vertices and faces to U,F, duDt.
    for (size_t i = 0; i < nv; i++)
    {
        V.row(i)<<m_st->pm_positions[i][0],m_st->pm_positions[i][1],m_st->pm_positions[i][2];

    }
    for (size_t i = 0; i < nt; i++)
    {
        F.row(i)<<m_st->m_mesh.m_tris[i][0],m_st->m_mesh.m_tris[i][1],m_st->m_mesh.m_tris[i][2];
    }

    class edge{
    public:
        int small_ind_ver;
        int large_ind_ver;
        edge(int i0,int i1)
        :small_ind_ver(i0<i1?i0:i1),large_ind_ver(i0>i1?i0:i1)
        {
            assert(small_ind_ver<large_ind_ver);
        }
        bool operator ==(const edge & e)const{
            return small_ind_ver==e.small_ind_ver and\
            large_ind_ver==e.large_ind_ver;
        }
        
        //For the comparing operation for sorting by C++ standard library.
        bool operator<(const edge &e) const{
            return small_ind_ver<e.small_ind_ver or (small_ind_ver==e.small_ind_ver and large_ind_ver<e.large_ind_ver);
        }
    };
    
    std::vector<edge> edges;
    std::vector<std::unique_ptr< std::vector<size_t>>> fe_vec(nt);

    for(int fi=0;fi<nt;++fi){
        int v0=F(fi,0);
        int v1=F(fi,1);
        int v2=F(fi,2);
        
        edge e=edge(v0,v1);
        auto eit=std::find(edges.begin(), edges.end(), e);
        if(eit==edges.end()){
            //True if not found.
            edges.push_back(e);
            fe_vec[fi]->push_back(edges.size()-1);
        }else{
            fe_vec[fi]->push_back(std::distance(edges.begin(), eit));
            
        }
        
        e=edge(v1,v2);
        eit=std::find(edges.begin(), edges.end(), e);
        if(eit==edges.end()){
            edges.push_back(e);
            fe_vec[fi]->push_back(edges.size()-1);
        }else{
            fe_vec[fi]->push_back(std::distance(edges.begin(), eit));
            
        }
        
        e=edge(v2,v0);
        eit=std::find(edges.begin(), edges.end(), e);
        if(eit==edges.end()){
            edges.push_back(e);
            fe_vec[fi]->push_back(edges.size()-1);
        }else{
            fe_vec[fi]->push_back(std::distance(edges.begin(), eit));
            
        }
        
    }
    
    FE=MatrixXi(3,nt);
    for (size_t fi = 0; fi < nt; fi++)
    {
        FE.row(fi)<<fe_vec[fi]->at(0),fe_vec[fi]->at(1),fe_vec[fi]->at(2);
    }
    
    int ne=edges.size();
    EV=MatrixXi(2,ne);
    for (size_t ei = 0; ei < ne; ei++)
    {
        EV.row(ei)<<edges[ei].large_ind_ver, edges[ei].small_ind_ver;
    }
}

void HGF::volumes_and_areas_sparse(const Eigen::MatrixXd &targetU , Eigen::VectorXd &volumes, Eigen::SparseMatrix<double>& area_matrix ){
    
    using namespace Eigen;
    
    volumes=VectorXd::Zero(num_bubbles);
    area_matrix.resize(num_bubbles,num_bubbles);
    
    std::vector<Triplet<double> > area_triplet;
    
    int nt=(int)mesh().nt();
    for (size_t t = 0; t < nt; t++)
    {

        auto f=mesh().get_triangle(t);
        Eigen::Vector3d p0 = targetU.row(f.v[0]);
        Eigen::Vector3d p1 = targetU.row(f.v[1]);;
        Eigen::Vector3d p2 = targetU.row(f.v[2]);

        Eigen::Vector3d n = (p1-p0).cross(p2-p0);
        
        int label0=mesh().get_triangle_label(t)[0];
        int label1=mesh().get_triangle_label(t)[1];
        
        double local_area=n.norm();

        if(label0<=AIR and label1>AIR){
            area_triplet.emplace_back(label1-region_offset,label1-region_offset,local_area);

        }else if(label1<=AIR and label0>AIR){
            area_triplet.emplace_back(label0-region_offset,label0-region_offset,local_area);
        }else if(label0>AIR and label1>AIR){
            area_triplet.emplace_back(label0-region_offset,label1-region_offset,local_area);
            area_triplet.emplace_back(label1-region_offset,label0-region_offset,local_area);
  
        }
        
        area_matrix.setFromTriplets(area_triplet.begin(), area_triplet.end());
        
        // total volume via divergence theorem: ∫ 1
        double local_volume=n.dot(p0);

        if(label0>AIR){
            volumes[label0-region_offset]+=local_volume;
        }
        if(label1>AIR){
            volumes[label1-region_offset]-=local_volume;
            
        }
    }

    volumes/=6.0;
    area_matrix/=2.0;

    for(int i=0;i<num_bubbles;++i){
        area_matrix.coeffRef(i, 0)+=-1;
        area_matrix.coeffRef(i, 1)*=-1;
        area_matrix.coeffRef(i, 1)*=-1;
        area_matrix.coeffRef(i, i)=std::abs(area_matrix.row(i).sum());
    }

}

void HGF::volumes_and_areas(const Eigen::MatrixXd &targetU , Eigen::VectorXd &volumes, Eigen::MatrixXd& area_matrix){
        using namespace Eigen;
    
    volumes=VectorXd::Zero(num_bubbles);
    area_matrix=MatrixXd::Zero(num_bubbles,num_bubbles);

    int nt=(int)mesh().nt();
    for (size_t t = 0; t < nt; t++)
    {

        auto f=mesh().get_triangle(t);
        Eigen::Vector3d p0 = targetU.row(f.v[0]);
        Eigen::Vector3d p1 = targetU.row(f.v[1]);;
        Eigen::Vector3d p2 = targetU.row(f.v[2]);

        Eigen::Vector3d n = (p1-p0).cross(p2-p0);
        
        int label0=mesh().get_triangle_label(t)[0];
        int label1=mesh().get_triangle_label(t)[1];
        
        double local_area=n.norm();

        if(label0<=AIR and label1>AIR){
            area_matrix(label1-region_offset,label1-region_offset)+=local_area;
        }else if(label1<=AIR and label0>AIR){
            area_matrix(label0-region_offset,label0-region_offset)+=local_area;
        }else if(label0>AIR and label1>AIR){
            area_matrix(label0-region_offset,label1-region_offset)+=local_area;
            area_matrix(label1-region_offset,label0-region_offset)+=local_area;
        }

        // total volume via divergence theorem: ∫ 1
        double local_volume=n.dot(p0);
        if(label0>AIR){
            volumes[label0-region_offset]+=local_volume;
        }
        if(label1>AIR){
            volumes[label1-region_offset]-=local_volume;
            
        }
    }

    volumes/=6.0;
    area_matrix/=2.0;

    for(int i=0;i<num_bubbles;++i){
        area_matrix.row(i)*=-1;
        area_matrix(i,i)=std::abs(area_matrix.row(i).sum());
    }

}

void HGF::volumes_and_areas( Eigen::VectorXd &volumes, Eigen::MatrixXd& area_matrix){

    using namespace Eigen;
    
    volumes=VectorXd::Zero(num_bubbles);
    area_matrix=MatrixXd::Zero(num_bubbles,num_bubbles);

    int nt=(int)mesh().nt();
    for (size_t t = 0; t < nt; t++)
    {

        auto f=mesh().get_triangle(t);

        LosTopos::Vec3d p0 = m_st->pm_positions[f.v[0]];
        LosTopos::Vec3d  p1 = m_st->pm_positions[f.v[1]];
        LosTopos::Vec3d  p2 = m_st->pm_positions[f.v[2]];
        
        LosTopos::Vec3d n=LosTopos::cross(p1-p0, p2-p0);

        int label0=mesh().get_triangle_label(t)[0];
        int label1=mesh().get_triangle_label(t)[1];

        double local_area=LosTopos::mag(n);

        if(label0<=AIR and label1>AIR){
            area_matrix(label1-region_offset,label1-region_offset)+=local_area;
        }else if(label1<=AIR and label0>AIR){
            area_matrix(label0-region_offset,label0-region_offset)+=local_area;
        }else if(label0>AIR and label1>AIR){
            area_matrix(label0-region_offset,label1-region_offset)+=local_area;
            area_matrix(label1-region_offset,label0-region_offset)+=local_area;
        }
        
        // total volume via divergence theorem: ∫ 1
        double local_volume=LosTopos::dot(n,p0);

        if(label0>AIR){
            volumes[label0-region_offset]+=local_volume;
        }
        if(label1>AIR){
            volumes[label1-region_offset]-=local_volume;
            
        }
        
    }

    volumes/=6.0;
    
    area_matrix/=2.0;

    for(int i=0;i<num_bubbles;++i){
        area_matrix.row(i)*=-1;
        area_matrix(i,i)=std::abs(area_matrix.row(i).sum());
    }

}

void HGF::test_update_mesh_via_LosTopos(){
    
    int nv=mesh().nv();
    for (size_t i = 0; i < nv; i++){
        m_st->pm_velocities[i] = vc(vel(i));
    }
    
    m_st->topology_changes();
    m_st->improve_mesh();
    
    //Necessary for volume correction.
    easy_orientation();
    m_st->defrag_mesh_from_scratch(m_constrained_vertices);
}

