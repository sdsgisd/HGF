//
//  Scenes.cpp
//
//
//  Created by Fang Da on 15/1/26.
//
//  Editted by Sadashige Ishida 2017.

#include <map>
#include <iostream>
#include <cassert>
#include <cmath>
#include <random>

#include <igl/read_triangle_mesh.h>
#include <igl/writeOBJ.h>

#include "eigenheaders.h"
#include "Scenes.h"
#include "Options.h"
#include "MeshIO.h"
#include "Sim.h"
#include "HGF.h"

//static double frame_center_x;
double Scenes::frame_center_x;
double Scenes::d_of_rings;
double Scenes::R_of_rings;
double Scenes::frame_out;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Scene-specific initialization
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace
{
    class OrderComp
    {
    public:
        bool operator () (const std::pair<int, double> & o1, const std::pair<int, double> & o2) const { return o1.second < o2.second; }
    };
    
    // subdivide every triangle in the mesh into four triangles, by introducing three new vertices at edge midpoints
    // after subdivision, vs will be expanded (retaining all original vertices), while fs and ls will be replaced (no original faces remain)
    // if r is positive, new vertices will be projected onto the sphere centered at center, with a radius of r (for ICO sphere generation)
    void subdivide(const Vec3d & center, double r, std::vector<Vec3d> & vs, std::vector<Vec3i> & fs, std::vector<Vec2i> & ls)
    {
        std::vector<Vec3i> new_fs;
        std::vector<Vec2i> new_ls;
        std::map<std::pair<int, int>, int> new_vs_map;
        
        // map edge-midpoint to new vertices
        for (size_t j = 0; j < fs.size(); j++)
        {
            Vec3i & f = fs[j];
            int v0 = f[0];
            int v1 = f[1];
            int v2 = f[2];
            
            std::pair<int, int> p01 = (v0 < v1 ? std::pair<int, int>(v0, v1) : std::pair<int, int>(v1, v0));
            std::pair<int, int> p12 = (v1 < v2 ? std::pair<int, int>(v1, v2) : std::pair<int, int>(v2, v1));
            std::pair<int, int> p20 = (v2 < v0 ? std::pair<int, int>(v2, v0) : std::pair<int, int>(v0, v2));
            
            new_vs_map[p01] = 0;
            new_vs_map[p12] = 0;
            new_vs_map[p20] = 0;
        }
        
        // create the new vertices
        for (std::map<std::pair<int, int>, int>::iterator j = new_vs_map.begin(); j != new_vs_map.end(); j++)
        {
            j->second = vs.size();
            if (r > 0)
                vs.push_back(((vs[j->first.first] + vs[j->first.second]) / 2 - center).normalized() * r + center);
            else
                vs.push_back((vs[j->first.first] + vs[j->first.second]) / 2);
        }
        
        // triangulate
        for (size_t j = 0; j < fs.size(); j++)
        {
            Vec3i & f = fs[j];
            int v0 = f[0];
            int v1 = f[1];
            int v2 = f[2];
            
            std::pair<int, int> p01 = (v0 < v1 ? std::pair<int, int>(v0, v1) : std::pair<int, int>(v1, v0));
            std::pair<int, int> p12 = (v1 < v2 ? std::pair<int, int>(v1, v2) : std::pair<int, int>(v2, v1));
            std::pair<int, int> p20 = (v2 < v0 ? std::pair<int, int>(v2, v0) : std::pair<int, int>(v0, v2));
            int nv01 = new_vs_map[p01];
            int nv12 = new_vs_map[p12];
            int nv20 = new_vs_map[p20];
            
            Vec2i & l = ls[j];
            new_fs.push_back(Vec3i(v0, nv01, nv20));   new_ls.push_back(l);
            new_fs.push_back(Vec3i(nv01, v1, nv12));   new_ls.push_back(l);
            new_fs.push_back(Vec3i(nv20, nv12, v2));   new_ls.push_back(l);
            new_fs.push_back(Vec3i(nv12, nv20, nv01)); new_ls.push_back(l);
        }
        
        fs = new_fs;
        ls = new_ls;
    }
    
    void createIcoSphere(const Vec3d & center, double r, int subdivision, std::vector<Vec3d> & vs, std::vector<Vec3i> & fs, std::vector<Vec2i> & ls, const Vec2i & label = Vec2i(1, 0))
    {
        // create the initial icosahedron
        double phi = (1.0 + sqrt(5.0)) / 2.0;
        double len = Vec3d(0, 1, phi).norm();
        
        vs.resize(12);
        vs[ 0] = center + r * Vec3d(0,  1,  phi) / len;
        vs[ 1] = center + r * Vec3d(0, -1,  phi) / len;
        vs[ 2] = center + r * Vec3d(0,  1, -phi) / len;
        vs[ 3] = center + r * Vec3d(0, -1, -phi) / len;
        vs[ 4] = center + r * Vec3d( 1,  phi, 0) / len;
        vs[ 5] = center + r * Vec3d(-1,  phi, 0) / len;
        vs[ 6] = center + r * Vec3d( 1, -phi, 0) / len;
        vs[ 7] = center + r * Vec3d(-1, -phi, 0) / len;
        vs[ 8] = center + r * Vec3d( phi, 0,  1) / len;
        vs[ 9] = center + r * Vec3d( phi, 0, -1) / len;
        vs[10] = center + r * Vec3d(-phi, 0,  1) / len;
        vs[11] = center + r * Vec3d(-phi, 0, -1) / len;
        
        fs.push_back(Vec3i( 0,  1,  8));    ls.push_back(label);
        fs.push_back(Vec3i( 1,  0, 10));    ls.push_back(label);
        fs.push_back(Vec3i( 2,  3, 11));    ls.push_back(label);
        fs.push_back(Vec3i( 3,  2,  9));    ls.push_back(label);
        
        fs.push_back(Vec3i( 4,  5,  0));    ls.push_back(label);
        fs.push_back(Vec3i( 5,  4,  2));    ls.push_back(label);
        fs.push_back(Vec3i( 6,  7,  3));    ls.push_back(label);
        fs.push_back(Vec3i( 7,  6,  1));    ls.push_back(label);
        
        fs.push_back(Vec3i( 8,  9,  4));    ls.push_back(label);
        fs.push_back(Vec3i( 9,  8,  6));    ls.push_back(label);
        fs.push_back(Vec3i(10, 11,  7));    ls.push_back(label);
        fs.push_back(Vec3i(11, 10,  5));    ls.push_back(label);
        
        fs.push_back(Vec3i( 0,  8,  4));    ls.push_back(label);
        fs.push_back(Vec3i( 1,  6,  8));    ls.push_back(label);
        fs.push_back(Vec3i( 0,  5, 10));    ls.push_back(label);
        fs.push_back(Vec3i( 1, 10,  7));    ls.push_back(label);
        
        fs.push_back(Vec3i(11,  3,  7));    ls.push_back(label);
        fs.push_back(Vec3i( 5,  2, 11));    ls.push_back(label);
        fs.push_back(Vec3i( 6,  3,  9));    ls.push_back(label);
        fs.push_back(Vec3i( 9,  2,  4));    ls.push_back(label);
        
        for (int i = 0; i < subdivision; i++)
            subdivide(center, r, vs, fs, ls);
        
    }
    
}

HGF * Scenes::sceneSphere(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("num_subdivision");
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    double r = 1.0;
    createIcoSphere(Vec3d(0, 0, 0), r, N, v, f, l);
    
    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d (v[i][0] * 0.4, v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
    
    return new HGF(vs, fs, ls, cv, cx);
}

HGF * Scenes::sceneTet(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("num_subdivision");
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    double r = 1.0;
    double len = sqrt(3.0);
    v.push_back(Vec3d(-1, -1, -1) * r / len);
    v.push_back(Vec3d(-1,  1,  1) * r / len);
    v.push_back(Vec3d( 1, -1,  1) * r / len);
    v.push_back(Vec3d( 1,  1, -1) * r / len);
    
    f.push_back(Vec3i(0, 2, 1));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 3, 2));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 1, 3));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(1, 2, 3));    l.push_back(Vec2i(1, 0));
    
    for (int i = 0; i < N; i++)
        subdivide(Vec3d(0, 0, 0), 0, v, f, l);
    
    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
    
    return new HGF(vs, fs, ls, cv, cx);
}

HGF * Scenes::sceneCube(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("num_subdivision");
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    double r = 1.0;//2.0;1.0;
    
    double len = sqrt(3.0);
    v.push_back(Vec3d(-1, -1, -1) * r / len);
    v.push_back(Vec3d( 1, -1, -1) * r / len);
    v.push_back(Vec3d( 1,  1, -1) * r / len);
    v.push_back(Vec3d(-1,  1, -1) * r / len);
    v.push_back(Vec3d(-1, -1,  1) * r / len);
    v.push_back(Vec3d( 1, -1,  1) * r / len);
    v.push_back(Vec3d( 1,  1,  1) * r / len);
    v.push_back(Vec3d(-1,  1,  1) * r / len);
    
    f.push_back(Vec3i(0, 3, 2));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(2, 1, 0));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 1, 5));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(5, 4, 0));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(1, 2, 6));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(6, 5, 1));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(2, 3, 7));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(7, 6, 2));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(3, 0, 4));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(4, 7, 3));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(4, 5, 6));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(6, 7, 4));    l.push_back(Vec2i(1, 0));
    
    for (int i = 0; i < N; i++)
        subdivide(Vec3d(0, 0, 0), 0, v, f, l);
    
    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
    
    return new HGF(vs, fs, ls, cv, cx);
}

HGF * Scenes::sceneSheet(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("num_subdivision");
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    // square sheet
    //    v.push_back(Vec3d(-1,  0, -1));
    //    v.push_back(Vec3d( 1,  0, -1));
    //    v.push_back(Vec3d(-1,  0,  1));
    //    v.push_back(Vec3d( 1,  0,  1));
    //
    //    f.push_back(Vec3i(0, 1, 3));    l.push_back(Vec2i(0, 1));
    //    f.push_back(Vec3i(3, 2, 0));    l.push_back(Vec2i(0, 1));
    
    // hexagonal sheet
    double h = std::sqrt(3.0) / 2;
    v.push_back(Vec3d( 0,    0,  0));
    v.push_back(Vec3d(-1,    0,  0));
    v.push_back(Vec3d(-0.5,  0, -h));
    v.push_back(Vec3d( 0.5,  0, -h));
    v.push_back(Vec3d( 1,    0,  0));
    v.push_back(Vec3d( 0.5,  0,  h));
    v.push_back(Vec3d(-0.5,  0,  h));
    
    f.push_back(Vec3i(0, 1, 2));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 2, 3));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 3, 4));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 4, 5));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 5, 6));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 6, 1));    l.push_back(Vec2i(1, 0));
    
    for (int i = 0; i < N; i++)
        subdivide(Vec3d(0, 0, 0), 0, v, f, l);
    
    double displacement = 0.1;
    for (size_t i = 0; i < v.size(); i++)
        if ((v[i] - Vec3d(0, 0, 0)).norm() < 1e-6)
            v[i].y() += displacement;
    
    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
    
    //return new HGF(vs, fs, ls, cv, cx);
    
    return new HGF(vs, fs, ls, cv, cx);
}

HGF * Scenes::sceneBarrel(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("mesh-size-n");
    int M = Options::intValue("mesh-size-m");
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    double r = 1;
    double L = 1;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
            v.push_back(Vec3d(r * cos(j * 2 * M_PI / N), ((double)i / (M - 1) - 0.5) * L, -r * sin(j * 2 * M_PI / N)));
    
    for (int i = 0; i < M - 1; i++)
        for (int j = 0; j < N; j++)
        {
            f.push_back(Vec3i(i * N + j, i * N + (j + 1) % N, (i + 1) * N + (j + 1) % N));  l.push_back(Vec2i(1, 0));
            f.push_back(Vec3i((i + 1) * N + (j + 1) % N, (i + 1) * N + j, i * N + j));      l.push_back(Vec2i(1, 0));
        }
    
    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
    
    return new HGF(vs, fs, ls, cv, cx);
}

HGF * Scenes::sceneDoubleBubble(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("mesh-size-n");
    int M = Options::intValue("mesh-size-m");
    
    double r = 1.0; // radius of the two bubbles
    
    int overlap_angle = N / 2;
    //int overlap_angle = N / 8;
    double d = 2 * r * cos(M_PI * overlap_angle / N);
    
    vs.push_back(LosTopos::Vec3d(-d / 2 - r, 0, 0));    // left cap
    for (int i = 1; i < N + 1 - overlap_angle; i++)
        for (int j = 0; j < M; j++)
            vs.push_back(LosTopos::Vec3d(-d / 2 - r * cos(M_PI * i / N), r * sin(M_PI * i / N) * cos(2 * M_PI * j / M), r * sin(M_PI * i / N) * sin(2 * M_PI * j / M)));
    for (int i = 1; i < N - overlap_angle; i++)
        for (int j = 0; j < M; j++)
            vs.push_back(LosTopos::Vec3d( d / 2 + r * cos(M_PI * i / N), r * sin(M_PI * i / N) * cos(2 * M_PI * j / M), r * sin(M_PI * i / N) * sin(2 * M_PI * j / M)));
    size_t rid = vs.size();
    vs.push_back(LosTopos::Vec3d(d / 2 + r, 0, 0));     // right cap
    size_t cid = vs.size();
    vs.push_back(LosTopos::Vec3d(0, 0, 0));             // center vertex for triangulating the planar wall separating the two bubbles
    
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(0, 1 + j, 1 + (j + 1) % M)), ls.push_back(LosTopos::Vec2i(0, 1));
    for (int i = 1; i < N - overlap_angle; i++)
        for (int j = 0; j < M; j++)
            fs.push_back(LosTopos::Vec3st(1 + (i - 1) * M + j, 1 + i * M + j, 1 + i * M + (j + 1) % M)),                 ls.push_back(LosTopos::Vec2i(0, 1)),
            fs.push_back(LosTopos::Vec3st(1 + i * M + (j + 1) % M, 1 + (i - 1) * M + (j + 1) % M, 1 + (i - 1) * M + j)), ls.push_back(LosTopos::Vec2i(0, 1));
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(rid, 1 + (N - overlap_angle) * M + j, 1 + (N - overlap_angle) * M + (j + 1) % M)), ls.push_back(LosTopos::Vec2i(2, 0));
    for (int i = 1; i < N - overlap_angle - 1; i++)
        for (int j = 0; j < M; j++)
            fs.push_back(LosTopos::Vec3st(1 + (N - overlap_angle) * M + (i - 1) * M + j, 1 + (N - overlap_angle) * M + i * M + j, 1 + (N - overlap_angle) * M + i * M + (j + 1) % M)),                 ls.push_back(LosTopos::Vec2i(2, 0)),
            fs.push_back(LosTopos::Vec3st(1 + (N - overlap_angle) * M + i * M + (j + 1) % M, 1 + (N - overlap_angle) * M + (i - 1) * M + (j + 1) % M, 1 + (N - overlap_angle) * M + (i - 1) * M + j)), ls.push_back(LosTopos::Vec2i(2, 0));
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(1 + (N - overlap_angle - 1) * 2 * M + j, 1 + (N - overlap_angle - 1) * M + j, 1 + (N - overlap_angle - 1) * M + (j + 1) % M)),               ls.push_back(LosTopos::Vec2i(2, 0)),
        fs.push_back(LosTopos::Vec3st(1 + (N - overlap_angle - 1) * M + (j + 1) % M, 1 + (N - overlap_angle - 1) * 2 * M + (j + 1) % M, 1 + (N - overlap_angle - 1) * 2 * M + j)), ls.push_back(LosTopos::Vec2i(2, 0));
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(cid, 1 + (N - overlap_angle - 1) * M + (j + 1) % M, 1 + (N - overlap_angle - 1) * M + j)), ls.push_back(LosTopos::Vec2i(2, 1));
    
    return new HGF(vs, fs, ls, cv, cx);
}

HGF * Scenes::sceneTwoBubbles(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("num_subdivision");
    
    std::vector<Vec3d> v_all;
    std::vector<Vec3i> f_all;
    std::vector<Vec2i> l_all;
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    double r = 1.0;
    
    v.clear();
    f.clear();
    l.clear();
    createIcoSphere(Vec3d(-1.02, 0, 0), r, N, v, f, l, Vec2i(1, 0));
    f_all.reserve(f_all.size() + f.size());
    for (size_t i = 0; i < f.size(); i++)
        f_all.push_back(Vec3i(f[i][0] + v_all.size(), f[i][1] + v_all.size(), f[i][2] + v_all.size()));
    v_all.insert(v_all.end(), v.begin(), v.end());
    l_all.insert(l_all.end(), l.begin(), l.end());
    
    v.clear();
    f.clear();
    l.clear();
    createIcoSphere(Vec3d(1.02, 0, 0), r, N, v, f, l, Vec2i(2, 0));
    f_all.reserve(f_all.size() + f.size());
    for (size_t i = 0; i < f.size(); i++)
        f_all.push_back(Vec3i(f[i][0] + v_all.size(), f[i][1] + v_all.size(), f[i][2] + v_all.size()));
    v_all.insert(v_all.end(), v.begin(), v.end());
    l_all.insert(l_all.end(), l.begin(), l.end());
    
    vs.resize(v_all.size());
    fs.resize(f_all.size());
    ls.resize(l_all.size());
    for (size_t i = 0; i < v_all.size(); i++)
        vs[i] = LosTopos::Vec3d (v_all[i][0] * 0.7, v_all[i][1], v_all[i][2]);
    for (size_t i = 0; i < f_all.size(); i++)
        fs[i] = LosTopos::Vec3st(f_all[i][0], f_all[i][1], f_all[i][2]);
    for (size_t i = 0; i < l_all.size(); i++)
        ls[i] = LosTopos::Vec2i (l_all[i][0], l_all[i][1]);
    
    return new HGF(vs, fs, ls, cv, cx);
}

HGF * Scenes::sceneTripleJunction(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("num_subdivision");
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    double len = 1.0;
    double r = 1.0;
    v.push_back(Vec3d( 0,         0,                0));
    v.push_back(Vec3d( 0,       len,                0));
    v.push_back(Vec3d( r * 2,         0,                0));
    v.push_back(Vec3d( r * 2,       len,                0));
    v.push_back(Vec3d(-r * 0.5,   0, -r * 0.866025404));
    v.push_back(Vec3d(-r * 0.5, len, -r * 0.866025404));
    v.push_back(Vec3d(-r * 0.5,   0,  r * 0.866025404));
    v.push_back(Vec3d(-r * 0.5, len,  r * 0.866025404));
    
    f.push_back(Vec3i(0, 1, 2));    l.push_back(Vec2i(0, 1));
    f.push_back(Vec3i(3, 2, 1));    l.push_back(Vec2i(0, 1));
    f.push_back(Vec3i(0, 1, 4));    l.push_back(Vec2i(1, 2));
    f.push_back(Vec3i(5, 4, 1));    l.push_back(Vec2i(1, 2));
    f.push_back(Vec3i(0, 1, 6));    l.push_back(Vec2i(2, 0));
    f.push_back(Vec3i(7, 6, 1));    l.push_back(Vec2i(2, 0));
    
    for (int i = 0; i < N; i++)
        subdivide(Vec3d(0, 0, 0), 0, v, f, l);
    
    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
    
    return new HGF(vs, fs, ls, cv, cx);
}

HGF * Scenes::sceneFoamInit(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("num_subdivision");   // subdiv levels
    int M = Options::intValue("mesh-size-m");   // number of bubbles
    
    std::vector<Vec3d> v_all;
    std::vector<Vec3i> f_all;
    std::vector<Vec2i> l_all;
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    std::vector<std::pair<Vec3d, double> > spheres;
    std::string sub_scene=Options::strValue("sub_scene");
    
    const bool two_side=Options::boolValue("two_side");

    double distance_between_chank=-1;
    if(two_side){distance_between_chank=10;};

    double xyz_max=3;//11;
    double x_max=xyz_max,y_max=xyz_max,z_max=xyz_max;//2.0;

    if(sub_scene=="triple"){
        spheres.resize(3);
        spheres[0].first = Vec3d(1, 0, 0);
        spheres[0].second = 0.5;
        spheres[1].first = Vec3d(-0.5, 0, -0.866025404);
        spheres[1].second = 0.5;
        spheres[2].first = Vec3d(-0.5, 0,  0.866025404);
        spheres[2].second = 0.5;
    }else  if(sub_scene=="double"){
        spheres.resize(2);
        spheres[0].first = Vec3d(-0.55, 0, 0);
        spheres[0].second = 0.5;
        spheres[1].first = Vec3d(0.55, 0,  0);
        spheres[1].second = 0.5;
    }
    else{
        
        double rmin = 0.75;//0.75; 0.4;0.2;//0.2;
        double rmax = 1.3 ;// 1.0;//0.5;
        
        srand(0);
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < 10000; j++)
            {
                Vec3d c = Vec3d(x_max * rand() / RAND_MAX , y_max * rand() / RAND_MAX , z_max * rand() / RAND_MAX );

                double r = rmin + (rmax - rmin) * rand() / RAND_MAX;
                bool collision = false;
                for (size_t k = 0; k < spheres.size(); k++)
                {
                    double d = (c - spheres[k].first).norm();
                    if (d < r + spheres[k].second + 0.001)
                    {
                        collision = true;
                        break;
                    }
                }
                
                if (!collision)
                {
                    spheres.push_back(std::pair<Vec3d, double>(c, r));
                    break;
                }
            }
        }
    }

    std::cout << M << " spheres requested; " << spheres.size() << " actually spheres generated." << std::endl;
    
    for (size_t i = 0; i < spheres.size(); i++)
    {
        Vec3d c = spheres[i].first;
        double r = spheres[i].second;
        
        v.clear();
        f.clear();
        l.clear();
        createIcoSphere(c, r, N, v, f, l, Vec2i(i + 1, 0));
        
        if(two_side){
            int sign=i%2==0?1:-1;
            for(int vi=0;vi<v.size();++vi){
                v[vi][0]+=sign*distance_between_chank;
            }
        }
        
        f_all.reserve(f_all.size() + f.size());
        for (size_t i = 0; i < f.size(); i++)
            f_all.push_back(Vec3i(f[i][0] + v_all.size(), f[i][1] + v_all.size(), f[i][2] + v_all.size()));
        v_all.insert(v_all.end(), v.begin(), v.end());
        l_all.insert(l_all.end(), l.begin( ), l.end());
    }
    
    vs.resize(v_all.size());
    fs.resize(f_all.size());
    ls.resize(l_all.size());
    for (size_t i = 0; i < v_all.size(); i++){
        
        vs[i] = LosTopos::Vec3d (v_all[i][0]-x_max/2, v_all[i][1]-y_max/2, v_all[i][2]-z_max/2);

    }
    for (size_t i = 0; i < f_all.size(); i++)
        fs[i] = LosTopos::Vec3st(f_all[i][0], f_all[i][1], f_all[i][2]);
    for (size_t i = 0; i < l_all.size(); i++)
        ls[i] = LosTopos::Vec2i (l_all[i][0], l_all[i][1]);
    
    return new HGF(vs, fs, ls, cv, cx);
}

HGF * Scenes::sceneQuadJunction(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("num_subdivision");
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    double r = 1.0;
    v.push_back(Vec3d(0, 0, 0));
    
    v.push_back(Vec3d(-1,  0, -0.70710678) * r);
    v.push_back(Vec3d( 1,  0, -0.70710678) * r);
    v.push_back(Vec3d( 0, -1,  0.70710678) * r);
    v.push_back(Vec3d( 0,  1,  0.70710678) * r);
    
    v.push_back(Vec3d( 0,  0, -1.41421356) * r);
    v.push_back(Vec3d(-1, -1,           0) * r);
    v.push_back(Vec3d(-1,  1,           0) * r);
    v.push_back(Vec3d( 1, -1,           0) * r);
    v.push_back(Vec3d( 1,  1,           0) * r);
    v.push_back(Vec3d( 0,  0,  1.41421356) * r);
    
    f.push_back(Vec3i(0, 1, 5));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(5, 2, 0));    l.push_back(Vec2i(1, 0));
    f.push_back(Vec3i(0, 1, 6));    l.push_back(Vec2i(0, 2));
    f.push_back(Vec3i(6, 3, 0));    l.push_back(Vec2i(0, 2));
    f.push_back(Vec3i(0, 1, 7));    l.push_back(Vec2i(2, 1));
    f.push_back(Vec3i(7, 4, 0));    l.push_back(Vec2i(2, 1));
    f.push_back(Vec3i(0, 2, 8));    l.push_back(Vec2i(3, 0));
    f.push_back(Vec3i(8, 3, 0));    l.push_back(Vec2i(3, 0));
    f.push_back(Vec3i(0, 2, 9));    l.push_back(Vec2i(1, 3));
    f.push_back(Vec3i(9, 4, 0));    l.push_back(Vec2i(1, 3));
    f.push_back(Vec3i(0, 3, 10));   l.push_back(Vec2i(3, 2));
    f.push_back(Vec3i(10, 4, 0));   l.push_back(Vec2i(3, 2));
    
    double disp = 0.7;
    v[0] += v[1] * disp;
    
    for (int i = 0; i < N; i++)
        subdivide(Vec3d(0, 0, 0), 0, v, f, l);
    
    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
    
    return new HGF(vs, fs, ls, cv, cx);
}

HGF * Scenes::sceneConstrainedSphere(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("num_subdivision");
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    double r = 1.0;
    createIcoSphere(Vec3d(0, 0, 0), r, N, v, f, l);
    
    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d (v[i][0] * 0.6, v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
    
    cv.push_back(96);
    cv.push_back(101);
    cv.push_back(154);
    cv.push_back(155);
    cv.push_back(161);
    cv.push_back(160);
    cv.push_back(635);
    cv.push_back(641);
    cv.push_back(565);
    cv.push_back(566);
    cv.push_back(580);
    cv.push_back(581);
    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vc(vs[cv[i]]));
    
    return new HGF(vs, fs, ls, cv, cx);
}

HGF * Scenes::sceneBubbleWand(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("mesh-size-n");
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    double r = 1.0;
    v.push_back(Vec3d(0, 0, 0));
    v.push_back(Vec3d(-1, 0, 0));
    for (int i = 0; i < N; i++)
        v.push_back(Vec3d(0, r * cos(i * 2 * M_PI / N), r * sin(i * 2 * M_PI / N)));
    
    for (int i = 0; i < N; i++)
    {
        f.push_back(Vec3i(0, i + 2, (i + 1) % N + 2));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(1, i + 2, (i + 1) % N + 2));    l.push_back(Vec2i(0, 1));
    }
    
    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
    
    cv.push_back(1);
    for (int i = 0; i < N; i++)
        cv.push_back(i + 2);
    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vc(vs[cv[i]]));
    
    int num_bubbles=1;//0;
    
    return new HGF(vs, fs, ls, cv, cx,num_bubbles);
}

HGF * Scenes::sceneTwoRingsPinching(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("mesh-size-n");
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    R_of_rings = 0.5;
    d_of_rings = 0.28;//0.22;
    
    v.push_back(Vec3d(-1 - d_of_rings, 0, 0));
    v.push_back(Vec3d( 1 + d_of_rings, 0, 0));
    for (int i = 0; i < N; i++)
        v.push_back(Vec3d(-d_of_rings, 0, 0) + Vec3d(0, cos(i * 2 * M_PI / N), sin(i * 2 * M_PI / N)) * R_of_rings);
    for (int i = 0; i < N; i++)
        v.push_back(Vec3d( 0, 0, 0) + Vec3d(0, cos(i * 2 * M_PI / N), sin(i * 2 * M_PI / N)) * R_of_rings);
    for (int i = 0; i < N; i++)
        v.push_back(Vec3d( d_of_rings, 0, 0) + Vec3d(0, cos(i * 2 * M_PI / N), sin(i * 2 * M_PI / N)) * R_of_rings);
    
    for (int i = 0; i < N; i++)
    {
        f.push_back(Vec3i(0, i + 2, (i + 1) % N + 2));                      l.push_back(Vec2i(0, 1));
        f.push_back(Vec3i(1, i + 2 + 2 * N, (i + 1) % N + 2 + 2 * N));      l.push_back(Vec2i(1, 0));
        
        f.push_back(Vec3i(i + 2, (i + 1) % N + 2, (i + 1) % N + 2 + N));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i((i + 1) % N + 2 + N, i + 2 + N, i + 2));          l.push_back(Vec2i(1, 0));
        
        f.push_back(Vec3i(i + 2 + N, (i + 1) % N + 2 + N, (i + 1) % N + 2 + 2 * N));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i((i + 1) % N + 2 + 2 * N, i + 2 + 2 * N, i + 2 + N));          l.push_back(Vec2i(1, 0));
    }
    
    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
    
    cv.push_back(0);
    cv.push_back(1);
    for (int i = 0; i < N; i++)
        cv.push_back(i + 2);
    for (int i = 0; i < N; i++)
        cv.push_back(i + 2 + 2 * N);
    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vc(vs[cv[i]]));
    
    //If cmc is true, the volume enclosed by the film and the imaginary regions is preserved.
    //Hence, films converge to constant mean curvature (CMC) surfaces.
    bool cmc=false;
    int num_bubbles=cmc?1:0;

    HGF * hgf=new HGF(vs, fs, ls, cv, cx,num_bubbles);
    Scenes::double_torus(sim, hgf);
    
    return hgf;
}

void Scenes::setPulling(HGF *hgf){
    
    int num_from_end=4;//2;
    //How far constrained vertices from the left/right ends of the given foam.
    
    auto&  vs = hgf->m_st->pm_positions;
    auto& fs = hgf->m_st->m_mesh.m_tris;
    auto&  ls = hgf->m_st->m_mesh.m_triangle_labels;
    
    auto &cv=hgf->m_constrained_vertices;
    auto&cx=hgf->m_constrained_positions;
    assert(cv.empty());
    
    //    LosTopos::Vec3d left_axis(-1, -1, 0);   // for foam10
    LosTopos::Vec3d left_axis(-1, 0, 0);    // for foam3
    std::vector<std::pair<int, double> > left_order;
    for (size_t i = 0; i < vs.size(); i++)
        left_order.push_back(std::pair<int, double>(i, -dot(vs[i], left_axis)));
    OrderComp comp;
    std::sort(left_order.begin(), left_order.end(), comp);
    
    //    LosTopos::Vec3d right_axis(1, 1, 0.5);  // for foam10
    LosTopos::Vec3d right_axis(1, 0, 0);    // for foam3
    std::vector<std::pair<int, double> > right_order;
    for (size_t i = 0; i < vs.size(); i++)
        right_order.push_back(std::pair<int, double>(i, -dot(vs[i], right_axis)));
    std::sort(right_order.begin(), right_order.end(), comp);

    size_t leftmost = left_order[0].first;
    std::set<size_t> floodfill;
    floodfill.insert(leftmost);
    std::set<size_t> floodfill_next = floodfill;

    for (int k = 0; k < num_from_end; k++)
    {
        floodfill = floodfill_next;
        for (std::set<size_t>::iterator ii = floodfill.begin(); ii != floodfill.end(); ii++)
        {
            size_t i = *ii;
            for (size_t j = 0; j < hgf->mesh().m_vertex_to_edge_map[i].size(); j++)
            {
                LosTopos::Vec2st e = hgf->mesh().m_edges[hgf->mesh().m_vertex_to_edge_map[i][j]];
                size_t vother = (e[0] == i ? e[1] : e[0]);
                floodfill_next.insert(vother);
            }
        }
    }
    
    std::vector<size_t> last_ring;
    for (std::set<size_t>::iterator ii = floodfill_next.begin(); ii != floodfill_next.end(); ii++)
        if (floodfill.find(*ii) == floodfill.end())
            last_ring.push_back(*ii);
    
    Vec3d lr_center(0, 0, 0);
    Vec3d lr_normal(0, 0, 0);
    Mat3d lr_normal_sum = Mat3d::Zero();
    double lr_radius = 0;
    for (size_t i = 0; i < last_ring.size(); i++)
        lr_center += hgf->pos(last_ring[i]);
    lr_center /= last_ring.size();
    for (size_t i = 0; i < last_ring.size(); i++)
    {
        for (size_t j = i + 1; j < last_ring.size(); j++)
        {
            Vec3d n = (hgf->pos(last_ring[i]) - lr_center).cross(hgf->pos(last_ring[j]) - lr_center);
            lr_normal_sum += n * n.transpose();
        }
        lr_radius += (hgf->pos(last_ring[i]) - lr_center).norm();
    }
    lr_radius /= last_ring.size();
    Eigen::SelfAdjointEigenSolver<Mat3d> eig(lr_normal_sum);
    lr_normal = eig.eigenvectors().col(2).normalized();
    
    for (size_t i = 0; i < last_ring.size(); i++)
    {
        Vec3d x = hgf->pos(last_ring[i]);
        x -= (x - lr_center).dot(lr_normal) * lr_normal;
        x = (x - lr_center).normalized() * lr_radius + lr_center;
        hgf->m_st->pm_positions[last_ring[i]] = hgf->m_st->pm_newpositions[last_ring[i]] = vc(x);
    }

    for (size_t i = 0; i < last_ring.size(); i++)
        cv.push_back(last_ring[i]);
    
    size_t rightmost = right_order[0].first;
    floodfill.clear();
    floodfill.insert(rightmost);
    floodfill_next = floodfill;
    
    for (int k = 0; k < num_from_end; k++)
    {
        floodfill = floodfill_next;
        for (std::set<size_t>::iterator ii = floodfill.begin(); ii != floodfill.end(); ii++)
        {
            size_t i = *ii;
            for (size_t j = 0; j < hgf->mesh().m_vertex_to_edge_map[i].size(); j++)
            {
                LosTopos::Vec2st e = hgf->mesh().m_edges[hgf->mesh().m_vertex_to_edge_map[i][j]];
                size_t vother = (e[0] == i ? e[1] : e[0]);
                floodfill_next.insert(vother);
            }
        }
    }
    
    last_ring.clear();
    for (std::set<size_t>::iterator ii = floodfill_next.begin(); ii != floodfill_next.end(); ii++)
        if (floodfill.find(*ii) == floodfill.end())
            last_ring.push_back(*ii);
    
    lr_center = Vec3d(0, 0, 0);
    lr_normal = Vec3d(0, 0, 0);
    lr_normal_sum = Mat3d::Zero();
    lr_radius = 0;
    for (size_t i = 0; i < last_ring.size(); i++)
        lr_center += hgf->pos(last_ring[i]);
    lr_center /= last_ring.size();
    for (size_t i = 0; i < last_ring.size(); i++)
    {
        for (size_t j = i + 1; j < last_ring.size(); j++)
        {
            Vec3d n = (hgf->pos(last_ring[i]) - lr_center).cross(hgf->pos(last_ring[j]) - lr_center);
            lr_normal_sum += n * n.transpose();
        }
        lr_radius += (hgf->pos(last_ring[i]) - lr_center).norm();
    }
    lr_radius /= last_ring.size();
    eig = Eigen::SelfAdjointEigenSolver<Mat3d>(lr_normal_sum);
    lr_normal = eig.eigenvectors().col(2).normalized();
    
    for (size_t i = 0; i < last_ring.size(); i++)
    {
        Vec3d x = hgf->pos(last_ring[i]);
        x -= (x - lr_center).dot(lr_normal) * lr_normal;
        x = (x - lr_center).normalized() * lr_radius + lr_center;
        hgf->m_st->pm_positions[last_ring[i]] = hgf->m_st->pm_newpositions[last_ring[i]] = vc(x);
    }
    
    for (size_t i = 0; i < last_ring.size(); i++)
        cv.push_back(last_ring[i]);

    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(hgf->pos(cv[i]));
    
    hgf->m_constrained_vertices = cv;
    hgf->m_constrained_positions = cx;
    
    for (size_t i = 0; i < cv.size(); i++)
        hgf->m_st->m_masses[cv[i]] = LosTopos::Vec3d(1, 1, 1) * std::numeric_limits<double>::infinity();
    
}

HGF * Scenes::scenePeanutBubble(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("mesh-size-n");
    int M = Options::intValue("mesh-size-m");
    
    double r = 1.0; // radius of the two bubbles
    int noverlap = N / 3;
    double d = 2 * r * cos(M_PI * noverlap / N);
    
    vs.push_back(LosTopos::Vec3d(-d / 2 - r, 0, 0));    // left cap
    for (int i = 1; i < N + 1 - noverlap; i++)
        for (int j = 0; j < M; j++)
            vs.push_back(LosTopos::Vec3d(-d / 2 - r * cos(M_PI * i / N), r * sin(M_PI * i / N) * cos(2 * M_PI * j / M), r * sin(M_PI * i / N) * sin(2 * M_PI * j / M)));
    double right_ratio = 0.999;
    for (int i = 1; i < N - noverlap; i++)
        for (int j = 0; j < M; j++)
            vs.push_back(LosTopos::Vec3d(d / 2, 0, 0) + LosTopos::Vec3d(cos(M_PI * i / N), sin(M_PI * i / N) * cos(2 * M_PI * j / M), sin(M_PI * i / N) * sin(2 * M_PI * j / M)) * r * (right_ratio + (1 - right_ratio) * i / (N - noverlap)));
    size_t rid = vs.size();
    vs.push_back(LosTopos::Vec3d(d / 2 + r * right_ratio, 0, 0));     // right cap
    size_t cid = vs.size();
    vs.push_back(LosTopos::Vec3d(0, 0, 0));             // center vertex for triangulating the planar wall separating the two bubbles
    
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(0, 1 + j, 1 + (j + 1) % M)), ls.push_back(LosTopos::Vec2i(0, 1));
    for (int i = 1; i < N - noverlap; i++)
        for (int j = 0; j < M; j++)
            fs.push_back(LosTopos::Vec3st(1 + (i - 1) * M + j, 1 + i * M + j, 1 + i * M + (j + 1) % M)),                 ls.push_back(LosTopos::Vec2i(0, 1)),
            fs.push_back(LosTopos::Vec3st(1 + i * M + (j + 1) % M, 1 + (i - 1) * M + (j + 1) % M, 1 + (i - 1) * M + j)), ls.push_back(LosTopos::Vec2i(0, 1));
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(rid, 1 + (N - noverlap) * M + j, 1 + (N - noverlap) * M + (j + 1) % M)), ls.push_back(LosTopos::Vec2i(1, 0));
    for (int i = 1; i < N - noverlap - 1; i++)
        for (int j = 0; j < M; j++)
            fs.push_back(LosTopos::Vec3st(1 + (N - noverlap) * M + (i - 1) * M + j, 1 + (N - noverlap) * M + i * M + j, 1 + (N - noverlap) * M + i * M + (j + 1) % M)),                 ls.push_back(LosTopos::Vec2i(1, 0)),
            fs.push_back(LosTopos::Vec3st(1 + (N - noverlap) * M + i * M + (j + 1) % M, 1 + (N - noverlap) * M + (i - 1) * M + (j + 1) % M, 1 + (N - noverlap) * M + (i - 1) * M + j)), ls.push_back(LosTopos::Vec2i(1, 0));
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(1 + (N - noverlap - 1) * 2 * M + j, 1 + (N - noverlap - 1) * M + j, 1 + (N - noverlap - 1) * M + (j + 1) % M)),               ls.push_back(LosTopos::Vec2i(1, 0)),
        fs.push_back(LosTopos::Vec3st(1 + (N - noverlap - 1) * M + (j + 1) % M, 1 + (N - noverlap - 1) * 2 * M + (j + 1) % M, 1 + (N - noverlap - 1) * 2 * M + j)), ls.push_back(LosTopos::Vec2i(1, 0));
    
    for (size_t i = 0; i < M; i++)
        cv.push_back(1 + (N - noverlap - 1) * M + i);
    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vc(vs[cv[i]]));
    
    return new HGF(vs, fs, ls, cv, cx);
}

HGF * Scenes::sceneStraw(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("mesh-size-n");
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    double r = 1.0;
    v.push_back(Vec3d(0, 0, 0));
    v.push_back(Vec3d(-1, 0, 0));
    for (int i = 0; i < N; i++)
        v.push_back(Vec3d(0, r * cos(i * 2 * M_PI / N), r * sin(i * 2 * M_PI / N)));
    
    for (int i = 0; i < N; i++)
    {
        f.push_back(Vec3i(0, i + 2, (i + 1) % N + 2));    l.push_back(Vec2i(1, 0));
        f.push_back(Vec3i(1, i + 2, (i + 1) % N + 2));    l.push_back(Vec2i(0, 1));
    }
    
    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
    
    cv.push_back(1);
    for (int i = 0; i < N; i++)
        cv.push_back(i + 2);
    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vc(vs[cv[i]]));
    
    int num_bubbles=1;
    return new HGF(vs, fs, ls, cv, cx,num_bubbles);
}

HGF * Scenes::sceneCarousel(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("mesh-size-n");
    int M = Options::intValue("mesh-size-m");
    
    // double bubble
    double r = 4.0; // radius of the two bubbles
    int noverlap = N / 3;
    double d = 2 * r * cos(M_PI * noverlap / N);
    
    vs.push_back(LosTopos::Vec3d(-d / 2 - r, 0, 0));    // left cap
    for (int i = 1; i < N + 1 - noverlap; i++)
        for (int j = 0; j < M; j++)
            vs.push_back(LosTopos::Vec3d(-d / 2 - r * cos(M_PI * i / N), r * sin(M_PI * i / N) * cos(2 * M_PI * j / M), r * sin(M_PI * i / N) * sin(2 * M_PI * j / M)));
    for (int i = 1; i < N - noverlap; i++)
        for (int j = 0; j < M; j++)
            vs.push_back(LosTopos::Vec3d( d / 2 + r * cos(M_PI * i / N), r * sin(M_PI * i / N) * cos(2 * M_PI * j / M), r * sin(M_PI * i / N) * sin(2 * M_PI * j / M)));
    size_t rid = vs.size();
    vs.push_back(LosTopos::Vec3d(d / 2 + r, 0, 0));     // right cap
    size_t cid = vs.size();
    vs.push_back(LosTopos::Vec3d(0, 0, 0));             // center vertex for triangulating the planar wall separating the two bubbles
    
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(0, 1 + j, 1 + (j + 1) % M)), ls.push_back(LosTopos::Vec2i(0, 1));
    for (int i = 1; i < N - noverlap; i++)
        for (int j = 0; j < M; j++)
            fs.push_back(LosTopos::Vec3st(1 + (i - 1) * M + j, 1 + i * M + j, 1 + i * M + (j + 1) % M)),                 ls.push_back(LosTopos::Vec2i(0, 1)),
            fs.push_back(LosTopos::Vec3st(1 + i * M + (j + 1) % M, 1 + (i - 1) * M + (j + 1) % M, 1 + (i - 1) * M + j)), ls.push_back(LosTopos::Vec2i(0, 1));
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(rid, 1 + (N - noverlap) * M + j, 1 + (N - noverlap) * M + (j + 1) % M)), ls.push_back(LosTopos::Vec2i(2, 0));
    for (int i = 1; i < N - noverlap - 1; i++)
        for (int j = 0; j < M; j++)
            fs.push_back(LosTopos::Vec3st(1 + (N - noverlap) * M + (i - 1) * M + j, 1 + (N - noverlap) * M + i * M + j, 1 + (N - noverlap) * M + i * M + (j + 1) % M)),                 ls.push_back(LosTopos::Vec2i(2, 0)),
            fs.push_back(LosTopos::Vec3st(1 + (N - noverlap) * M + i * M + (j + 1) % M, 1 + (N - noverlap) * M + (i - 1) * M + (j + 1) % M, 1 + (N - noverlap) * M + (i - 1) * M + j)), ls.push_back(LosTopos::Vec2i(2, 0));
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(1 + (N - noverlap - 1) * 2 * M + j, 1 + (N - noverlap - 1) * M + j, 1 + (N - noverlap - 1) * M + (j + 1) % M)),               ls.push_back(LosTopos::Vec2i(2, 0)),
        fs.push_back(LosTopos::Vec3st(1 + (N - noverlap - 1) * M + (j + 1) % M, 1 + (N - noverlap - 1) * 2 * M + (j + 1) % M, 1 + (N - noverlap - 1) * 2 * M + j)), ls.push_back(LosTopos::Vec2i(2, 0));
    for (int j = 0; j < M; j++)
        fs.push_back(LosTopos::Vec3st(cid, 1 + (N - noverlap - 1) * M + (j + 1) % M, 1 + (N - noverlap - 1) * M + j)), ls.push_back(LosTopos::Vec2i(2, 1));

    for (size_t i = 0; i < M; i++)
        cv.push_back(1 + (N / 6) * M + i);
    for (size_t i = 0; i < M; i++)
        cv.push_back(1 + ((N - noverlap) + N / 6) * M + i);
    
    // straw
    int nstart_straw = vs.size();
    N /= 2;
    
    Mat3d R;
    R << 0, -1, 0, 1, 0, 0, 0, 0, 1;
    Vec3d t(0, -6, 0);
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    r = 1.0;
    v.push_back(Vec3d(0, 0, 0));
    v.push_back(Vec3d(-1, 0, 0));
    for (int i = 0; i < N; i++)
        v.push_back(Vec3d(0, r * cos(i * 2 * M_PI / N), r * sin(i * 2 * M_PI / N)));
    
    for (size_t i = 0; i < v.size(); i++)
        v[i] = R * v[i] + t;
    
    for (int i = 0; i < N; i++)
    {
        f.push_back(Vec3i(0, i + 2, (i + 1) % N + 2) + Vec3i(nstart_straw, nstart_straw, nstart_straw));    l.push_back(Vec2i(3, 0));
        f.push_back(Vec3i(1, i + 2, (i + 1) % N + 2) + Vec3i(nstart_straw, nstart_straw, nstart_straw));    l.push_back(Vec2i(0, 3));
    }
    
    for (size_t i = 0; i < v.size(); i++)
        vs.push_back(LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]));
    for (size_t i = 0; i < f.size(); i++)
        fs.push_back(LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]));
    for (size_t i = 0; i < l.size(); i++)
        ls.push_back(LosTopos::Vec2i (l[i][0], l[i][1]));
    
    cv.push_back(nstart_straw + 1);
    for (int i = 0; i < N; i++)
        cv.push_back(nstart_straw + i + 2);

    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vc(vs[cv[i]]));
    
    return new HGF(vs, fs, ls, cv, cx);
}

HGF * Scenes::sceneOctahedron(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    int N = Options::intValue("num_subdivision");
    
    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    Vec3d perturbation = 0.1 * Vec3d(2.0 * rand() / RAND_MAX - 1.0, 2.0 * rand() / RAND_MAX - 1.0, 2.0 * rand() / RAND_MAX - 1.0);
    
    double r = 1.0;
    v.push_back(Vec3d( 0,  0,  0) * r + perturbation);
    
    v.push_back(Vec3d( 1,  0,  0) * r);
    v.push_back(Vec3d(-1,  0,  0) * r);
    v.push_back(Vec3d( 0,  1,  0) * r);
    v.push_back(Vec3d( 0, -1,  0) * r);
    v.push_back(Vec3d( 0,  0,  1) * r);
    v.push_back(Vec3d( 0,  0, -1) * r);
    
    double a = 0.3;
    v.push_back(Vec3d( a,  a,  a) * r);
    v.push_back(Vec3d(-a, -a,  a) * r);
    v.push_back(Vec3d(-a,  a, -a) * r);
    v.push_back(Vec3d( a, -a, -a) * r);
    
    double rp = r * 0.5;
    v.push_back(Vec3d( 1,  1,  1) * rp);
    v.push_back(Vec3d( 1,  1, -1) * rp);
    v.push_back(Vec3d( 1, -1,  1) * rp);
    v.push_back(Vec3d( 1, -1, -1) * rp);
    v.push_back(Vec3d(-1,  1,  1) * rp);
    v.push_back(Vec3d(-1,  1, -1) * rp);
    v.push_back(Vec3d(-1, -1,  1) * rp);
    v.push_back(Vec3d(-1, -1, -1) * rp);
    
    f.push_back(Vec3i(0,  7,  8));    l.push_back(Vec2i(1, 2));
    f.push_back(Vec3i(0,  7,  9));    l.push_back(Vec2i(3, 1));
    f.push_back(Vec3i(0,  7, 10));    l.push_back(Vec2i(2, 3));
    f.push_back(Vec3i(0,  8,  9));    l.push_back(Vec2i(1, 4));
    f.push_back(Vec3i(0,  8, 10));    l.push_back(Vec2i(4, 2));
    f.push_back(Vec3i(0,  9, 10));    l.push_back(Vec2i(3, 4));
    
    f.push_back(Vec3i(5,  7,  8));    l.push_back(Vec2i(2, 1));
    f.push_back(Vec3i(3,  7,  9));    l.push_back(Vec2i(1, 3));
    f.push_back(Vec3i(1,  7, 10));    l.push_back(Vec2i(3, 2));
    f.push_back(Vec3i(2,  8,  9));    l.push_back(Vec2i(4, 1));
    f.push_back(Vec3i(4,  8, 10));    l.push_back(Vec2i(2, 4));
    f.push_back(Vec3i(6,  9, 10));    l.push_back(Vec2i(4, 3));
    
    f.push_back(Vec3i(1,  3,  7));    l.push_back(Vec2i(3, 5));
    f.push_back(Vec3i(2,  4,  8));    l.push_back(Vec2i(4, 6));
    f.push_back(Vec3i(2,  3,  9));    l.push_back(Vec2i(1, 7));
    f.push_back(Vec3i(1,  4, 10));    l.push_back(Vec2i(2, 8));
    
    f.push_back(Vec3i(1,  5,  7));    l.push_back(Vec2i(5, 2));
    f.push_back(Vec3i(3,  5,  7));    l.push_back(Vec2i(1, 5));
    f.push_back(Vec3i(2,  5,  8));    l.push_back(Vec2i(6, 1));
    f.push_back(Vec3i(4,  5,  8));    l.push_back(Vec2i(2, 6));
    f.push_back(Vec3i(2,  6,  9));    l.push_back(Vec2i(7, 4));
    f.push_back(Vec3i(3,  6,  9));    l.push_back(Vec2i(3, 7));
    f.push_back(Vec3i(1,  6, 10));    l.push_back(Vec2i(8, 3));
    f.push_back(Vec3i(4,  6, 10));    l.push_back(Vec2i(4, 8));
    
    f.push_back(Vec3i(1,  3, 11));    l.push_back(Vec2i(5, 0));
    f.push_back(Vec3i(3,  5, 11));    l.push_back(Vec2i(5, 0));
    f.push_back(Vec3i(5,  1, 11));    l.push_back(Vec2i(5, 0));
    f.push_back(Vec3i(1,  3, 12));    l.push_back(Vec2i(0, 3));
    f.push_back(Vec3i(3,  6, 12));    l.push_back(Vec2i(0, 3));
    f.push_back(Vec3i(6,  1, 12));    l.push_back(Vec2i(0, 3));
    f.push_back(Vec3i(1,  4, 13));    l.push_back(Vec2i(0, 2));
    f.push_back(Vec3i(4,  5, 13));    l.push_back(Vec2i(0, 2));
    f.push_back(Vec3i(5,  1, 13));    l.push_back(Vec2i(0, 2));
    f.push_back(Vec3i(1,  4, 14));    l.push_back(Vec2i(8, 0));
    f.push_back(Vec3i(4,  6, 14));    l.push_back(Vec2i(8, 0));
    f.push_back(Vec3i(6,  1, 14));    l.push_back(Vec2i(8, 0));
    f.push_back(Vec3i(2,  3, 15));    l.push_back(Vec2i(0, 1));
    f.push_back(Vec3i(3,  5, 15));    l.push_back(Vec2i(0, 1));
    f.push_back(Vec3i(5,  2, 15));    l.push_back(Vec2i(0, 1));
    f.push_back(Vec3i(2,  3, 16));    l.push_back(Vec2i(7, 0));
    f.push_back(Vec3i(3,  6, 16));    l.push_back(Vec2i(7, 0));
    f.push_back(Vec3i(6,  2, 16));    l.push_back(Vec2i(7, 0));
    f.push_back(Vec3i(2,  4, 17));    l.push_back(Vec2i(6, 0));
    f.push_back(Vec3i(4,  5, 17));    l.push_back(Vec2i(6, 0));
    f.push_back(Vec3i(5,  2, 17));    l.push_back(Vec2i(6, 0));
    f.push_back(Vec3i(2,  4, 18));    l.push_back(Vec2i(0, 4));
    f.push_back(Vec3i(4,  6, 18));    l.push_back(Vec2i(0, 4));
    f.push_back(Vec3i(6,  2, 18));    l.push_back(Vec2i(0, 4));
    
    for (int i = 0; i < N; i++)
        subdivide(Vec3d(0, 0, 0), 0, v, f, l);
    
    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++)
        vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
    for (size_t i = 0; i < f.size(); i++)
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);
    
    for (int i =  1; i <  7; i++) cv.push_back(i);
    for (int i = 11; i < 19; i++) cv.push_back(i);
    
    for (size_t i = 0; i < cv.size(); i++)
        cx.push_back(vc(vs[cv[i]]));
    
    int num_bubbles=0;
    return new HGF(vs, fs, ls, cv, cx,num_bubbles);
}

HGF * Scenes::sceneInputData(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx,const std::string inputdata_dir)
{
    
    if(not (sim->m_scene=="inputrec")){//equally, sim->m_scene=="inputmesh".
        int num_bubbles=-1;

        int N= Options::intValue("num_subdivision");

        using namespace Eigen;
        MatrixXd V;
        MatrixXi F;
        
        std::string meshfile;
        std::string labelfile;
        std::string constrained_file;
        
        std::string sub_scene=Options::strValue("sub_scene");
        
        meshfile=sub_scene+".obj";
        labelfile=sub_scene+"_flabel.txt";
        constrained_file=sub_scene+"_constvertex.txt";

        std::ifstream mesh_ifs(inputdata_dir+meshfile);
        if(!mesh_ifs){
            std::cout<<"No mesh file specified in \"sub_scene\" option."<<std::endl;
            exit(-1);
        }
        igl::read_triangle_mesh(inputdata_dir+meshfile, V, F);

        std::vector<Vec3d> v;
        std::vector<Vec3i> f;
        std::vector<Vec2i> l;
        
        int nv=V.rows();
        for(int vi=0;vi<nv;++vi){
            v.push_back(V.row(vi));
            
        }
        int nt=F.rows();

        for(int fi=0;fi<nt;++fi){
            f.push_back(F.row(fi));
        }
        
        std::ifstream label_ifs(inputdata_dir+labelfile);
        if(label_ifs){

            std::string str;
            getline(label_ifs,str);
            
            getline(label_ifs,str);
            int ntl=std::stoi(str);
            assert(nt=ntl);

            Vec2i flabel;
            int fcount=0;
            while(getline(label_ifs,str)){
                std::string region_str;
                std::istringstream stream(str);
                
                //seperate line by ' '.
                int region_num=0;
                while(getline(stream,region_str,' ')){
                    
                    int region=std::stoi(region_str);
                    //cout<<region<<",";
                    // FLabel(fcount,region_num)=region-1;
                    flabel[region_num]=region;
                    ++region_num;
                }
                
                l.push_back(flabel);
                ++fcount;
            }
        }else{
            std::cout<<"No label file."<<std::endl;
            
            for(int fi=0;fi<nt;++fi){
                l.push_back(Vec2i(1,0));
            }
        }
        
        std::ifstream constrained_ifs(inputdata_dir+constrained_file);
        
        if(constrained_ifs){
            
            if(N!=0){
                std::cout<<"Currently num_subdivision must be 0 with constrained vertices."<<std::endl;
                exit(-1);
            }

            std::string str;
            getline(constrained_ifs,str);
            int nt_cv=std::stoi(str);

            std::string constrained_vertex_str;
            
            while(getline(constrained_ifs,constrained_vertex_str)){

                int constrained_vertex=std::stoi(constrained_vertex_str);
                
                cv.push_back(constrained_vertex);
                cx.push_back(v[constrained_vertex]);

            }
        }else{
            std::cout<<"No constvertex file."<<std::endl;
        }

        for (int i = 0; i < N; i++)
            subdivide(Vec3d(0, 0, 0), 0, v, f, l);
        
        vs.resize(v.size());
        fs.resize(f.size());
        ls.resize(l.size());
        for (size_t i = 0; i < v.size(); i++)
            vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
        for (size_t i = 0; i < f.size(); i++)
            fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
        for (size_t i = 0; i < l.size(); i++)
            ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);

        return new HGF(vs,fs,ls,cv,cx,num_bubbles );

    }else{

        HGF * temp_hgf = new HGF(std::vector<LosTopos::Vec3d>(), std::vector<LosTopos::Vec3st>(), std::vector<LosTopos::Vec2i>());
        
        //        sim->m_load_directory = Options::strValue("load-dir");
        //        assert(sim->m_load_directory != "");

        std::string sub_scene=Options::strValue("sub_scene");
        MeshIO::load(*temp_hgf, inputdata_dir + sub_scene+".rec");

        vs = temp_hgf->m_st->pm_positions;
        fs = temp_hgf->m_st->m_mesh.m_tris;
        ls = temp_hgf->m_st->m_mesh.m_triangle_labels;

        return new HGF(vs, fs, ls, cv, cx);
        
    }

}

HGF * Scenes::sceneBrakke(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx,const std::string inputdata_dir)
{

    assert ( Options::intValue("num_subdivision")==0);
    HGF *hgf= sceneInputData(sim, vs, fs, ls, cv, cx,inputdata_dir);
    
    const bool transition_and_rotation=1;
    
    if(!transition_and_rotation){
        return hgf;
    }else{
        
        int nv=hgf->mesh().nv();
        int nt=hgf->mesh().nt();
        
        using namespace Eigen;
        MatrixXd V(nv,3);
        MatrixXi F(nt,3);
        
        for (int vi=0;vi<nv;++vi){
            V.row(vi)=vc(vs[vi]);
        }
        
        for (int ti=0;ti<nt;++ti){
            F(ti,0)=fs[ti][0];
            F(ti,1)=fs[ti][1];
            F(ti,2)=fs[ti][2];
            
        }

        using namespace Eigen;
        
        Translation<double, 3> translation = Translation<double, 3>(0.0, -1.0, -0.5);
        
        // Rotation by quarternion.
        Eigen::Quaterniond rotate(Eigen::AngleAxisd(0.5*M_PI, Eigen::Vector3d::UnitZ()));
        
        Eigen::Quaterniond rotate2(Eigen::AngleAxisd(0., Vector3d(1,1,1)));

        DiagonalMatrix<double, 3> scaling = Scaling(1.0, 1., 1.0);

        Eigen::Affine3d affine_matrix;

        affine_matrix =   scaling * rotate2*rotate*translation;
        
        MatrixXd newV(V.rows(),4);
        newV<<V,VectorXd::Constant(V.rows(), 1);

        auto &mat=affine_matrix.matrix();
        auto &tempMatrix=newV*mat.transpose();

        newV=tempMatrix;
        newV.conservativeResize(newV.rows(), 3);
        V=newV;

        vs.resize(V.rows());
        fs.resize(F.rows());
        ls.resize(F.rows());
        for (size_t i = 0; i < vs.size(); i++)
            vs[i]=vc(V.row(i));
        for (size_t i = 0; i < fs.size(); i++)
            fs[i]= LosTopos::Vec3st(F(i,0), F(i,1), F(i,2));
        for (size_t i = 0; i < ls.size(); i++)
            ls[i] = LosTopos::Vec2i (1,0);

        return new HGF(vs,fs,ls,cv,cx );
        
    }
}

HGF * Scenes::sceneCubicFrame(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx){

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;

    frame_center_x=0.0;
    frame_out=2.0;
    double frame_in=1.;
    double imag_vert_dist=5.0;
    
    //v0~7. vertices of the cubic frame.
    v.emplace_back(-frame_out,-frame_out,-frame_out);
    v.emplace_back(frame_out,-frame_out,-frame_out);
    v.emplace_back(frame_out,frame_out,-frame_out);
    v.emplace_back(-frame_out,frame_out,-frame_out);
    v.emplace_back(-frame_out,-frame_out,frame_out);
    v.emplace_back(frame_out,-frame_out,frame_out);
    v.emplace_back(frame_out,frame_out,frame_out);
    v.emplace_back(-frame_out,frame_out,frame_out);
    
    //v8~15. vertices of the film inside.
    v.emplace_back(-frame_in,-frame_in,-frame_in);
    v.emplace_back(frame_in,-frame_in,-frame_in);
    v.emplace_back(frame_in,frame_in,-frame_in);
    v.emplace_back(-frame_in,frame_in,-frame_in);
    v.emplace_back(-frame_in,-frame_in,frame_in);
    v.emplace_back(frame_in,-frame_in,frame_in);
    v.emplace_back(frame_in,frame_in,frame_in);
    v.emplace_back(-frame_in,frame_in,frame_in);
    
    //v16~v21. imaginary vertices for LosTopos.
    v.emplace_back(0,-imag_vert_dist,0);
    v.emplace_back(imag_vert_dist,0,0);
    v.emplace_back(0,imag_vert_dist,0);
    v.emplace_back(-imag_vert_dist,0,0);
    v.emplace_back(0,0,-imag_vert_dist);
    v.emplace_back(0,0,imag_vert_dist);

    //Inner cubic film
    f.push_back(Vec3i(8, 9, 12));    l.push_back(Vec2i(7, 1));
    f.push_back(Vec3i(12, 9, 13));    l.push_back(Vec2i(7, 1));
    f.push_back(Vec3i(9, 10, 14));    l.push_back(Vec2i(7, 2));
    f.push_back(Vec3i(9, 14, 13));    l.push_back(Vec2i(7, 2));
    f.push_back(Vec3i(10, 11, 15));    l.push_back(Vec2i(7, 3));
    f.push_back(Vec3i(10, 15, 14));    l.push_back(Vec2i(7, 3));
    f.push_back(Vec3i(11, 8, 12));    l.push_back(Vec2i(7, 4));
    f.push_back(Vec3i(11, 12, 15));    l.push_back(Vec2i(7, 4));
    f.push_back(Vec3i(9, 8, 11));    l.push_back(Vec2i(7, 5));
    f.push_back(Vec3i(9, 11, 10));    l.push_back(Vec2i(7, 5));
    f.push_back(Vec3i(14, 15, 13));    l.push_back(Vec2i(7, 6));
    f.push_back(Vec3i(13, 15, 12));    l.push_back(Vec2i(7, 6));
    
    //triangles connecting vertices of inside films and vertices of cubic frames.
    //Lower  part
    f.emplace_back(0,1,8);l.emplace_back(5,1);
    f.emplace_back(1,9,8);l.emplace_back(5,1);
    
    f.emplace_back(1,2,9);l.emplace_back(5,2);
    f.emplace_back(2,10,9);l.emplace_back(5,2);
    
    f.emplace_back(2,3,10);l.emplace_back(5,3);
    f.emplace_back(3,11,10);l.emplace_back(5,3);
    
    f.emplace_back(3,0,11);l.emplace_back(5,4);
    f.emplace_back(0,8,11);l.emplace_back(5,4);
    
    //Upper part
    f.emplace_back(5,4,12);l.emplace_back(6,1);
    f.emplace_back(13,5,12);l.emplace_back(6,1);
    
    f.emplace_back(6,5,13);l.emplace_back(6,2);
    f.emplace_back(14,6,13);l.emplace_back(6,2);

    f.emplace_back(7,6,14);l.emplace_back(6,3);
    f.emplace_back(15,7,14);l.emplace_back(6,3);
    
    f.emplace_back(4,7,15);l.emplace_back(6,4);
    f.emplace_back(12,4,15);l.emplace_back(6,4);
    
    //Side part
    f.emplace_back(1,5,13);l.emplace_back(2,1);
    f.emplace_back(13,9,1);l.emplace_back(2,1);
    
    f.emplace_back(2,6,14);l.emplace_back(3,2);
    f.emplace_back(14,10,2);l.emplace_back(3,2);
    
    f.emplace_back(3,7,15);l.emplace_back(4,3);
    f.emplace_back(15,11,3);l.emplace_back(4,3);
    
    f.emplace_back(4,0,12);l.emplace_back(4,1);
    f.emplace_back(8,12,0);l.emplace_back(4,1);

    //Triangles connecting real vertices and imaginary vertices.
    f.emplace_back(16,0,1);l.emplace_back(1,0);
    f.emplace_back(16,1,5);l.emplace_back(1,0);
    f.emplace_back(16,5,4);l.emplace_back(1,0);
    f.emplace_back(16,4,0);l.emplace_back(1,0);
    
    f.emplace_back(17,1,2);l.emplace_back(2,0);
    f.emplace_back(17,2,6);l.emplace_back(2,0);
    f.emplace_back(17,6,5);l.emplace_back(2,0);
    f.emplace_back(17,5,1);l.emplace_back(2,0);
    
    f.emplace_back(18,2,3);l.emplace_back(3,0);
    f.emplace_back(18,3,7);l.emplace_back(3,0);
    f.emplace_back(18,7,6);l.emplace_back(3,0);
    f.emplace_back(18,6,2);l.emplace_back(3,0);
    
    f.emplace_back(19,3,0);l.emplace_back(4,0);
    f.emplace_back(19,0,4);l.emplace_back(4,0);
    f.emplace_back(19,4,7);l.emplace_back(4,0);
    f.emplace_back(19,7,3);l.emplace_back(4,0);
    //
    f.emplace_back(20,1,0);l.emplace_back(5,0);
    f.emplace_back(20,0,3);l.emplace_back(5,0);
    f.emplace_back(20,3,2);l.emplace_back(5,0);
    f.emplace_back(20,2,1);l.emplace_back(5,0);
    //
    f.emplace_back(21,4,5);l.emplace_back(6,0);
    f.emplace_back(21,5,6);l.emplace_back(6,0);
    f.emplace_back(21,6,7);l.emplace_back(6,0);
    f.emplace_back(21,7,4);l.emplace_back(6,0);

    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++){
        vs[i] = LosTopos::Vec3d (v[i][0], v[i][1], v[i][2]);
    }
    for (size_t i = 0; i < f.size(); i++){
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);
    }
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);

    int nv=v.size();
    for(int vi=0;vi<nv;++vi){
        double x=vs[vi][0];
        double y=vs[vi][1];
        double z=vs[vi][2];
        
        const double eps=1e-6;
        if(std::abs(x)+eps>frame_out or std::abs(y)+eps>frame_out or std::abs(z)+eps>frame_out){
            cv.push_back(vi);
            
        }

    }

    for (size_t i = 0; i < cv.size(); i++)
    {cx.push_back(vc(vs[cv[i]]));}

    int num_bubbles=1;
    
    HGF * hgf=new HGF(vs, fs, ls, cv, cx,num_bubbles);
    Scenes::wire_frame(sim, hgf);
    
    return hgf;

}

HGF * Scenes::sceneLattice(Sim *sim, std::vector<LosTopos::Vec3d> &vs, std::vector<LosTopos::Vec3st> &fs, std::vector<LosTopos::Vec2i> &ls, std::vector<size_t> &cv, std::vector<Vec3d> &cx){
    
    int N = Options::intValue("num_subdivision");
    int M = Options::intValue("mesh-size-m");

    std::vector<Vec3d> v;
    std::vector<Vec3i> f;
    std::vector<Vec2i> l;
    
    double dx=0.1;
    double r = 1.0;
    double len = sqrt(3.0);
    dx=2*r/len;
    
    int num_x=5,num_y=5,num_z=5;

    int num_per_axis=M;
    
    num_x=num_y=num_z=num_per_axis;

    int array_x=num_x+1;
    int array_y=num_y+1;
    int array_z=num_z+1;

    for(int z=0;z<array_z;++z){
        for(int y=0;y<array_y;++y){
            for(int x=0;x<array_x;++x){
                v.emplace_back(x *dx,y*dx,z*dx);
            }
        }
    }
    
    auto convertXYZ2int=[=](int x,int y,int z)->int{
        return z* array_x*array_y+y*array_x+x;
    };
    
    auto convertXYZ2intForLabel=[=](int x,int y,int z)->int{
        return z* (array_x-1)*(array_y-1)+y*(array_x-1)+x;
    };
    
    //faces and labels pararell to z-axis.
    for(int z=0;z<array_z;++z){
        for(int y=0;y<array_y-1;++y){
            for(int x=0;x<array_x-1;++x){
                
                int v0=convertXYZ2int(x,y,z);
                int v1=convertXYZ2int(x+1,y,z);
                int v2=convertXYZ2int(x+1,y+1,z);
                int v3=convertXYZ2int(x,y+1,z);

                if(z==array_z-1){
                    f.emplace_back(v0,v1,v2);
                    f.emplace_back(v0,v2,v3);

                }else{
                    f.emplace_back(v1,v0,v2);
                    f.emplace_back(v2,v0,v3);
                }
                
                int region0,region1;
                
                if(z==0){
                    region0=1+convertXYZ2intForLabel(x,y,z);//0 is reserved for AIR, so +1.
                    region1=0;
                }else if(z==array_z-1){
                    
                    region0=1+convertXYZ2intForLabel(x,y,z-1);
                    region1=0;
                }else{
                    region0=1+convertXYZ2intForLabel(x,y,z);
                    region1=1+convertXYZ2intForLabel(x,y,z-1);
                    
                }
                
                l.emplace_back(region0,region1);
                l.emplace_back(region0,region1);
                
            }
        }
    }
    
    //faces and labels pararell to y-axis.
    for(int z=0;z<array_z-1;++z){
        for(int y=0;y<array_y;++y){
            for(int x=0;x<array_x-1;++x){
                
                int v0=convertXYZ2int(x,y,z);
                int v1=convertXYZ2int(x+1,y,z);
                int v2=convertXYZ2int(x+1,y,z+1);
                int v3=convertXYZ2int(x,y,z+1);
                
                if(y==array_y-1){
                    f.emplace_back(v1,v0,v2);
                    f.emplace_back(v2,v0,v3);
                    
                }else{

                    f.emplace_back(v0,v1,v2);
                    f.emplace_back(v0,v2,v3);
                }
                int region0,region1;
                if(y==0){
                    region0=1+convertXYZ2intForLabel(x,y,z);
                    region1=0;
                }else if(y==array_y-1){
                    region0=1+convertXYZ2intForLabel(x,y-1,z);
                    region1=0;
                    
                }else{
                    region0=1+convertXYZ2intForLabel(x,y,z);
                    region1=1+convertXYZ2intForLabel(x,y-1,z);
                    
                }
                
                l.emplace_back(region0,region1);
                l.emplace_back(region0,region1);
                
            }
        }
    }
    
    //faces and labels pararell to x-axis.
    for(int z=0;z<array_z-1;++z){
        for(int y=0;y<array_y-1;++y){
            for(int x=0;x<array_x;++x){
                
                int v0=convertXYZ2int(x,y,z);
                int v1=convertXYZ2int(x,y+1,z);
                int v2=convertXYZ2int(x,y+1,z+1);
                int v3=convertXYZ2int(x,y,z+1);
                
                if(x==array_x-1){                    f.emplace_back(v0,v1,v2);
                    f.emplace_back(v0,v2,v3);
                    
                }else{
                    f.emplace_back(v1,v0,v2);
                    f.emplace_back(v2,v0,v3);

                }
                int region0,region1;
                if(x==0){
                    region0=1+convertXYZ2intForLabel(x,y,z);
                    region1=0;
                }else if(x==array_x-1){
                    region0=1+convertXYZ2intForLabel(x-1,y,z);
                    region1=0;
                    
                }else{
                    region0=1+convertXYZ2intForLabel(x,y,z);
                    region1=1+convertXYZ2intForLabel(x-1,y,z);
                    
                }
                
                l.emplace_back(region0,region1);
                l.emplace_back(region0,region1);
                
            }
        }
    }
    
    for (int i = 0; i < N; i++)
        subdivide(Vec3d(0, 0, 0), 0, v, f, l);

    int nv=v.size();
    int nt=f.size();
    Eigen::MatrixXd V(nv,3);
    Eigen::MatrixXi F(nt,3);
    for(int vi=0;vi<nv;++vi){
        V.row(vi)<<v[vi][0],v[vi][1],v[vi][2];
    }
    
    for(int ti=0;ti<nt;++ti){
        F.row(ti)<<f[ti][0],f[ti][1],f[ti][2];
    }
    
    const bool rotation_and_transition=true;
    
    if(rotation_and_transition){
        
        using namespace Eigen;
        Translation<double, 3> translation = Translation<double, 3>(1.0, 0.0, 0.0);
        
        Eigen::Quaterniond rotate(Eigen::AngleAxisd(1.0, Eigen::Vector3d::UnitZ()));
        Eigen::Quaterniond rotate2(Eigen::AngleAxisd(0.2, Vector3d(1,1,1)));

        DiagonalMatrix<double, 3> scaling = Scaling(1.0, 1., 1.0);

        Eigen::Affine3d affine_matrix;
        
        affine_matrix = translation * scaling * rotate2*rotate;

        MatrixXd newV(V.rows(),4);
        newV<<V,VectorXd::Constant(V.rows(), 1);

        auto &mat=affine_matrix.matrix();
        
        std::cout<<mat.rows()<<":"<<mat.cols()<<std::endl;

        auto &tempMatrix=newV*mat.transpose();

        newV=tempMatrix;
        newV.conservativeResize(newV.rows(), 3);
        V=newV;

    }

    vs.resize(v.size());
    fs.resize(f.size());
    ls.resize(l.size());
    for (size_t i = 0; i < v.size(); i++){
        vs[i]=vc(V.row(i));

    }
    for (size_t i = 0; i < f.size(); i++){
        fs[i] = LosTopos::Vec3st(f[i][0], f[i][1], f[i][2]);

    }
    for (size_t i = 0; i < l.size(); i++)
        ls[i] = LosTopos::Vec2i (l[i][0], l[i][1]);

    int num_bubbles=(array_z-1)*(array_y-1)*(array_x-1);
    return new HGF(vs, fs, ls, cv, cx);

}

HGF * Scenes::sceneCubeOverFilm(Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx)
{
    
    std::vector<Vec3d> v_cubes;
    std::vector<Vec3i> f_cubes;
    std::vector<Vec2i> l_cubes;

    std::vector<Eigen::Vector3d> centers;
    std::vector<double> radi;

    centers.emplace_back(0,0,0);
    radi.push_back(1.);
    
    centers.emplace_back(0,0.5,2);
    radi.push_back(0.8);
    
    centers.emplace_back(1.5,1.5,1.5);
    radi.push_back(0.5);
    
    centers.emplace_back(-2,-2,1.4);
    radi.push_back(0.7);
    
    int num_bubbles=centers.size();
    
    int subdiv=Options::intValue("num_subdivision");

    for(int bi=0;bi<num_bubbles;++bi){
        std::vector<Vec3d> v_single_cube;
        std::vector<Vec3i> f_single_cube;
        std::vector<Vec2i> l_single_cube;
        createIcoSphere(centers[bi], radi[bi], subdiv, v_single_cube, f_single_cube, l_single_cube, Vec2i(2+bi, 0));

        f_cubes.reserve(f_cubes.size() + f_cubes.size());
        for (size_t i = 0; i < f_single_cube.size(); i++)
            f_cubes.push_back(Vec3i(f_single_cube[i][0] + v_cubes.size(), f_single_cube[i][1] + v_cubes.size(), f_single_cube[i][2] + v_cubes.size()));
        v_cubes.insert(v_cubes.end(), v_single_cube.begin(), v_single_cube.end());
        l_cubes.insert(l_cubes.end(), l_single_cube.begin( ), l_single_cube.end());
        
    }

    std::vector<Vec3d> v_film;
    std::vector<Vec3i> f_film;
    std::vector<Vec2i> l_film;
    
    double out=4;//5;
    double height=-2.2;//-4.0;//-2.0;
    double imag_vertex_height=height-5.0;
    
    v_film.emplace_back(-out,-out,height);
    v_film.emplace_back(out,-out,height);
    v_film.emplace_back(out,out,height);
    v_film.emplace_back(-out,out,height);
    
    v_film.emplace_back(0,0,height);//center
    v_film.emplace_back(0,0,imag_vertex_height);//imag_vertex
    
    for(int i=0;i<4;++i){
        f_film.emplace_back(i,(i+1)%4,4);
        l_film.emplace_back(1,0);
        
        f_film.emplace_back(i,(i+1)%4,5);
        l_film.emplace_back(1,0);
    }

    vs.resize(v_cubes.size()+v_film.size());
    fs.resize(f_cubes.size()+f_film.size());
    ls.resize(l_cubes.size()+l_film.size());
    
    for (size_t i = 0; i < v_cubes.size(); i++){
        vs[i]=vc(v_cubes[i]);
    }
    for (size_t i = 0; i < f_cubes.size(); i++){
        fs[i]=LosTopos::Vec3st(f_cubes[i][0], f_cubes[i][1], f_cubes[i][2]);
    }
    
    for (size_t i = 0; i < l_cubes.size(); i++){
        ls[i]=LosTopos::Vec2i(l_cubes[i][0], l_cubes[i][1]);
    }
    
    int v_cube_offset=v_cubes.size();
    int f_cube_offset=f_cubes.size();
    int l_cube_offset=l_cubes.size();
    
    for (size_t i = 0; i < v_film.size(); i++){
        vs[i+v_cube_offset]=vc(v_film[i]);
    }
    for (size_t i = 0; i < f_film.size(); i++){
        fs[i+f_cube_offset]=LosTopos::Vec3st(f_film[i][0]+v_cube_offset, f_film[i][1]+v_cube_offset, f_film[i][2]+v_cube_offset);
    }
    
    for (size_t i = 0; i < l_film.size(); i++){
        ls[i+l_cube_offset]=LosTopos::Vec2i(l_film[i][0], l_film[i][1]);
    }

    int nv=vs.size();
    const double eps=1e-5;
    for(int i=0;i<nv;++i){
        if(vs[i][2]-eps<imag_vertex_height){
            cv.push_back(i );
        }
        else if(vs[i][2]-eps<height and (std::abs(vs[i][0])+eps>out or std::abs(vs[i][1])+eps>out)){
            
            cv.push_back(i );
        }
    }
    
    for (size_t i = 0; i < cv.size(); i++)
    {
        cx.push_back(vc(vs[cv[i]]));
    }

    HGF * hgf=new HGF(vs, fs, ls, cv, cx,num_bubbles);
    Scenes::square_frame(sim,hgf,out,height);
    return hgf;
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Scene-specific time stepping
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Scenes::stepSphere(double dt, Sim * sim, HGF * hgf)
{
    
}

void Scenes::stepTet(double dt, Sim * sim, HGF * hgf)
{
    
}

void Scenes::stepCube(double dt, Sim * sim, HGF * hgf)
{
    
}

void Scenes::stepSheet(double dt, Sim * sim, HGF * hgf)
{
    
}

void Scenes::stepBarrel(double dt, Sim * sim, HGF * hgf)
{
    
}

void Scenes::stepDoubleBubble(double dt, Sim * sim, HGF * hgf)
{
    
}

void Scenes::stepTwoBubbles(double dt, Sim * sim, HGF * hgf)
{
    
}

void Scenes::stepTripleJunction(double dt, Sim * sim, HGF * hgf)
{
    
}

void Scenes::stepFoamInit(double dt, Sim * sim, HGF * hgf)
{
    
    //Grow bubbles either via evolving the surface toward the normals
    //evolving by increasing enclosed volumes.
    const bool growBubblesViaNormal=true;
    const bool growBubblesViaVolume=!growBubblesViaNormal;
    
    if(growBubblesViaNormal){
        static const double alpha = 0.1;
        
        std::vector<Vec3d> n(hgf->mesh().nt(), Vec3d(0, 0, 0));
        for (size_t i = 0; i < hgf->mesh().nt(); i++)
        {
            LosTopos::Vec3st t = hgf->mesh().get_triangle(i);
            Vec3d x0 = hgf->pos(t[0]);
            Vec3d x1 = hgf->pos(t[1]);
            Vec3d x2 = hgf->pos(t[2]);
            LosTopos::Vec2i l = hgf->mesh().get_triangle_label(i);
            if (l[0] == 0)
                n[i] = -(x1 - x0).cross(x2 - x0).normalized();
            else if (l[1] == 0)
                n[i] = (x1 - x0).cross(x2 - x0).normalized();
        }
        
        std::vector<Vec3d> displacement(hgf->mesh().nv(), Vec3d(0, 0, 0));
        for (size_t i = 0; i < hgf->mesh().nv(); i++)
        {
            int counter = 0;
            for (size_t j = 0; j < hgf->mesh().m_vertex_to_triangle_map[i].size(); j++)
            {
                LosTopos::Vec2i l = hgf->mesh().get_triangle_label(hgf->mesh().m_vertex_to_triangle_map[i][j]);
                if (l[0] == 0 || l[1] == 0)
                {
                    displacement[i] += n[hgf->mesh().m_vertex_to_triangle_map[i][j]];
                    counter++;
                }
            }
            if (counter == 0)
                displacement[i] = Vec3d(0, 0, 0);
            else
                displacement[i] /= counter;
            
            hgf->surfTrack()->pm_newpositions[i] = hgf->surfTrack()->pm_positions[i] + dt * alpha * vc(displacement[i]);
        }
        
        double actual_dt;
        hgf->surfTrack()->integrate(dt, actual_dt);

        Eigen::MatrixXd area_matrix;
        hgf->volumes_and_areas( hgf->initialVolumes, area_matrix);

        hgf->m_st->topology_changes();
        hgf->m_st->improve_mesh();
        hgf->m_st->defrag_mesh_from_scratch(hgf->m_constrained_vertices);
    }else if(growBubblesViaVolume){
        
        int num_bubbles=hgf->num_bubbles;
        
        for(int bi=0;bi<num_bubbles;++bi){
            hgf->initialVolumes[bi]*=1.01;
        }
        
    }

}

void Scenes::burstBubbles(double dt, Sim * sim, HGF * hgf)
{
    
    auto burst_a_bubble=[hgf](){
        std::cout << "Bursting a random bubble now." << std::endl;
        
        std::set<int> burstable_regions_set;
        for (size_t i = 0; i < hgf->mesh().nt(); i++)
        {
            LosTopos::Vec2i l = hgf->mesh().get_triangle_label(i);
            if (l[0] == 0 || l[1] == 0)
                burstable_regions_set.insert(l[0] == 0 ? l[1] : l[0]);
        }
        
        std::vector<int> burstable_regions;
        burstable_regions.assign(burstable_regions_set.begin(), burstable_regions_set.end());
        std::cout << "Eligible regions: ";
        for (size_t i = 0; i < burstable_regions.size(); i++)
            std::cout << burstable_regions[i] << " ";
        std::cout << std::endl;

        int region_to_burst = burstable_regions[rand() % burstable_regions.size()];
        std::cout << "Region to be bursted: " << region_to_burst << std::endl;

        // first, find all the triple junctions that would become manifold curves because of the deletion of the bursted region. because
        //  the other two wings around those triple junctions which were originally unrelated (i.e. having different arbitrary constants)
        //  would come together and form one continuous potential field, they need to be brought to terms with each other.
        // specifically, bursting region A will remove faces with label (0, A), and for each region B that A is adjacent to, faces with labels
        //  (0, B) and (A, B) will be connected into one manifold patch.
        int A = region_to_burst;
        std::set<int> Bs;
        for (size_t i = 0; i < hgf->mesh().nt(); i++)
        {
            LosTopos::Vec2i l = hgf->mesh().get_triangle_label(i);
            if (l[0] == region_to_burst || l[1] == region_to_burst)
                Bs.insert(l[0] == region_to_burst ? l[1] : l[0]);
        }

        for (size_t i = 0; i < hgf->mesh().nt(); i++)
        {
            LosTopos::Vec2i l = hgf->mesh().get_triangle_label(i);
            if (l[0] == region_to_burst || l[1] == region_to_burst)
            {
                int lother = (l[0] == region_to_burst ? l[1] : l[0]);
                if (lother == 0)
                    hgf->surfTrack()->remove_triangle(i);
                else
                    hgf->mesh().set_triangle_label(i, (l[0] == region_to_burst ? LosTopos::Vec2i(0, l[1]) : LosTopos::Vec2i(l[0], 0)));
            }
        }
        
        hgf->surfTrack()->defrag_mesh_from_scratch(hgf->m_constrained_vertices);
        
        std::cout << "Bubble bursted." << std::endl;
        
    };
    
    if(Options::boolValue("auto-burst")){
        static double s_next_burst = Options::doubleValue("auto-burst-start");
        if (sim->m_time >= s_next_burst)
        {
            burst_a_bubble();
            s_next_burst += Options::doubleValue("auto-burst-interval");
            
        }
        
    }
    
    if(hgf->bursting){
        burst_a_bubble();
    }
    
}

void Scenes::stepQuadJunction(double dt, Sim * sim, HGF * hgf)
{
    
}

void Scenes::stepConstrainedSphere(double dt, Sim * sim, HGF * hgf)
{
    
}

void Scenes::stepBubbleWand(double dt, Sim * sim, HGF * hgf)
{
    
    for (size_t i = 0; i < hgf->m_constrained_vertices.size(); i++)
        hgf->m_constrained_positions[i] += Vec3d(-1, 0, 0) * dt;
    
    double r=1.0;
    double circle_area=r*r*M_PI;
    
    if(sim->m_time>0.5){
        hgf->initialVolumes[0]+=circle_area*dt;
    }
    return;

    static double wand_velocity=0;
    
    const double wand_acceleration=10;
    wand_velocity-=dt*wand_acceleration;
    
    const double move_velocity=-1;wand_velocity;-10;-1;
    
    //    for (size_t i = 0; i < hgf->m_constrained_vertices.size(); i++)
    //    {hgf->m_constrained_positions[i].y() *=0.99;
    //        hgf->m_constrained_positions[i].z() *=0.99;
    //    }

    for (size_t i = 0; i < hgf->m_constrained_vertices.size(); i++)
    {hgf->m_constrained_positions[i] += Vec3d(move_velocity, 0, 0) * dt;
        
    }
}

void Scenes::stepTwoRingsPinching(double dt, Sim * sim, HGF * hgf)
{
    double pulling_velocity=0.04;//0.02;0.04;
    double move_duration=5;//3.;//2;//5;10;5;
    if (sim->m_time < move_duration){
        for (size_t i = 0; i < hgf->m_constrained_vertices.size(); i++){
            hgf->m_constrained_positions[i] += Vec3d((hgf->m_constrained_positions[i].x() > frame_center_x ? 1 : -1), 0, 0) * pulling_velocity * dt;
        }
        
        auto &CV=hgf->Constrained_V;
        for(int cvi=0,n_cv=CV.rows();cvi<n_cv;++cvi){
            CV(cvi,0)+=CV(cvi,0)>frame_center_x?pulling_velocity * dt:-pulling_velocity * dt;
        }
        
        d_of_rings+=pulling_velocity * dt;
    }
    
    //Twist
    bool twist =false;
    if(twist){
        double max_x=d_of_rings;
        const double max_angle=0.001*2*M_PI;
        
        //auto &CV=hgf->Constrained_V;
        int n_cv=hgf->constrainedVertices().size();
        auto & cp=hgf->constrainedPositions();
        for(int i=0;i<n_cv;++i){
            double x=cp[i].x();
            double y=cp[i].y();
            double z=cp[i].z();
            double angle=max_angle*(x-frame_center_x)/max_x;
            
            cp[i]<<x,cos(angle)*y-sin(angle)*z,sin(angle)*y+cos(angle)*z;
            //double angle=max_angle*CV(i,0)/max_x;
            
            //double y=CV(i,1);
            //double z=CV(i,2);
            
            //CV(i,1)=cos(angle)*y-sin(angle)*z;
            //CV(i,2)=sin(angle)*y+cos(angle)*z;
        }
    }
}

void Scenes::pullBubbles(double dt, Sim * sim, HGF * hgf)
{
    Vec3d center0(0, 0, 0);
    int counter0 = 0;
    Vec3d center1(0, 0, 0);
    int counter1 = 0;
    
    for (size_t i = 0; i < hgf->m_constrained_positions.size(); i++)
    {
        if (hgf->m_constrained_positions[i].x() > 0)
        {
            center0 += hgf->m_constrained_positions[i];
            counter0++;
        } else
        {
            center1 += hgf->m_constrained_positions[i];
            counter1++;
        }
    }
    
    center0 /= counter0;
    center1 /= counter1;
    
    Vec3d dir = (center0 - center1).normalized();
    
    double pulling_velocity=0.1;//0.02;
    
    for (size_t i = 0; i < hgf->m_constrained_vertices.size(); i++)
        hgf->m_constrained_positions[i] += (hgf->m_constrained_positions[i].x() > 0 ? 1 : -1) * dir * pulling_velocity * dt;
}

void Scenes::stepPeanutBubble(double dt, Sim * sim, HGF * hgf)
{
    
}

void Scenes::stepStraw(double dt, Sim * sim, HGF * hgf)
{
    for (size_t i = 0; i < hgf->m_constrained_vertices.size(); i++)
    {
        size_t cv = hgf->m_constrained_vertices[i];
        
    }
}

void Scenes::stepCarousel(double dt, Sim * sim, HGF * hgf)
{

    for (size_t i = 0; i < hgf->m_constrained_vertices.size(); i++)
    {
        size_t cv = hgf->m_constrained_vertices[i];
        if (hgf->m_constrained_positions[i][1] < -4)
            hgf->m_constrained_positions[i][1] -= 0.2 * dt;
    }
    
}

void Scenes::stepOctahedron(double dt, Sim * sim, HGF * hgf)
{
    
}

void Scenes::stepCubicFrame(double dt, Sim * sim, HGF * hgf)
{

    //Pulling;
    bool pulling =0;
    if(pulling){
        double pulling_velocity=0.05;

        if (sim->m_time < 4){
            for (size_t i = 0; i < hgf->m_constrained_vertices.size(); i++){
                
                if(hgf->m_constrained_positions[i].x()>frame_center_x){
                    hgf->m_constrained_positions[i] +=Vec3d(1,0,0)*pulling_velocity*dt;

                }else if
                    (hgf->m_constrained_positions[i].x()<frame_center_x){
                        hgf->m_constrained_positions[i] -=Vec3d(1,0,0)*pulling_velocity*dt;
                    }
                
            }
            
            auto &CV=hgf->Constrained_V;
            for(int cvi=0,n_cv=CV.rows();cvi<n_cv;++cvi){
                CV(cvi,0)+=CV(cvi,0)>0?pulling_velocity * dt:-pulling_velocity * dt;
            }
        }
        
        frame_out+=pulling_velocity * dt;

    }

    //Parallel Transition
    bool transitioning=0;
    if(transitioning){
        
        double transition_velocity=.1;//0.1;0.05;
        Vec3d direction;
        Vec3d transition;
        if (sim->m_time < 100){
            direction<<1,0,1;
            transition=transition_velocity*direction;
        }else if (sim->m_time < 10){
            direction<<2,0,-2;
            transition=transition_velocity*direction;
        }else{
            transition<<0,0,0;
        }
        
        for (size_t i = 0; i < hgf->m_constrained_vertices.size(); i++){
            hgf->m_constrained_positions[i] += transition * dt;}
        
    }
    
    //Twist
    bool twist=true;
    if(twist and sim->m_time < 7.5){
        double max_x=frame_out;
        const double max_angle=0.0001*2*M_PI;
        
        //auto &CV=hgf->Constrained_V;
        int n_cv=hgf->constrainedVertices().size();
        auto & cp=hgf->constrainedPositions();
        for(int i=0;i<n_cv;++i){
            double x=cp[i].x();
            double y=cp[i].y();
            double z=cp[i].z();
            double angle=max_angle*(x-frame_center_x)/max_x;
            
            cp[i]<<x,cos(angle)*y-sin(angle)*z,sin(angle)*y+cos(angle)*z;
            
        }
        
        auto &CV=hgf->Constrained_V;
        int N_CV=CV.rows();
        for(int I=0;I<N_CV;++I){
            double angle=max_angle*(CV(I,0)-frame_center_x)/max_x;
            
            double Y=CV(I,1);
            double Z=CV(I,2);
            
            CV(I,1)=cos(angle)*Y-sin(angle)*Z;
            CV(I,2)=sin(angle)*Y+cos(angle)*Z;
        }
    }
    
}

void Scenes::stepCubeOverFilm(double dt, Sim * sim, HGF * hgf)
{

    static bool air_resistance=false;
    double reduceScale=.2;

    if(air_resistance){
        using namespace Eigen;
        
        int nt=hgf->mesh().nt();
        for(int ti=0;ti<nt;++ti){
            
            size_t label0=hgf->mesh().m_triangle_labels[ti][0];
            size_t label1=hgf->mesh().m_triangle_labels[ti][1];
            
            if((label0==1 and label1==2) or (label0==2 and label1==1))
            {
                air_resistance=false;
                return;
            }

        }

        int nv=hgf->mesh().nv();
        std::vector<LosTopos::Vec2i> Vlabels(nv);

        for(int vi=0;vi<nv;++vi){

            assert(hgf->mesh().m_vertex_to_triangle_map[vi].size()>0);
            size_t tri=hgf->mesh().m_vertex_to_triangle_map[vi][0];
            Vlabels[vi]=hgf->mesh().get_triangle_label(tri);
            
            size_t label0=Vlabels[vi][0];
            size_t label1=Vlabels[vi][1];
            
            if((label0==0 and label1==2) or (label0==2 and label1==0))
            {
                //Reduce velocity by air resistance.
                auto& velocity=hgf->vel(vi);
                Vector3d wholeVelocity(0,0,-1);

                Vector3d vertex_normal=vc(hgf->m_st->get_vertex_normal(vi));
                
                double spatial_innerProduct=vertex_normal.dot(wholeVelocity);

                if(spatial_innerProduct<=0){
                    continue;
                }
                
                double innerProduct=velocity.dot(wholeVelocity);
                
                velocity-=reduceScale*innerProduct* wholeVelocity;

            }
            
        }
        
    }

}

void Scenes::stepBrakke             (double dt, Sim * sim, HGF * hgf){
    std::cout<<"m_scene brakke is only for checking the geometric properties of Brakke's surface evolver."<<std::endl;
    exit(0);
}

void Scenes::moveLeftOrRight(double dt, Sim * sim, HGF * hgf){

    double move_velocity=0.1;
    
    if(hgf->move_left){
        for (size_t i = 0; i < hgf->m_constrained_vertices.size(); i++)
        {hgf->m_constrained_positions[i] -= Vec3d(move_velocity, 0, 0) *  dt;
            
        }
        
        hgf->Constrained_V.col(0)-=Eigen::VectorXd::Constant(hgf->Constrained_V.size(), dt*move_velocity);
        
        frame_center_x-=dt*move_velocity;
        
    }else if(hgf->move_right){
        for (size_t i = 0; i < hgf->m_constrained_vertices.size(); i++)
        {hgf->m_constrained_positions[i] += Vec3d(move_velocity, 0, 0) *  dt;

        }
        hgf->Constrained_V.col(0)+=Eigen::VectorXd::Constant(hgf->Constrained_V.size(), dt*move_velocity);
        
        frame_center_x+=dt*move_velocity;
        
    }
    
}

void Scenes::volume_change(double dt, Sim * sim, HGF * hgf){
    
    bool blow_and_absorb=false;
    if(blow_and_absorb and hgf->absorbing_bubble0){
        
        hgf->initialVolumes[0]-=0.02;
        hgf->initialVolumes[1]+=0.02;
        return;
    }

    if(hgf->blowing_bubble0){
        
        hgf->initialVolumes[0]+=0.01;
        
    }else if(hgf->absorbing_bubble0){
        hgf->initialVolumes[0]-=0.01;
    }
}

void Scenes::set_initial_bubble_velocities( Sim *sim, HGF *hgf){
    using namespace Eigen;
    
    hgf->velocity_per_bubble.clear();
    
    using namespace Eigen;
    
    std::random_device rd;
    
    std::mt19937 mt(rd());
    double rand_scale=0.00025;//0.0001;0.0005;//0.03;//0.02;
    std::uniform_real_distribution<double> rand_v(-1.0*rand_scale,1.0*rand_scale);
    double z_offset=0.0005;//0.001;0.05;
    
    for(int ri=0;ri<hgf->num_region;++ri){
        if(ri<=hgf->AIR){
            hgf->velocity_per_bubble.push_back(std::make_unique<Vector3d>(0,0,0));
        }else{
            
            double vx=rand_v(mt);
            double vy=rand_v(mt);
            double vz=rand_v(mt);
            vz+=z_offset;
            
            hgf->velocity_per_bubble.push_back(std::make_unique<Vector3d>(vx,vy,vz));
        }
    }
    
}

void Scenes::advect_vertices(double dt, Sim *sim, HGF *hgf){
    using namespace Eigen;
    
    int nv=hgf->mesh().nv();
    
    for(int i=0;i<nv;++i){
        double & x_position=hgf->surfTrack()->pm_positions[i][0];
        const double force_scale=0.0001;
        const double z_offset = 0.2;
        const double z_force_scale = 0.5;
        
        const double cylinder_end=0;
        if(std::abs(x_position)>cylinder_end){
            if(x_position > 0){
                x_position -= (x_position)*(x_position*force_scale);
            }else{
                x_position += (x_position)*(x_position*force_scale);
                //x_position *= 0.99;
                //hgf->vel(i).x() += 10*force_scale;
            }
        }
        
        double & z_position = hgf->surfTrack()->pm_positions[i][2];
        //hgf->vel(i).z() += std::min(0.04, force_scale + z_offset);
        z_position += std::min(0.05, z_force_scale/(10.0+std::abs(x_position)));
        
    }
    return;
    
}

void Scenes::add_velocities_to_bubbles(double dt, Sim *sim, HGF *hgf){
    using namespace Eigen;
    
    bool toCenter=0;
    if(toCenter){
        int nv=hgf->mesh().nv();
        for(int i=0;i<nv;++i){
            double x_position=hgf->surfTrack()->pm_positions[i][0];
            const double force_scale=0.01;
            
            const double cylinder_end=4;
            if(std::abs(x_position)>cylinder_end){
                hgf->vel(i).x()-=force_scale;
            }else{
                
            }

        }
        return;
    }
    
    //Run through triangles.
    //Give velocity to each vertex according to the labels of the triangles containing the vertex.
    
    int nt=hgf->mesh().nt();
    
    for(size_t ti=0;ti<nt;++ti){
        
        auto labels=hgf->mesh().get_triangle_label(ti);
        int label0=labels[0];
        int label1=labels[1];
        
        auto tri=hgf->mesh().get_triangle(ti);

        size_t v0 = tri[0];
        size_t v1 = tri[1];
        size_t v2 = tri[2];

        int AIR=hgf->AIR;
        
        Vector3d added_velocity;
        if(1){//Add the velocity of the smaller region (excluding AIR).
            if(label0>AIR and label1>AIR){
                continue;
            }
            int target_region=label0>AIR?label0:label1;
            added_velocity=*(hgf->velocity_per_bubble)[target_region];
            
        }
        
        else if(0){//Add the sum of the velocities.
            added_velocity=*(hgf->velocity_per_bubble)[label0]+*(hgf->velocity_per_bubble)[label1];
            
        }else if(0){//Add the mean of the velocities.
            added_velocity=*(hgf->velocity_per_bubble)[label0]+*(hgf->velocity_per_bubble)[label1];
            added_velocity*=0.5;
            
        }
        
        hgf->mesh().m_vertex_to_edge_map[v0].size();
        
        hgf->vel(v0)+=added_velocity/hgf->mesh().m_vertex_to_edge_map[v0].size();
        hgf->vel(v1)+=added_velocity/hgf->mesh().m_vertex_to_edge_map[v1].size();;
        hgf->vel(v2)+=added_velocity/hgf->mesh().m_vertex_to_edge_map[v2].size();;
    }
}

void Scenes::give_large_velocities_to_bubbles(double dt, Sim *sim, HGF *hgf){
    
    //Give large velocities to vertices according to the label of the triangles containing them.
    using namespace Eigen;
    std::vector<Vector3d> velocities(hgf->num_region);
    
    std::random_device rd;
    
    std::mt19937 mt(rd());
    double rand_scale=0.05;//0.03;//0.02;
    std::uniform_real_distribution<double> rand_v(-1.0*rand_scale,1.0*rand_scale);
    double z_offset=0.1;//0.05;
    
    for(int ri=0;ri<hgf->num_region;++ri){
        if(ri<=hgf->AIR){
            velocities[ri]<<0,0,0;
            
        }else{
            
            double vx=rand_v(mt);
            double vy=rand_v(mt);
            double vz=rand_v(mt);
            vz+=z_offset;
            velocities[ri]<<vx,vy,vz;
        }
    }

    //Run through triangles.
    //Give velocity to each vertex according to the labels of the triangles containing the vertex.
    int nt=hgf->mesh().nt();
    for(size_t ti=0;ti<nt;++ti){
        
        auto labels=hgf->mesh().get_triangle_label(ti);
        int label0=labels[0];
        int label1=labels[1];
        
        auto tri=hgf->mesh().get_triangle(ti);

        size_t v0 = tri[0];
        size_t v1 = tri[1];
        size_t v2 = tri[2];
        
        hgf->vel(v0)+=velocities[label0];
        hgf->vel(v0)+=velocities[label1];
        
        hgf->vel(v1)+=velocities[label0];
        hgf->vel(v1)+=velocities[label1];
        
        hgf->vel(v2)+=velocities[label0];
        hgf->vel(v2)+=velocities[label1];
        
    }

}

// Clamp the first argument to be greater than or equal to the second
// and less than or equal to the third.
double Clamp(double val, double min, double max) {
    if (val < min) {
        val = min;
    }
    if (val > max) {
        val = max;
    }
    return val;
}

// Return true if the first value is within epsilon of the second value.
bool NearByMargin(double actual, double expected) {
    double diff = actual - expected;
    if (diff < 0.0) {
        diff = -diff;
    }
    // 5 bits of error in mantissa (source of '32 *')
    return diff < 32 * std::numeric_limits<double>::epsilon();
}

void ToSphericalCoords(const Eigen::Vector3d& dir, double* phi, double* theta) {

    //assert(NearByMargin(dir.squaredNorm(), 1.0));//check if dir is unit vector.
    // Explicitly clamp the z coordinate so that numeric errors don't cause it
    // to fall just outside of acos' domain.
    *theta = acos(Clamp(dir.z(), -1.0, 1.0));
    // We don't need to divide dir.y() or dir.x() by sin(theta) since they are
    // both scaled by it and atan2 will handle it appropriately.
    *phi = atan2(dir.y(), dir.x());
}

template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl_cylinder_arranged(
                                      const double height,
                                      const int axis_devisions,
                                      const int height_devisions,
                                      Eigen::PlainObjectBase<DerivedV> & V,
                                      Eigen::PlainObjectBase<DerivedF> & F)
{
    const double radius=.02;
    //const double height=1.0;
    
    V.resize(axis_devisions*height_devisions+2,3);
    F.resize(2*(axis_devisions*(height_devisions-1))+2*axis_devisions,3);
    int f = 0;
    typedef typename DerivedV::Scalar Scalar;
    for(int th = 0;th<axis_devisions;th++)
    {
        Scalar x = radius*cos(2.*M_PI*Scalar(th)/Scalar(axis_devisions));
        Scalar y = radius*sin(2.*M_PI*Scalar(th)/Scalar(axis_devisions));
        for(int h = 0;h<height_devisions;h++)
        {
            Scalar z = height*Scalar(h)/Scalar(height_devisions-1);
            V(th+h*axis_devisions,0) = x;
            V(th+h*axis_devisions,1) = y;
            V(th+h*axis_devisions,2) = z;
            if(h > 0)
            {
                F(f,0) = ((th+0)%axis_devisions)+(h-1)*axis_devisions;
                F(f,1) = ((th+0)%axis_devisions)+(h+0)*axis_devisions;
                F(f,2) = ((th+1)%axis_devisions)+(h-1)*axis_devisions;
                f++;
                F(f,0) = ((th+1)%axis_devisions)+(h-1)*axis_devisions;
                F(f,1) = ((th+0)%axis_devisions)+(h+0)*axis_devisions;
                F(f,2) = ((th+1)%axis_devisions)+(h+0)*axis_devisions;
                f++;
            }

        }
    }
    
    const double eps=0.01;
    V.row(V.rows()-2)<<0,0,-eps;
    V.row(V.rows()-1)<<0,0,height+eps;
    for(int th = 0;th<axis_devisions;th++)
    {
        
        F(f,0) =V.rows()-2 ;
        F(f,1) =th;
        F(f,2) = (th+1)%axis_devisions;
        f++;
        
        F(f,0) =V.rows()-1 ;
        F(f,1) = (th+1)%axis_devisions+(height_devisions-1)*axis_devisions;;
        F(f,2) =th+(height_devisions-1)*axis_devisions;
        f++;
        
    }
    assert(f == F.rows());
    
}

void cylinder(const Eigen::Vector3d from, const Eigen::Vector3d till, Eigen::MatrixXd &V_cyl,Eigen::MatrixXi &F_cyl){
    using namespace Eigen;

    Vector3d directed_edge=till-from;
    
    const int axis_devision=10;
    const int height_devisions=5;
    const double height=directed_edge.norm();
    igl_cylinder_arranged(height,axis_devision, height_devisions,V_cyl,F_cyl);
    
    Vector3d dir=directed_edge.normalized();
    double theta;
    double phi;
    ToSphericalCoords(dir, &phi, &theta);

    MatrixXd newV(V_cyl.rows(),4);
    newV<<V_cyl,VectorXd::Constant(V_cyl.rows(), 1);

    Eigen::Quaterniond rotateY(Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitY()));//Rotation for theta
    
    Eigen::Quaterniond rotateZ(Eigen::AngleAxisd(phi, Eigen::Vector3d::UnitZ()));//Rotation for phi

    DiagonalMatrix<double, 3> scaling = Scaling(1.0, 1.0, 1.0);
    
    Translation<double, 3> translation = Translation<double, 3>(from);

    Eigen::Affine3d affine_matrix;
    
    affine_matrix = translation  * rotateZ*rotateY *scaling;
    
    auto mat=affine_matrix.matrix();

    auto tempMatrix=newV*mat.transpose();
    
    newV=tempMatrix;
    newV.conservativeResize(newV.rows(), 3);
    V_cyl=newV;
}

void Scenes::wire_frame(Sim *sim, HGF *hgf){
    using namespace Eigen;
    using namespace std;
    class edge{
    public: int v0,v1;
        edge(int v0,int v1):v0(v0),v1(v1){
            
        }
    };

    vector<unique_ptr<Vector3d>> vertices;
    vector<unique_ptr<edge>> edges;
    
    vertices.push_back(make_unique<Vector3d>(-frame_out,-frame_out,-frame_out));
    vertices.push_back(make_unique<Vector3d>(frame_out,-frame_out,-frame_out));
    vertices.push_back(make_unique<Vector3d>(frame_out,frame_out,-frame_out));
    vertices.push_back(make_unique<Vector3d>(-frame_out,frame_out,-frame_out));
    vertices.push_back(make_unique<Vector3d>(-frame_out,-frame_out,frame_out));
    vertices.push_back(make_unique<Vector3d>(frame_out,-frame_out,frame_out));
    vertices.push_back(make_unique<Vector3d>(frame_out,frame_out,frame_out));
    vertices.push_back(make_unique<Vector3d>(-frame_out,frame_out,frame_out));
    
    edges.push_back(make_unique<edge>(0,1));
    edges.push_back(make_unique<edge>(1,2));
    edges.push_back(make_unique<edge>(2,3));
    edges.push_back(make_unique<edge>(3,0));
    
    edges.push_back(make_unique<edge>(4,5));
    edges.push_back(make_unique<edge>(5,6));
    edges.push_back(make_unique<edge>(6,7));
    edges.push_back(make_unique<edge>(7,4));
    
    edges.push_back(make_unique<edge>(0,4));
    edges.push_back(make_unique<edge>(1,5));
    edges.push_back(make_unique<edge>(2,6));
    edges.push_back(make_unique<edge>(3,7));
    
    auto &V=hgf->Constrained_V;
    auto &F=hgf->Constrained_F;
    
    for(int ei=0;ei<edges.size();++ei){
        MatrixXd V_cyl;
        MatrixXi F_cyl;
        cylinder(*vertices[edges[ei]->v0],*vertices[edges[ei]->v1],V_cyl,F_cyl);
        
        if(ei==0){
            V=V_cyl;
            F=F_cyl;
        }else{
            // Concatenate (VA,FA) and (VB,FB) into (V,F)
            MatrixXd tempV;
            MatrixXi tempF;
            tempV.resize(V.rows()+V_cyl.rows(),V.cols());
            tempV<<V,V_cyl;
            tempF.resize(F.rows()+F_cyl.rows(),F.cols());
            tempF<<F,(F_cyl.array()+V.rows());
            //
            V=tempV;
            F=tempF;
        }
    }
    
}

void Scenes::square_frame(Sim *sim, HGF *hgf,double out,double height){
    using namespace Eigen;
    using namespace std;
    class edge{
    public: int v0,v1;
        edge(int v0,int v1):v0(v0),v1(v1){
            
        }
    };

    vector<unique_ptr<Vector3d>> vertices;
    vector<unique_ptr<edge>> edges;
    
    vertices.push_back(make_unique<Vector3d>(-out,-out,height));
    vertices.push_back(make_unique<Vector3d>(out,-out,height));
    vertices.push_back(make_unique<Vector3d>(out,out,height));
    vertices.push_back(make_unique<Vector3d>(-out,out,height));

    edges.push_back(make_unique<edge>(0,1));
    edges.push_back(make_unique<edge>(1,2));
    edges.push_back(make_unique<edge>(2,3));
    edges.push_back(make_unique<edge>(3,0));

    auto &V=hgf->Constrained_V;
    auto &F=hgf->Constrained_F;
    
    for(int ei=0;ei<edges.size();++ei){
        MatrixXd V_cyl;
        MatrixXi F_cyl;
        cylinder(*vertices[edges[ei]->v0],*vertices[edges[ei]->v1],V_cyl,F_cyl);
        
        if(ei==0){
            V=V_cyl;
            F=F_cyl;
        }else{
            // Concatenate (VA,FA) and (VB,FB) into (V,F)
            MatrixXd tempV;
            MatrixXi tempF;
            tempV.resize(V.rows()+V_cyl.rows(),V.cols());
            tempV<<V,V_cyl;
            tempF.resize(F.rows()+F_cyl.rows(),F.cols());
            tempF<<F,(F_cyl.array()+V.rows());
            //
            V=tempV;
            F=tempF;
        }
    }
    
}

// Globals.
static int p = 72; // Number of grid columns.
static int q = 16; // Number of grid rows

// Fuctions to map the grid vertex (u_i,v_j) to the mesh vertex (f(u_i,v_j), g(u_i,v_j), h(u_i,v_j)) on the torus.
float f(double R, double r,int i, int j)
{
    return ( ( R + r * cos( (-1 + 2*(float)j/q) * M_PI ) ) * cos( (-1 + 2*(float)i/p) * M_PI ) );
}

float g(double R, double r,int i, int j)
{
    return ( ( R + r * cos( (-1 + 2*(float)j/q) * M_PI ) ) * sin( (-1 + 2*(float)i/p) * M_PI ) );
}

float h(double R, double r,int i, int j)
{
    return ( r * sin( (-1 + 2*(float)j/q) * M_PI ) );
}

void torus(const double R, const double r, Eigen::MatrixXd &V,Eigen::MatrixXi &F){
    // Routine to fill the vertex array with co-ordinates of the mapped sample points.
    int i, j, k;
    
    k = 0;

    V.resize((q+1)*(p+1),3);
    //F.resize(2*(axis_devisions*(height_devisions-1))+2*axis_devisions,3);
    for (j = 0; j <= q; j++){
        for (i = 0; i <= p; i++)
        {
            V.row(k)<<h(R,r,i,j),f(R,r,i,j),g(R,r,i,j);
            ++k;
        }
    }
    
    F.resize(2*p*q,3);
    // Make the approximating triangular mesh.
    for(j = 0; j < q; j++)
    {
        int index=2*p*j;

        int v0=-1,v1=-1,v2=-1;
        
        auto push=[&](int v){
            v2=v1;
            v1=v0;
            v0=v;
        };
        for(i = 0; i <= p; i++)
        {
            int v=(j+1)*(p+1) + i;

            push(v);
            if(i>0){
                F.row(index)<<v0,v1,v2;
                ++index;
                
            }

            v=j*(p+1) + i;

            push(v);
            if(i>0){
                F.row(index)<<v1,v0,v2;
                ++index;
                
            }

        }
    }
    //
    
}

void Scenes::double_torus(Sim *sim, HGF *hgf){
    using namespace Eigen;
    auto &V=hgf->Constrained_V;
    auto &F=hgf->Constrained_F;
    
    MatrixXd V_torus0,V_torus1;
    MatrixXi F_torus0,F_torus1;
    double R=0.5;
    
    double r=0.01;
    torus(R, r, V_torus0, F_torus0);
    torus(R, r, V_torus1, F_torus1);
    
    double d=d_of_rings;
    int num_v=V_torus1.rows();
    MatrixXd transition=MatrixXd::Zero(num_v, 3);
    transition<<VectorXd::Constant(num_v, d) ,MatrixXd::Zero(num_v,2);
    
    V_torus0-=transition;
    V_torus1+=transition;

    V.resize(V_torus0.rows()+V_torus1.rows(),V_torus0.cols());
    V<<V_torus0,V_torus1;
    F.resize(F_torus0.rows()+F_torus1.rows(),F_torus0.cols());
    F<<F_torus0,(F_torus1.array()+V_torus0.rows());
}

