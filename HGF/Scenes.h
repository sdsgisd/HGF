//
//  Scenes.h
//
//
//  Created by Fang Da on 15/1/26.
//
//  Eddited by Sadashige Ishida 2017.

#ifndef __MultiTracker__Scenes__
#define __MultiTracker__Scenes__

#include <iostream>
#include "HGF.h"

class Scenes
{
public:
    // scene-specific initialization
    static HGF * sceneSphere            (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneTet               (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneCube              (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneSheet             (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneBarrel            (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneDoubleBubble      (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneTwoBubbles        (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneTripleJunction    (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneFoamInit          (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);

    static HGF * sceneQuadJunction      (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneConstrainedSphere (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneBubbleWand        (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneTwoRingsPinching  (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);

    static HGF * scenePeanutBubble      (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneStraw             (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneCarousel          (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneOctahedron        (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    
    static HGF * sceneInputData        (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx,const std::string inputdata_dir="");
    static HGF * sceneBrakke       (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx,const std::string inputdata_dir="");
    static HGF * sceneCubicFrame        (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneLattice        (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);
    static HGF * sceneCubeOverFilm        (Sim * sim, std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, std::vector<size_t> & cv, std::vector<Vec3d> & cx);

    static void setPulling (HGF*hgf);

    // scene-specific time evolution
    static void stepSphere             (double dt, Sim * sim, HGF * hgf);
    static void stepTet                (double dt, Sim * sim, HGF * hgf);
    static void stepCube               (double dt, Sim * sim, HGF * hgf);
    static void stepSheet              (double dt, Sim * sim, HGF * hgf);
    static void stepBarrel             (double dt, Sim * sim, HGF * hgf);
    static void stepDoubleBubble       (double dt, Sim * sim, HGF * hgf);
    static void stepTwoBubbles         (double dt, Sim * sim, HGF * hgf);
    static void stepTripleJunction     (double dt, Sim * sim, HGF * hgf);
    static void stepFoamInit           (double dt, Sim * sim, HGF * hgf);
    static void stepQuadJunction       (double dt, Sim * sim, HGF * hgf);
    static void stepConstrainedSphere  (double dt, Sim * sim, HGF * hgf);
    static void stepBubbleWand         (double dt, Sim * sim, HGF * hgf);
    static void stepTwoRingsPinching   (double dt, Sim * sim, HGF * hgf);
    static void stepPeanutBubble       (double dt, Sim * sim, HGF * hgf);
    static void stepStraw              (double dt, Sim * sim, HGF * hgf);
    static void stepCarousel           (double dt, Sim * sim, HGF * hgf);
    static void stepOctahedron         (double dt, Sim * sim, HGF * hgf);
    static void stepCubicFrame        (double dt, Sim * sim, HGF * hgf);
    static void stepCubeOverFilm       (double dt, Sim * sim, HGF * hgf);
    static void stepBrakke             (double dt, Sim * sim, HGF * hgf);

    static void burstBubbles               (double dt, Sim * sim, HGF * hgf);
    static void pullBubbles             (double dt, Sim * sim, HGF * hgf);

    static void set_initial_bubble_velocities (Sim * sim, HGF * hgf);
    static void add_velocities_to_bubbles (double dt, Sim * sim, HGF * hgf);
    static void give_large_velocities_to_bubbles           (double dt, Sim * sim, HGF * hgf);
    static void advect_vertices (double dt, Sim * sim, HGF * hgf);
    
    static void moveLeftOrRight(double dt, Sim * sim, HGF * hgf);
    static void volume_change(double dt, Sim * sim, HGF * hgf);
    
    //Set up constrained vertices for rendering purposes.
    static void wire_frame           (Sim * sim, HGF * hgf);
    static void square_frame           (Sim * sim, HGF * hgf,double frame_out, double height);
    static void double_torus           (Sim * sim, HGF * hgf);
    static double frame_center_x;
    
    static double d_of_rings;//for catenoid.
    static double R_of_rings;//for catenoid.
    static double frame_out;//for cubicframe.

};

#endif /* defined(__MultiTracker__Scenes__) */
