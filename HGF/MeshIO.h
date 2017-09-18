//
//  MeshIO.h
//
//  Christopher Batty, Fang Da 2014
//
//

#ifndef __MeshIO__
#define __MeshIO__

#include <iostream>
#include <vector>
#include "surftrack.h"
#include "HGF.h"

class MeshIO
{
public:
    static bool save(HGF & hgf, const std::string & filename, bool binary = true);
    static bool load(HGF & hgf, const std::string & filename, bool binary = true);
    
    static bool loadIntoRaw(std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, const std::string & filename, bool binary = true);
    
    static bool saveOBJ(HGF & hgf, const std::string & filename,const bool with_imaginary_vertices=true, const bool with_normal=true, const bool change_y_z=false);

};

#endif /* defined(__MeshIO__) */
