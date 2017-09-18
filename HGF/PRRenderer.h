//
//  PRRenderer.h
//  
//
//  Created by Fang Da on 1/15/15.
//
//

#ifndef __MultiTracker__PRRenderer__
#define __MultiTracker__PRRenderer__

#include <iostream>
#include "HGF.h"
#include <fstream>
#include <sstream>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#include "Shader.h"

class PRRenderer
{
public:
    PRRenderer(HGF * hgf,const std::string env_map_path="");
    std::string env_map_path;
    
public:
    void render();
    
protected:
    LosTopos::NonDestructiveTriMesh & mesh() { return hgf->mesh(); }

    bool create_cube_map(const std::string & name, GLuint & tex);
    bool load_cube_map_side(GLuint tex, GLenum side, const std::string & filename);
    
protected:
    HGF * hgf;
    Shader m_shader_bubble;
    Shader m_shader_env;
    GLuint m_tex_env;
};

#endif /* defined(__MultiTracker__PRRenderer__) */
