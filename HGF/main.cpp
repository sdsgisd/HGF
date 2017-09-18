//  main.cpp
//
//  Created by Fang Da on 2014.
//
//  Modified by Sadashige Ishida in 2017.

#include <iostream>
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

#include "Options.h"
#include "Sim.h"
#include "MeshIO.h"

Sim g_sim(false);

struct SimControl
{
    bool wo_visualization;
    
    int win_w;
    int win_h;
    
    bool step;
    bool run;
    bool autoload;
    
    double view_theta;
    double view_alpha;
    double view_dist;
    
    int mouse_x;
    int mouse_y;
    
    bool ldrag;
    int ldrag_start_x;
    int ldrag_start_y;
    bool rdrag;
    int rdrag_start_x;
    int rdrag_start_y;
    Sim::RenderMode render_mode;
    
    int selection_mode;
    
} g_sc;

void renderBitmapString(float x, float y, float z, std::string s)
{
    glColor3f(0, 0, 0);
    glRasterPos3f(x, y, z);
    for (size_t i = 0; i < s.size(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, s[i]);
}

void display()
{
    // object center and zoom
    Vec3d center(0, 0, 0);
    for (size_t i = 0; i < g_sim.get_hgf()->mesh().nv(); i++)
        center += g_sim.get_hgf()->pos(i);
    center /= g_sim.get_hgf()->mesh().nv();
    Vec3d radius(0, 0, 0);
    for (size_t i = 0; i < g_sim.get_hgf()->mesh().nv(); i++)
        for (size_t j = 0; j < 3; j++)
            radius[0] = std::max(radius[0], (g_sim.get_hgf()->pos(i) - center)[0]);
    double min_d = std::max(std::max(radius[0], radius[1]), radius[2]) * 2.2;
    g_sc.view_dist = std::max(min_d, g_sc.view_dist);
    
    glClearColor(1, 1, 1, 1);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, (double)g_sc.win_w / g_sc.win_h, 0.001, 50);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0, -g_sc.view_dist, 0, 0, 0, 0, 0, 0, 1);
    glRotated(g_sc.view_alpha, 1, 0, 0);
    glRotated(g_sc.view_theta, 0, 0, 1);
    
    glBegin(GL_LINES);
    glColor3d(1, 0, 0);     glVertex3d(0, 0, 0);    glVertex3d(2, 0, 0);
    glColor3d(0, 1, 0);     glVertex3d(0, 0, 0);    glVertex3d(0, 2, 0);
    glColor3d(0, 0, 1);     glVertex3d(0, 0, 0);    glVertex3d(0, 0, 2);
    glEnd();
    
    g_sim.render(g_sc.render_mode, Vec2d((double)g_sc.mouse_x / g_sc.win_w * 2 - 1, 1 - (double)g_sc.mouse_y / g_sc.win_h * 2), g_sc.selection_mode);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, g_sc.win_w, 0, g_sc.win_h, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    std::stringstream ss;
    ss << "T = " << g_sim.time();
    std::string s = ss.str();
    renderBitmapString(20, g_sc.win_h - 20, 0, s);
    
    glutSwapBuffers();
}

void idle()
{
    if (g_sc.run || g_sc.step)
    {
        g_sc.step = false;
        g_sim.step();
        std::cout << "Finished step: T = " << g_sim.time() << std::endl;
        g_sim.stepOutput(g_sc.wo_visualization);
        if (g_sim.isFinished()){
            g_sim.get_hgf()->writeObj_FaceLabel_constrainedVertices();
            exit(0);
        }
        
        if (!g_sc.wo_visualization)
            glutPostRedisplay();
    }
    
    if (g_sc.autoload)
    {
        if (!g_sim.load(1))
            exit(0);
        
        g_sim.stepOutput(g_sc.wo_visualization);
        
        if (!g_sc.wo_visualization)
            glutPostRedisplay();
    }
}

void keyboard(unsigned char k, int x, int y)
{
    if (k == 27 || k == 'q' || k == 'Q')
    {
        exit(0);
    } else if (k == ' ')
    {
        g_sc.run = !g_sc.run;
    } else if (k == 's' || k == 'S')
    {
        g_sc.step = true;
    } else if (k == 'm' || k == 'M')
    {
        g_sc.render_mode = (Sim::RenderMode)(((int)g_sc.render_mode + (k == 'm' ? 1 : -1)) % ((int)Sim::RM_COUNT));
        std::cout << "Render mode: " << (int)g_sc.render_mode << std::endl;
        glutPostRedisplay();
    } else if (k == 'v' || k == 'V')
    {
        g_sc.selection_mode = (k == 'v' ? (g_sc.selection_mode | Sim::SM_VERTEX) : (g_sc.selection_mode & ~Sim::SM_VERTEX));
        std::cout << "Mouse cursor selecting" << ((g_sc.selection_mode & Sim::SM_VERTEX) ? " vertices" : "") << ((g_sc.selection_mode & Sim::SM_EDGE) ? " edges" : "") << ((g_sc.selection_mode & Sim::SM_FACE) ? " faces" : "") << "." << std::endl;
        glutPostRedisplay();
    } else if (k == 'e' || k == 'E')
    {
        g_sc.selection_mode = (k == 'e' ? (g_sc.selection_mode | Sim::SM_EDGE) : (g_sc.selection_mode & ~Sim::SM_EDGE));
        std::cout << "Mouse cursor selecting" << ((g_sc.selection_mode & Sim::SM_VERTEX) ? " vertices" : "") << ((g_sc.selection_mode & Sim::SM_EDGE) ? " edges" : "") << ((g_sc.selection_mode & Sim::SM_FACE) ? " faces" : "") << "." << std::endl;
        glutPostRedisplay();
    } else if (k == 'f' || k == 'F')
    {
        g_sc.selection_mode = (k == 'f' ? (g_sc.selection_mode | Sim::SM_FACE) : (g_sc.selection_mode & ~Sim::SM_FACE));
        std::cout << "Mouse cursor selecting" << ((g_sc.selection_mode & Sim::SM_VERTEX) ? " vertices" : "") << ((g_sc.selection_mode & Sim::SM_EDGE) ? " edges" : "") << ((g_sc.selection_mode & Sim::SM_FACE) ? " faces" : "") << "." << std::endl;
        glutPostRedisplay();
    }
    else if (k == 'n' || k == 'N')
    {
      g_sim.showPrimitiveInfo();
    }

    else if (k == 'o' or k=='O' )
    {
        const bool write_imaginary_vertices= k=='o'?false:true;
        g_sim.get_hgf()->writeObj_FaceLabel_constrainedVertices(write_imaginary_vertices);
        g_sim.get_hgf()->write_constrained_mesh("./constrained_mesh.obj");
    }

    else if (k == 'i' || k == 'I')
    {
        g_sim.get_hgf()->test_update_mesh_via_LosTopos();
        
    }
    else if (k == '+')
    {
        //Increasing the volume of the 0-th bubble.
        g_sim.get_hgf()->blowing_bubble0=!g_sim.get_hgf()->blowing_bubble0;
    }
    else if (k == '-')
    {
        //Decreasing the volume of the 0-th bubble.
        g_sim.get_hgf()->absorbing_bubble0=!g_sim.get_hgf()->absorbing_bubble0;
    }

    else if (k == '/')
    {
        g_sim.get_hgf()->do_stepScenes=!g_sim.get_hgf()->do_stepScenes;
    }
    else if (k == 'b')
    {
        //Burst a randomly chosen bubble.
        g_sim.get_hgf()->bursting=true;
        g_sim.step();
        g_sim.get_hgf()->bursting=false;
        
    }
    else if (k == 'r')
    {
        //Move constrained vertices to right.
        g_sim.get_hgf()->move_right=!g_sim.get_hgf()->move_right;
        g_sim.get_hgf()->move_left=false;
    }
    else if (k == 'l')
    {
        //Move constrained vertices to left.
        g_sim.get_hgf()->move_left=!g_sim.get_hgf()->move_left;
        g_sim.get_hgf()->move_right=false;
    }
    else if (k == 'g')
    {
        Scenes::give_large_velocities_to_bubbles(0, &g_sim, g_sim.get_hgf());
    }
    else if (k == 'L')
    {
        g_sim.camera_information=true;
        
    }
    else if (k == 'T')
    {
        //For adaptive time steps.
        static bool fast=true;
        static bool save_mesh_on=g_sim.hgf->save_mesh;
        if(fast){
            g_sim.m_dt=0.001;
            fast=false;
            if(save_mesh_on){
                g_sim.hgf->save_mesh=false;
            }
        }else{
            g_sim.m_dt=0.01;
            fast=true;
            //write_mesh on
            if(save_mesh_on){
                g_sim.hgf->save_mesh=true;
            }
        }
    }
    else if(k == 'Y'){
        //For adaptive time steps.
        g_sim.m_dt=0.003;
        g_sim.hgf->damp=0.998;
        
    }
    
    //FIXME: Loading functions are curretnly unavailable.
    //    else if (k == ']' || k == '}')
    //    {
    //        g_sim.load(k == ']' ? 1 : 10);
    //        //Currently, not available.
    //
    //        glutPostRedisplay();
    //    }
    //    else if (k == '[' || k == '{')
    //    {
    //        g_sim.load(k == '[' ? -1 : -10);
    //        //Currently, not available.
    //
    //
    //        glutPostRedisplay();
    //    }
    //    else if (k == '.' || k == '>')
    //    {
    //        g_sim.load(k == '.' ? 100 : 1000);
    //        //Currently, not available.
    //
    //        glutPostRedisplay();
    //    }
    //    else if (k == 'a' || k == 'A')
    //    {
    //        g_sc.autoload = true;
    //    }

}

void mouse(int b, int s, int x, int y)
{
    if (b == GLUT_LEFT_BUTTON && s == GLUT_DOWN)
    {
        g_sc.ldrag = true;
        g_sc.ldrag_start_x = x;
        g_sc.ldrag_start_y = y;
    } else if (b == GLUT_LEFT_BUTTON && s == GLUT_UP)
    {
        g_sc.ldrag = false;
    } else if (b == GLUT_RIGHT_BUTTON && s == GLUT_DOWN)
    {
        g_sc.rdrag = true;
        g_sc.rdrag_start_x = x;
        g_sc.rdrag_start_y = y;
    } else if (b == GLUT_RIGHT_BUTTON && s == GLUT_UP)
    {
        g_sc.rdrag = false;
    }
    
    glutPostRedisplay();
}

void motion(int x, int y)
{
    if (g_sc.ldrag)
    {
        g_sc.view_theta += (x - g_sc.ldrag_start_x) * 1.0;
        g_sc.view_alpha += (y - g_sc.ldrag_start_y) * 1.0;
        
        g_sc.ldrag_start_x = x;
        g_sc.ldrag_start_y = y;
    }
    if (g_sc.rdrag)
    {
        g_sc.view_dist *= pow(2.0, (y - g_sc.rdrag_start_y) * 0.01);
        
        g_sc.rdrag_start_x = x;
        g_sc.rdrag_start_y = y;
    }
    
    g_sc.mouse_x = x;
    g_sc.mouse_y = y;
    
    glutPostRedisplay();
}

void passiveMotion(int x, int y)
{
    g_sc.mouse_x = x;
    g_sc.mouse_y = y;
    
    glutPostRedisplay();
}

void reshape(int w, int h)
{
    g_sc.win_w = w;
    g_sc.win_h = h;
    
    glViewport(0, 0, w, h);
    
    glutPostRedisplay();
}

void parse_arguments(int argc, char **argv,bool *output_ptr,bool*wo_visualization_ptr,std::string*env_map_path_ptr,std::string*inputdata_dir){
    
    using std::cout;
    using std::endl;
    using std::string;
    using std::stringstream;
    
    for(int i=1;i<argc;++i){
        
        stringstream stream(argv[i]);
        string attribute_name;
        string attribute_value;
        
        getline(stream,attribute_name,'=');
        getline(stream,attribute_value,'=');

        //        if(attribute_value.empty()){
        //            continue;
        //        }
        //        assert(!attribute_value.empty());

        if(attribute_name=="output"){
            *output_ptr=1;
        }
        else if(attribute_name=="wo_visualization"){
            *wo_visualization_ptr=1;
        }
        else if(attribute_name=="env_map_path"){
            *env_map_path_ptr= attribute_value;
        }
        else if(attribute_name=="inputdata_dir"){
            *inputdata_dir= attribute_value;
        }

    }

}

int main(int argc, char * argv[])
{

    std::cout << "[Usage] specify a configuration file as the command line arguments. e.g. doublebubble.txt" << std::endl;
    std::cout<< "[Other options]"<<std::endl;
    std::cout << "output: for outputting images and meshes." << std::endl;
    std::cout << "wo_visualization: Simulating without visualization." << std::endl;
    std::cout << "env_map_path=\"ENV_MAP_PATH\": Specify the directory of the environment map." << std::endl;
    std::cout << "inputdata_dir=\"INPUT_DATA_DIR\": Specify the directory of the input data such as obj or rec." << std::endl;
    std::cout<<std::endl;

    if(argc==1){exit(1);}
    
    // simulation setup
    g_sc.run = false;
    g_sc.step = false;
    g_sc.autoload = false;
    
    g_sc.win_w = 1280;
    g_sc.win_h = 720;
    
    g_sc.view_theta = 0;
    g_sc.view_alpha = 0;
    g_sc.view_dist = 4;
    
    g_sc.mouse_x = 0;
    g_sc.mouse_y = 0;
    
    g_sc.ldrag = false;
    g_sc.ldrag_start_x = 0;
    g_sc.ldrag_start_y = 0;
    g_sc.rdrag = false;
    g_sc.rdrag_start_x = 0;
    g_sc.rdrag_start_y = 0;
    
    g_sc.render_mode = Sim::RM_TRANSPARENT;
    
    g_sc.selection_mode = Sim::SM_VERTEX | Sim::SM_EDGE | Sim::SM_FACE;
    
    if (!g_sc.wo_visualization)
    {
        // glut setup
        glutInit(&argc, argv);
        
        glutInitWindowPosition(0, 0);
        glutInitWindowSize(g_sc.win_w, g_sc.win_h);
        glutCreateWindow("A Hyperbolic Gometric Flow for Evloving Soap Films and Foams");
        
        glutKeyboardFunc(keyboard);
        glutMouseFunc(mouse);
        glutMotionFunc(motion);
        glutPassiveMotionFunc(passiveMotion);
        glutIdleFunc(idle);
        glutDisplayFunc(display);
        glutReshapeFunc(reshape);
    }
    
    bool output=false;
    std::string env_map_path="";
    std::string inputdata_dir="";

    parse_arguments(argc, argv, &output,&(g_sc.wo_visualization),&env_map_path,&inputdata_dir);

    bool success = g_sim.init(argv[1], output, g_sc.wo_visualization,env_map_path,inputdata_dir);
    if (!success)
        return 1;
    
    std::cout << "Initialization complete. Starting the simulation..." << std::endl;
    
    // main loop
    if (g_sc.wo_visualization)
    {
        g_sc.run = true;
        while (true)
            idle();
    } else
    {
        glutMainLoop();
    }
    
}

