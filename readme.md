# A Hyperbolic Geometric Flow for Evolving Films and Foams

This is an implementation of the paper “A Hyperbolic Geometric Flow for Evolving Films and Foams” Sadashige Ishida, Masafumi Yamamoto, Ryoichi Ando, and Toshiya Hachisuka, Transactions on Graphics (SIGGRAPH Asia 2017).

[[Project site]][P]  
<a href="https://sadashigeishida.bitbucket.io/hgf/">  <img src="https://sadashigeishida.bitbucket.io/hgf/teaser_representative_image_small.jpg" height="150px"> </a>  
[[Youtube video]][Y]  
<a href="https://www.youtube.com/watch?v=sqywQK7i4L4"><img src="http://i.ytimg.com/vi/sqywQK7i4L4/0.jpg" width="200px"></a>

[Y]:https://www.youtube.com/watch?v=sqywQK7i4L4
[P]:https://sadashigeishida.bitbucket.io/hgf/
This program is built by standard procedures using CMAKE (http://www.cmake.org).
The following external libraries are required:   
Eigen (http://eigen.tuxfamily.org)  
OpenGL and GLUT (http://www.opengl.org/resources/libraries/glut/)  
LAPACK (http://www.netlib.org/lapack/)  
libigl (http://libigl.github.io/libigl/)  

Author: Sadashige Ishida  
Lisence: BSD

## Basic Usage
Running the executables without command line arguments will display the usage. Data files are located in the hgf_assets folder.

[KEY SETTING]  
Space: Turn on / off the clock.  
s: Proceed one time step.  

b: burst a randomly chosen bubble.  
l: Move constrained vertices to left.  
r: Move constrained vertices to right.  

m: Change rendering mode. 
h: Switch flow to heyperbolic mean curvature flow / volume preserving hyperbolic mean curvature flow
shift+w: save the screenshot.  
w: Turn on/off image save mode.  
etc

### Acknowledgement
This program is built on SoapFilm3D, which is a simulation program accompanying the paper “Double Bubbles Sans Toil and Trouble: Discrete Circulation-Preserving Vortex Sheets for Soap Films and Foams”, Fang Da, Christopher Batty, Chris Wojtan, and Eitan Grinspun, Transactions on Graphics (SIGGRAPH 2015).



