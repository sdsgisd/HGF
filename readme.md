# A Hyperbolic Geometric Flow for Evolving Films and Foams

This is an implementation of the paper “A Hyperbolic Geometric Flow for Evolving Films and Foams” Ishida et al., Transactions on Graphics (SIGGRAPH Asia 2017).

[[Project site]][P]  
<a href="https://sadashigeishida.bitbucket.io/hgf/">  <img src="https://sadashigeishida.bitbucket.io/hgf/teaser_representative_image_small.jpg" height="150px"> </a>  
[[Youtube video]][Y]  
<a href="https://www.youtube.com/watch?v=sqywQK7i4L4"><img src="http://i.ytimg.com/vi/sqywQK7i4L4/0.jpg" width="200px"></a>

[Y]:https://www.youtube.com/watch?v=sqywQK7i4L4
[P]:https://sadashigeishida.bitbucket.io/hgf/  
Authors: Sadashige Ishida and Masafumi Yamamoto  
Lisence: MPL-2.0

## Basic Usage
Running the executables without command line arguments will display the usage. Data files are located in the hgf_assets folder.

[KEY SETTING]  
Space: Turn on/off the clock.  
s: Proceed one time step.  

Slash: Turn on/off the scene specific update.  
b: Burst a randomly chosen bubble.  
l: Move constrained vertices to left.  
r: Move constrained vertices to right.   
+: Inflate a bubble.  
-: Deflate a bubble.   

m: Change rendering mode.  
o: Save the state as files containing information of mesh, labels, and constrained vertices.  
O: Save the state as above, but with ghost vertices and faces.  
and etc.

## Dependencies
This program is built by standard procedures using CMAKE (http://www.cmake.org).
The following external libraries are required:   
Eigen (http://eigen.tuxfamily.org)  
LAPACK (http://www.netlib.org/lapack/)  
libigl (http://libigl.github.io/libigl/)  
OpenGL (https://www.opengl.org/)  
GLEW (http://glew.sourceforge.net/) for non-mac OS

### Acknowledgement
This program is built on SoapFilm3D, which is a simulation program accompanying the paper “Double Bubbles Sans Toil and Trouble: Discrete Circulation-Preserving Vortex Sheets for Soap Films and Foams”, Da et al., Transactions on Graphics (SIGGRAPH 2015).



