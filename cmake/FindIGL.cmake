# - Find IGL
# Find the IGL headers
#
# IGL_INCLUDE_DIR - where to find the Eigen headers
# IGL_FOUND       - True if Eigen is found

if (IGL_INCLUDE_DIR)
  # already in cache, be silent
  set (IGL_FIND_QUIETLY TRUE)
endif (IGL_INCLUDE_DIR)


# find the headers
find_path(
    IGL_INCLUDE_PATH
    igl/cotmatrix.h
    HINTS IGL_INC_DIR ENV IGL_INC_DIR
    DOC "IGL include directory (http://libigl.github.io/libigl/)"
    PATH_SUFFIXES "include"
    PATHS
    /usr/local/include/igl
    /local/include/igl
)


# handle the QUIETLY and REQUIRED arguments and set EIGEN_FOUND to
# TRUE if all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (IGL "IGL (http://libigl.github.io/libigl/) could not be found. Set IGL_INC_DIR to point to the include directory." IGL_INCLUDE_PATH)

if (IGL_FOUND)
  set (IGL_INCLUDE_DIR ${IGL_INCLUDE_PATH})
endif (IGL_FOUND)
