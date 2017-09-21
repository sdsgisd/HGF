cmake_minimum_required (VERSION 2.6)

cmake_policy(SET CMP0042 OLD)

#Use gcc/g++-5
#set(CMAKE_CXX_COMPILER /usr/local/bin/g++-5)
#set(CMAKE_C_COMPILER /usr/local/bin/gcc-5)

#Use C++14.
#set(CMAKE_CXX_FLAGS “${CMAKE_CXX_FLAGS} -std=c++14 -Wall”)
set(CMAKE_CXX_STANDARD 14)


project (HGF)
set (MULTITRACKER_MAJOR_VERSION "1")
set (MULTITRACKER_MINOR_VERSION "0")
set (MULTITRACKER_BUGFIX_VERSION "0")
set (MULTITRACKER_VERSION_NUMBER "${MULTITRACKER_MAJOR_VERSION}.${MULTITRACKER_MINOR_VERSION}.${MULTITRACKER_BUGFIX_VERSION}")
set (MULTITRACKER_VERSION "${MULTITRACKER_VERSION_NUMBER}")
add_definitions (-DMULTITRACKER_VERSION="${MULTITRACKER_VERSION}")

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# add directory with macros
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
# add cmake file
include (common)

# set options

if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions (-DDEBUG)
endif (CMAKE_BUILD_TYPE MATCHES Debug)


#Eigen library is required
find_package (Eigen REQUIRED)
include_directories (${EIGEN_INCLUDE_DIR})

option (USE_PNG "Builds in support for exporting to png" ON)
if (USE_PNG)
  find_package (PNG)
  if (PNG_FOUND)
    add_definitions (-DHAVE_PNG)
    add_definitions (${PNG_DEFINITIONS})
    include_directories (${PNG_INCLUDE_DIR})
    set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${PNG_LIBRARIES})
  endif (PNG_FOUND)
endif (USE_PNG)

include_directories(/usr/local/include)
link_directories(/usr/X11/lib)

#For Lapack
find_package (LAPACK REQUIRED)
include_directories(${LAPACK_INCLUDE_DIR})

add_definitions (-DHAVE_LAPACK)

find_package (BLAS REQUIRED)
add_definitions (-DHAVE_BLAS)

# Rendering dependencies
find_package (OpenGL REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS} )

find_package (GLEW REQUIRED)
include_directories( ${GLEW_INCLUDE_DIRS} )

# LIBIGL
find_package (IGL REQUIRED)
include_directories( ${IGL_INCLUDE_DIR} )

find_package (Threads)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/LosTopos3D)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common/newsparse)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common/tunicate)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory (LosTopos)
add_subdirectory (HGF)

