if (LAPACK_LIBRARIES)
  set (LAPACK_FIND_QUIETLY TRUE)
endif (LAPACK_LIBRARIES)

find_library (LAPACK_LIBRARIES
  lapack
  HINTS
  LAPACK_LIB
  $ENV{LAPACK_LIB}
  ${LIB_INSTALL_DIR}
  )

##############
#For including clpack.h

set(LAPACK_SEARCH_PATHS
  ${LAPACK_DIR}
  $ENV{LAPACK_DIR}
  $ENV{CMAKE_PREFIX_PATH}
  ${CMAKE_PREFIX_PATH}
  /usr
  /usr/local
  /usr/local/opt/lapack  ## Mac Homebrew install path
  /opt/LAPACK
  ${CMAKE_OSX_SYSROOT}/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers ##Mac default installed lapack path
)

find_path(LAPACK_INCLUDE_DIR 
             NAMES clapack.h 
             PATHS ${LAPACK_SEARCH_PATHS} 
             PATH_SUFFIXES include)
##############

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACK DEFAULT_MSG LAPACK_LIBRARIES)

mark_as_advanced (LAPACK_LIBRARIES)




if (LAPACK_FOUND)
  set (LAPACK_INCLUDE_DIR ${LAPACK_INCLUDE_PATH})

endif (LAPACK_FOUND)


