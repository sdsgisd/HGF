append_files (Headers "h" .)
append_files (Sources "cpp" .)

include_directories (${EIGEN_INCLUDE_DIR})

add_executable (HGF ${Headers} ${Sources})
target_link_libraries (HGF LosTopos)