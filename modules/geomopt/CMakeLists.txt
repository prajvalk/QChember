project(newscf_geomopt)

add_library(newscf_geomopt pso.cpp)
target_link_libraries(newscf_geomopt)

add_executable(geomopt_rastrigin tests/rastrigin.cpp)
add_executable(geomopt_lj tests/lj_clusters.cpp)
target_link_libraries(geomopt_rastrigin newscf_geomopt)
target_link_libraries(geomopt_lj newscf_geomopt)