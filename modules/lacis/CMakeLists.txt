project(lacis)

add_library(lacis SHARED lacis.cpp matrixmarketio.cpp)

add_executable(lacis_basic_test tests/basictest.cpp)
target_link_libraries(lacis_basic_test lacis)

enable_testing()
add_test(NAME lacis.basicop COMMAND lacis_basic_test)