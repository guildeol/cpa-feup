add_executable(original matrix_mult_original.cpp)
target_link_libraries(original papi)

add_executable(matrix_mult matrix_mult.cpp)
target_link_libraries(matrix_mult papi)