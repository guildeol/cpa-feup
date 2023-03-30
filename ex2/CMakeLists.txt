add_executable(sieve_seq sieve_seq.cpp)
target_include_directories(sieve_seq PUBLIC "../papi_helper")
target_link_libraries(sieve_seq "papi_helper" "papi")

add_executable(sieve_parallel sieve_parallel.cpp)
target_include_directories(sieve_parallel PUBLIC "../papi_helper")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(sieve_parallel PUBLIC "papi_helper" "papi" OpenMP::OpenMP_CXX)    
endif()

add_executable(sieve_parallel_tasks sieve_parallel_tasks.cpp)
target_include_directories(sieve_parallel_tasks PUBLIC "../papi_helper")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(sieve_parallel_tasks PUBLIC "papi_helper" "papi" OpenMP::OpenMP_CXX)    
endif()