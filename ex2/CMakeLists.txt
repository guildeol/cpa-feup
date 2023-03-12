add_executable(sieve_seq sieve_seq.cpp)
target_link_libraries(sieve_seq papi)

target_include_directories(sieve_seq PUBLIC "../papi_helper")
target_link_libraries(sieve_seq "papi_helper" "papi")