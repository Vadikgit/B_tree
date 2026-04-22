cmake_minimum_required (VERSION 3.27)
project(B_tree CXX)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")

#find_package(GoogleTest REQUIRED)

add_library(btree SHARED btree.cpp btprocessable.cpp)


add_executable(main main.cpp)
add_executable(test test/test.cpp)
add_executable(benchmarks benchmark/benchmark.cpp)

target_link_libraries(main btree)
target_link_libraries(test btree gtest gmock pthread)
target_link_libraries(benchmarks btree benchmark pthread)

file(COPY ${CMAKE_SOURCE_DIR}/btree.h ${CMAKE_SOURCE_DIR}/btprocessable.h DESTINATION ${PROJECT_BINARY_DIR}/include)