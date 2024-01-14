cmake_minimum_required(VERSION 3.16)

project(qrot LANGUAGES CXX)
set(PROJECT_NAMESPACE qrot)

option(QROT_VERBOSE "Print debug info" OFF)

include(cmake/deps.cmake)
add_subdirectory(src)

add_subdirectory(examples)

include(GoogleTest)
enable_testing()
add_subdirectory(tests)
