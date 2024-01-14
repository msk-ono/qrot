add_library(
  qrot
  qrot/decomposition.cpp
  qrot/diophantine.cpp
  qrot/gate.cpp
  qrot/geometry.cpp
  qrot/grid_solver.cpp
  qrot/matrix.cpp
  qrot/number.cpp
  qrot/parser.cpp)
target_include_directories(qrot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(qrot PUBLIC Boost::boost)
# target_link_libraries(qrot PUBLIC Boost::boost OpenMP::OpenMP_CXX)
set_target_properties(
  qrot
  PROPERTIES CXX_STANDARD 20
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF)
target_compile_options(qrot PRIVATE -Wall -Wextra)
if(QROT_VERBOSE)
  target_compile_definitions(qrot PUBLIC QROT_VERBOSE)
endif()
