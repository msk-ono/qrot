add_executable(gridsynth_cpp gridsynth.cpp)
target_link_libraries(gridsynth_cpp PUBLIC Boost::boost Boost::program_options
                                           qrot)
set_target_properties(
  gridsynth_cpp
  PROPERTIES CXX_STANDARD 20
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF)
target_compile_options(gridsynth_cpp PRIVATE -Wall -Wextra)
if(QROT_VERBOSE)
  target_compile_definitions(gridsynth_cpp PUBLIC QROT_VERBOSE)
endif()
