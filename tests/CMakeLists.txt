function(add_test filename)
  set(target test_${filename})
  add_executable(${target} ${filename}.cpp)
  target_link_libraries(${target} PRIVATE GTest::gtest GTest::gtest_main qrot)
  set_target_properties(
    ${target}
    PROPERTIES CXX_STANDARD 20
               CXX_STANDARD_REQUIRED ON
               CXX_EXTENSIONS OFF)
  gtest_discover_tests(${target})
endfunction()
add_test(decomposition)
add_test(diophantine)
add_test(gate)
add_test(geometry)
add_test(grid_solver)
add_test(matrix)
add_test(number)
add_test(parser)
