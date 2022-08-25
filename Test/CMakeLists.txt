if(NOT skip_google_test)

  # Download and unpack googletest at configure time
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt.in googletest-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/googletest-download )
  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/googletest-download )
  if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif()

  # Prevent overriding the parent project's compiler/linker
  # settings on Windows
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  # Add googletest directly to our build. This defines
  # the gtest and gtest_main targets.
  add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                   ${CMAKE_BINARY_DIR}/googletest-build
                   EXCLUDE_FROM_ALL)

endif()

include_directories(${CMAKE_BINARY_DIR}/googletest-src/googletest/include)

include(GoogleTest)

add_library(testqrservice STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/Details/examplenetwork.cpp
)

set(LIBS
  testqrservice
  uiiitqr
  uiiitsupport
  gtest
  ${Boost_LIBRARIES}
)

add_executable(testcapacitynetwork testmain.cpp testcapacitynetwork.cpp)
target_link_libraries(testcapacitynetwork ${LIBS})
gtest_discover_tests(testcapacitynetwork)

add_executable(testgraphml testmain.cpp testgraphml.cpp)
target_link_libraries(testgraphml ${LIBS})
gtest_discover_tests(testgraphml)

add_executable(testqrutils testmain.cpp testqrutils.cpp)
target_link_libraries(testqrutils ${LIBS})
gtest_discover_tests(testqrutils)

add_executable(testpoissonpointprocess testmain.cpp testpoissonpointprocess.cpp)
target_link_libraries(testpoissonpointprocess ${LIBS})
gtest_discover_tests(testpoissonpointprocess)

add_executable(testyen testmain.cpp testyen.cpp)
target_link_libraries(testyen ${LIBS})
gtest_discover_tests(testyen)

