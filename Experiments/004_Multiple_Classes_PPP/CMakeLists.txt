add_executable(main-004
  ${CMAKE_CURRENT_SOURCE_DIR}/main-004.cpp
)

target_link_libraries(main-004
  uiiitqr
  uiiitsupport
  ${GLOG}
  ${Boost_LIBRARIES}
)

