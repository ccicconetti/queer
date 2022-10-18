add_executable(main-005
  ${CMAKE_CURRENT_SOURCE_DIR}/main-005.cpp
)

target_link_libraries(main-005
  uiiitqr
  uiiitsupport
  ${GLOG}
  ${Boost_LIBRARIES}
)

