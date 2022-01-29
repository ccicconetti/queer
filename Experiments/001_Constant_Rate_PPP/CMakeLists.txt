add_executable(main-001
  ${CMAKE_CURRENT_SOURCE_DIR}/main-001.cpp
)

target_link_libraries(main-001
  uiiitqr
  uiiitsupport
  ${GLOG}
  ${Boost_LIBRARIES}
)

