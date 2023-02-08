add_executable(main-006
  ${CMAKE_CURRENT_SOURCE_DIR}/main-006.cpp
)

target_link_libraries(main-006
  uiiitqr
  uiiitsupport
  ${GLOG}
  ${Boost_LIBRARIES}
)

