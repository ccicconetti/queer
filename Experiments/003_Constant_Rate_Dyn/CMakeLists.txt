add_executable(main-003
  ${CMAKE_CURRENT_SOURCE_DIR}/main-003.cpp
)

target_link_libraries(main-003
  uiiitqr
  uiiitsupport
  ${GLOG}
  ${Boost_LIBRARIES}
)

