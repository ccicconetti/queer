add_executable(main-002
  ${CMAKE_CURRENT_SOURCE_DIR}/main-002.cpp
)

target_link_libraries(main-002
  uiiitqr
  uiiitsupport
  ${GLOG}
  ${Boost_LIBRARIES}
)

