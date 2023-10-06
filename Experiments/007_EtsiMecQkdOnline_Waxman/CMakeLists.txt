add_executable(main-007
  ${CMAKE_CURRENT_SOURCE_DIR}/main-007.cpp
)

target_link_libraries(main-007
  uiiitqr
  uiiitsupport
  ${GLOG}
  ${Boost_LIBRARIES}
)

