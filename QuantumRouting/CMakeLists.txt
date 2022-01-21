add_library(uiiitqr SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/network.cpp
)

target_link_libraries(uiiitqr
  uiiitsupport
)
