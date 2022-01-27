add_library(uiiitqr SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/network.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/poissonpointprocess.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/qrutils.cpp
)

target_link_libraries(uiiitqr
  uiiitsupport
)
