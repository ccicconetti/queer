add_library(uiiitqr STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/capacitynetwork.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/esnetwork.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mecqkdnetwork.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mecqkdonlinenetwork.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mecqkdworkload.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/peerassignment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/poissonpointprocess.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/qrutils.cpp
)

target_link_libraries(uiiitqr
  hungarian
  uiiitsupport
)
