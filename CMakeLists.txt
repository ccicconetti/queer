cmake_policy(SET CMP0048 NEW)

project(quantum-routing VERSION 0.0.1)
cmake_minimum_required(VERSION 3.2)

set (CMAKE_CXX_STANDARD 17)

# Cmake find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

if(APPLE)
  find_package(glog REQUIRED)
  set(GLOG glog::glog)
else()
  set(GLOG glog)
endif()
set(Boost_USE_MULTITHREADED TRUE)
find_package(Boost
  REQUIRED
  COMPONENTS chrono filesystem graph program_options system thread
)

set(DISABLE_WARNINGS "-Wno-missing-field-initializers -Wno-unused-parameter -Wnon-virtual-dtor")
set(COMPILER_COMMON_FLAGS "-Wall -Wextra -Werror -DGTEST_CREATE_SHARED_LIBRARY=1 -DGTEST_LINKED_AS_SHARED_LIBRARY=1 -fPIC ${DISABLE_WARNINGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${COMPILER_COMMON_FLAGS} -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${COMPILER_COMMON_FLAGS} -O2 -DNDEBUG")

MESSAGE("============CONFIGURATION SUMMARY================")
MESSAGE("")
MESSAGE("CMAKE_SOURCE_DIR:         ${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE("COMPILER:                 ${CMAKE_CXX_COMPILER}")
MESSAGE("COMPILER FLAGS DEBUG:     ${CMAKE_CXX_FLAGS_DEBUG}")
MESSAGE("COMPILER FLAGS RELEASE:   ${CMAKE_CXX_FLAGS_RELEASE}")
MESSAGE("CMAKE_BUILD_TYPE:         ${CMAKE_BUILD_TYPE}")

# header of local libraries
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# local of header of external repos
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/support)

# external dependencies
include_directories(${Boost_INCLUDE_DIRS})

# executables
# add_subdirectory(Executables)

# local libraries
add_subdirectory(QuantumRouting)

# add experiments
add_subdirectory(Experiments)

# external repos
add_subdirectory(support/Support)

# unit tests
string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_LOWER)
if (${CMAKE_BUILD_TYPE_LOWER} STREQUAL "debug")
  add_subdirectory(Test)
  set(skip_google_test TRUE)
  add_subdirectory(support/Test)
endif()
