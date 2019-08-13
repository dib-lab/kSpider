# ----------------------------------------------------------------------------
# kCluster project
# ----------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.4)
project (kCluster C CXX)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -Ofast -Wall -Wextra")

# ----------------------------------------------------------------------------
# kProcessor Setup
# ----------------------------------------------------------------------------

add_subdirectory(lib/kProcessor)

set(kProcessor_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/lib/kProcessor/include")
include_directories("${PROJECT_SOURCE_DIR}/lib/kProcessor/ThirdParty/MQF/include")

# Add include directories.
include_directories (${kProcessor_INCLUDE_PATH})

# ----------------------------------------------------------------------------
# END kProcessor Setup
# ----------------------------------------------------------------------------

include_directories("${PROJECT_SOURCE_DIR}/lib/argh")

set(PHMAP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib/parallel-hashmap")
include_directories("${PHMAP_INCLUDE_DIRS}")

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++14 -fPIC -fopenmp -W -pedantic -O3 -Ofast")

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

file(GLOB SOURCES
        "${PROJECT_SOURCE_DIR}/src/*"
        "${PROJECT_SOURCE_DIR}/main.cpp"
        )

set(kCluster_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")


include_directories(${kCluster_INCLUDE_DIRS})

add_executable (kCluster ${SOURCES} src/combinations.cpp)

target_link_libraries (kCluster kProcessor)
target_include_directories(kCluster INTERFACE ${PHMAP_INCLUDE_DIRS} ${kProcessor_INCLUDE_PATH})