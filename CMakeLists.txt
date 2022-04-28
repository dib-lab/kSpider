# ----------------------------------------------------------------------------
# kCluster project
# ----------------------------------------------------------------------------

# Minimum cmake version
cmake_minimum_required (VERSION 3.14...3.20)
# Name of project and that it is C/C++.
project (
        kSpider
        LANGUAGES CXX C
        VERSION 2.0.0
)

set(default_build_type "Release")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -lstdc++fs -fPIC -lgomp -lrt -fopenmp -pedantic -O3 -Ofast")
# ----------------------------------------------------------------------------
# kProcessor Setup
# ----------------------------------------------------------------------------

add_subdirectory(lib/kProcessor)

include_directories("${PROJECT_SOURCE_DIR}/lib/kProcessor/include/kProcessor")
include_directories("${PROJECT_SOURCE_DIR}/lib/kProcessor/include/kProcessor/kDataframes")
include_directories("${PROJECT_SOURCE_DIR}/lib/kProcessor/ThirdParty/MQF/include")

set(kProcessor_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/lib/kProcessor/include" "${PROJECT_SOURCE_DIR}/lib/kProcessor/include/kProcessor")
include_directories(${kProcessor_INCLUDE_PATH})

# ----------------------------------------------------------------------------
# END kProcessor Setup
# ----------------------------------------------------------------------------

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/lib/argh")
include_directories("${PROJECT_SOURCE_DIR}/lib/kProcessor/ThirdParty/caches/include")

set(PHMAP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib/kProcessor/ThirdParty/kmerDecoder/lib/parallel-hashmap")
include_directories("${PHMAP_INCLUDE_DIRS}")

add_library(kSpider STATIC src/pairwise.cpp src/index.cpp src/fastx_to_kf.cpp  ${PROJECT_SOURCE_DIR}/lib/kProcessor/include/kProcessor)
set_target_properties(kSpider PROPERTIES POSITION_INDEPENDENT_CODE 1 CXX_STANDARD 17)
target_link_libraries (kSpider kProcessor z)
target_include_directories(kSpider INTERFACE ${PHMAP_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/lib/kProcessor/include/kProcessor)

# add_executable(index_kframes apps/index_kframes.cpp)
# target_link_libraries(index_kframes kSpider kProcessor z)

# add_executable(normalize apps/normalize_pe.cpp)
# target_link_libraries(normalize kSpider kProcessor z)

# add_executable(dump apps/dump_kframe.cpp)
# target_link_libraries(dump kSpider kProcessor z)

# add_executable(remove_percentile apps/remove_percentile.cpp)
# target_link_libraries(remove_percentile kSpider kProcessor z)

# add_executable(remove_singletones apps/remove_singletones.cpp)
# target_link_libraries(remove_singletones kSpider kProcessor z)

# add_executable(dump_count apps/dump_kmer_count.cpp)
# target_link_libraries(dump_count kSpider kProcessor z)

# add_executable(kmc_to_kfs apps/kmc_to_kfs.cpp)
# target_link_libraries(kmc_to_kfs kSpider kProcessor z)

# add_executable(percentile apps/percentile.cpp)
# target_link_libraries(percentile kSpider kProcessor z)