cmake_minimum_required (VERSION 3.11.1)

### To use gcc/g++ on a Macintosh, you must set the Compilers
### here, not inside the project
##if(APPLE)
##       set(CMAKE_C_COMPILER "/usr/local/bin/gcc-7")
##       set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-7")
##endif()
### TODO: for now, we use CLang for Mac
###
### In order to create OpenFHE's static libraries you should enable
### the BUILD_STATIC option. For that, you run "cmake .. -DBUILD_STATIC=ON".
### After having your link completed you will find static libs
### with the suffix "_static" in ./build/libs/.
### Examples: OPENFHEpke_static.a, OPENFHEcore_static.a, etc.
### After you run "make install" in your build directory, you can build your custom application.
### If you need your application to be linked statically, then run "cmake .. -DBUILD_STATIC=ON"

project(demo CXX)
set(CMAKE_CXX_STANDARD 17)
option( BUILD_STATIC "Set to ON to include static versions of the library" OFF)

find_package(OpenFHE CONFIG REQUIRED)
if (OpenFHE_FOUND)
    message(STATUS "FOUND PACKAGE OpenFHE")
    message(STATUS "OpenFHE Version: ${BASE_OPENFHE_VERSION}")
    message(STATUS "OpenFHE installed as shared libraries: ${OpenFHE_SHARED}")
    message(STATUS "OpenFHE include files location: ${OpenFHE_INCLUDE}")
    message(STATUS "OpenFHE lib files location: ${OpenFHE_LIBDIR}")
    message(STATUS "OpenFHE Native Backend size: ${OpenFHE_NATIVE_SIZE}")
else()
    message(FATAL_ERROR "PACKAGE OpenFHE NOT FOUND")
endif ()

set( CMAKE_CXX_FLAGS ${OpenFHE_CXX_FLAGS} )

include_directories( ${OPENMP_INCLUDES} )
include_directories( ${OpenFHE_INCLUDE} )
include_directories( ${OpenFHE_INCLUDE}/third-party/include )
include_directories( ${OpenFHE_INCLUDE}/core )
include_directories( ${OpenFHE_INCLUDE}/pke )
include_directories( ${OpenFHE_INCLUDE}/binfhe )
### add directories for other OpenFHE modules as needed for your project

include_directories(${googletest_SOURCE_DIR}/googletest/include)

link_directories( ${OpenFHE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
if(BUILD_STATIC)
    set( CMAKE_EXE_LINKER_FLAGS "${OpenFHE_EXE_LINKER_FLAGS} -static")
    link_libraries( ${OpenFHE_STATIC_LIBRARIES} )
else()
    set( CMAKE_EXE_LINKER_FLAGS ${OpenFHE_EXE_LINKER_FLAGS} )
    link_libraries( ${OpenFHE_SHARED_LIBRARIES} )
endif()

# include(FetchContent)
# FetchContent_Declare(
#     googletest
#     GIT_REPOSITORY https://github.com/google/googletest.git
#     GIT_TAG        release-1.11.0
# )
# FetchContent_MakeAvailable(googletest)

# # Fetch and populate Google Benchmark
# FetchContent_Declare(
#     googlebenchmark
#     GIT_REPOSITORY https://github.com/google/benchmark.git
#     GIT_TAG v1.5.4
# )
# FetchContent_Populate(googlebenchmark)

# add_subdirectory(${googlebenchmark_SOURCE_DIR} ${googlebenchmark_BINARY_DIR})

# include_directories(${googletest_SOURCE_DIR}/googletest/include)
# include_directories(${googletest_SOURCE_DIR}/googlemock/include)
# include_directories(${googlebenchmark_SOURCE_DIR}/include)

include(FetchContent)
FetchContent_Declare(
    nlohman_json 3.2.0
    GIT_REPOSITORY https://github.com/nlohmann/json.git
)
FetchContent_MakeAvailable(nlohman_json)

include_directories(include)
include_directories(include/func)
include_directories(include/util)
include_directories(${nlohmann_json_SOURCE_DIR}/include/nlohmann)

file(GLOB_RECURSE SRC_FILES src/*.cpp)
list(FILTER SRC_FILES EXCLUDE REGEX "src/main.cpp")
file(GLOB_RECURSE TEST_FILES test/*.cpp)

add_executable(pquant src/main.cpp ${SRC_FILES})
target_link_libraries(pquant PRIVATE nlohmann_json::nlohmann_json)

# add_executable(test_all ${TEST_FILES} ${SRC_FILES})
# target_link_libraries(test_all gtest gtest_main)

# set(benchmark_ROOT_DIR "" CACHE PATH "Folder containing benchmark")

# find_path(benchmark_INCLUDE_DIR "benchmark/benchmark.h"
#   PATHS ${benchmark_ROOT_DIR}
#   PATH_SUFFIXES include
#   NO_DEFAULT_PATH)
# find_path(benchmark_INCLUDE_DIR "benchmark/benchmark.h")

# find_library(benchmark_LIBRARY NAMES "benchmark"
#   PATHS ${benchmark_ROOT_DIR}
#   PATH_SUFFIXES lib lib64
#   NO_DEFAULT_PATH)
# find_library(benchmark_LIBRARY NAMES "benchmark")

# include(FindPackageHandleStandardArgs)
# # handle the QUIETLY and REQUIRED arguments and set benchmark_FOUND to TRUE
# # if all listed variables are TRUE
# find_package_handle_standard_args(benchmark FOUND_VAR benchmark_FOUND
#   REQUIRED_VARS benchmark_LIBRARY
#   benchmark_INCLUDE_DIR)

# if(benchmark_FOUND)
#   set(benchmark_LIBRARIES ${benchmark_LIBRARY})
#   set(benchmark_INCLUDE_DIRS ${benchmark_INCLUDE_DIR})
# endif()

# mark_as_advanced(benchmark_INCLUDE_DIR benchmark_LIBRARY)

# add_executable(benchmark benchmark.cpp ${SRC_FILES})
# include_directories(${benchmark_INCLUDE_DIR})
# target_link_libraries(benchmark ${benchmark_LIBRARY})