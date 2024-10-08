cmake_minimum_required (VERSION 3.11.1)

if(APPLE)
      set(CMAKE_C_COMPILER "/usr/local/bin/gcc-7")
      set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-7")
endif()

project(pquant CXX)
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
#     nlohman_json 3.2.0
#     GIT_REPOSITORY https://github.com/nlohmann/json.git
# )
# FetchContent_MakeAvailable(nlohman_json)

include_directories(include)
include_directories(include/func)
include_directories(include/util)
# include_directories(${nlohmann_json_SOURCE_DIR}/include/nlohmann)

file(GLOB_RECURSE SRC_FILES src/*.cpp)
list(FILTER SRC_FILES EXCLUDE REGEX "src/main.cpp")
file(GLOB_RECURSE TEST_FILES test/*.cpp)

add_executable(pquant src/main.cpp ${SRC_FILES})
# target_link_libraries(pquant PRIVATE nlohmann_json::nlohmann_json)