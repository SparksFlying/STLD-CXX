cmake_minimum_required(VERSION 3.20)


# DEBUG
set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb -w")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")
set(CMAKE_CXX_FLAGS "-Werror=narrowing")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
MESSAGE(STATUS ${CMAKE_CXX_STANDARD})

project(STLD)

# find dependencies
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_SOURCE_DIR}/cmake")
find_package(NTL REQUIRED)
find_package(GMP REQUIRED)
find_package(GMPXX REQUIRED)
find_package(Threads REQUIRED)
find_package(GTest REQUIRED)
find_package(OpenMP)
find_package(Boost REQUIRED COMPONENTS log locale)

if(OPENMP_FOUND)
    message(STATUS "OPENMP FOUND")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(Boost_FOUND)
    message( STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}.")
    message( STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}.")
    message( STATUS "Boost_LIB_VERSION = ${Boost_LIB_VERSION}.")
endif()

# sync submodules
include(UpdateSubmodules)
sync_submodules()

file(GLOB_RECURSE CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/*.cpp)

message(STATUS ${CPP_FILES})

add_executable(main ${CPP_FILES} test/main.cpp )

target_include_directories(main PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_include_directories(main PUBLIC "/usr/local/include/rest_rpc/include")
target_include_directories(main PUBLIC ${GTEST_INCLUDE_DIRS})
target_include_directories(main PUBLIC ${Boost_INCLUDE_DIRS})

target_link_libraries(main
                      PUBLIC
                      "/usr/local/lib/libophelib.a"
                      ${LIBNTL_LIBRARIES}
                      ${LIBGMP_LIBRARIES}
                      ${LIBGMPXX_LIBRARIES}
                      Threads::Threads
                      ${GTEST_BOTH_LIBRARIES}
                      ${Boost_LIBRARIES}
                      )