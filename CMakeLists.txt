cmake_minimum_required(VERSION 2.8.4)
project(num_method)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-static")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(SOURCE_FILES main.cpp)
add_library(method STATIC method.cpp)

################################BOOST###################################
find_package(Boost COMPONENTS program_options REQUIRED)
if(NOT Boost_FOUND)
    message(SEND_ERROR "Failed to find boost::program_options.")
    return()
else()
    include_directories(${Boost_INCLUDE_DIRS})
endif()

add_executable(num_method ${SOURCE_FILES})
target_link_libraries(num_method method ${Boost_LIBRARIES})