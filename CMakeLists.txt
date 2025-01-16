cmake_minimum_required(VERSION 3.24)

project(XML_FAERS)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT DEFINED INCLUDE_DIR)
    if(APPLE)
        set(INCLUDE_DIR "/opt/local/include;/opt/homebrew/include") 
    else()
        set(INCLUDE_DIR "/usr/local/include")
    endif()
endif()
include_directories(
    ${INCLUDE_DIR}
)

find_package(pugixml REQUIRED)

add_executable(FAERSParser main.cpp)

target_link_libraries(FAERSParser pugixml)
