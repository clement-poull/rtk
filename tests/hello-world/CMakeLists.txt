cmake_minimum_required(VERSION 3.28)

project(rtk-test-hello-world VERSION 1.0.0 LANGUAGES CXX)

# the following line is commented because the rtk library is already present (this file was included by rtk's CMakeLists.txt)
# in your client program, you need it
# find_package(rtk REQUIRED)

add_executable(rtk-test-hello-world sources/main.cpp)
target_link_libraries(rtk-test-hello-world PUBLIC rtk::rtk)
