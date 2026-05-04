cmake_minimum_required(VERSION 3.28)

project(rtk VERSION 1.0.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(RTK_HEADERS_CONFIGURED
    include/rtk/constants.hpp
)

set(RTK_HEADERS
    include/rtk/rtk.hpp

    include/rtk/core/core.hpp
    include/rtk/core/types.hpp
    include/rtk/core/math.hpp
    include/rtk/core/projection.hpp
    include/rtk/core/vertices.hpp

    include/rtk/util/util.hpp
    include/rtk/util/log.hpp
    include/rtk/util/assert.hpp
    include/rtk/util/exporter.hpp

    include/rtk/prelude/prelude.hpp
    include/rtk/prelude/types.hpp

    include/rtk/model/model.hpp
    include/rtk/model/fif.hpp
    include/rtk/model/pifs.hpp

    include/rtk/dif/dif.hpp
    include/rtk/dif/utils.hpp
    include/rtk/dif/dcf.hpp
    include/rtk/dif/sdcf.hpp

    include/rtk/measures/d2/utils.hpp
    include/rtk/measures/d3/utils.hpp
)

set(RTK_SOURCES
    sources/rtk/rtk.cpp

    sources/rtk/model/pifs.cpp
)

find_package(scion REQUIRED CONFIG)

find_package(Eigen3 REQUIRED)

list(TRANSFORM RTK_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
list(TRANSFORM RTK_HEADERS PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

foreach(HEADER_HPP IN LISTS RTK_HEADERS_CONFIGURED)
    string(REGEX REPLACE "\\.hpp$" ".ipp" HEADER_IPP "${HEADER_HPP}")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${HEADER_IPP}" "${CMAKE_CURRENT_BINARY_DIR}/${HEADER_HPP}" @ONLY)
    list(APPEND RTK_HEADERS "${CMAKE_CURRENT_BINARY_DIR}/${HEADER_HPP}")
endforeach()

add_library(rtk STATIC)
add_library(rtk::rtk ALIAS rtk)
target_sources(rtk PRIVATE ${RTK_SOURCES} PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" "${CMAKE_CURRENT_BINARY_DIR}/include" FILES ${RTK_HEADERS})
target_include_directories(rtk PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(rtk PUBLIC scion::scion Eigen3::Eigen)
target_compile_features(rtk PUBLIC cxx_std_23)
set_target_properties(rtk PROPERTIES CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
target_compile_options(rtk PRIVATE $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra> $<$<CXX_COMPILER_ID:MSVC>:/W4>)
target_compile_definitions(rtk PRIVATE $<$<PLATFORM_ID:Windows>:SCION_PLATFORM_WINDOWS> $<$<PLATFORM_ID:Darwin>:SCION_PLATFORM_MAC> $<$<PLATFORM_ID:Linux>:SCION_PLATFORM_LINUX> $<$<CXX_COMPILER_ID:MSVC>:SCION_COMPILER_MSVC> $<$<CXX_COMPILER_ID:Clang>:SCION_COMPILER_CLANG> $<$<CXX_COMPILER_ID:AppleClang>:SCION_COMPILER_CLANG> $<$<CXX_COMPILER_ID:GNU>:SCION_COMPILER_GCC> $<$<CONFIG:Debug>:SCION_CONFIGURATION_DEBUG> $<$<CONFIG:Release>:SCION_CONFIGURATION_RELEASE> $<$<CONFIG:RelWithDebInfo>:SCION_CONFIGURATION_PROFILE>)

install(TARGETS rtk EXPORT rtk-targets FILE_SET HEADERS)

install(EXPORT rtk-targets NAMESPACE rtk:: DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/rtk")

configure_package_config_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/rtk-config.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/rtk-config.cmake" INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/rtk")
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/rtk-config-version.cmake" VERSION ${PROJECT_VERSION} COMPATIBILITY AnyNewerVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/rtk-config.cmake" "${CMAKE_CURRENT_BINARY_DIR}/rtk-config-version.cmake" DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/rtk")

if(PROJECT_IS_TOP_LEVEL)
    add_subdirectory(tests)
endif()
