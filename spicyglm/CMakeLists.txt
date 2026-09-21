cmake_minimum_required(VERSION 3.18)
project(spicyglm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

find_package(Eigen3 3.4 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  include(FetchContent)
  # download headers only; Eigen's own CMake project is not configured
  FetchContent_Declare(eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    SOURCE_SUBDIR do-not-configure)
  FetchContent_MakeAvailable(eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${eigen_SOURCE_DIR})
endif()

pybind11_add_module(_core
  cpp/src/bindings.cpp
  cpp/src/spatial.cpp
  cpp/src/model.cpp
  cpp/src/diagnostics.cpp
  cpp/src/binomial.cpp)
target_include_directories(_core PRIVATE cpp/include)
target_link_libraries(_core PRIVATE Eigen3::Eigen)

install(TARGETS _core DESTINATION spicyglm)
