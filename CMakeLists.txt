cmake_minimum_required(VERSION 3.14)
project(archimede_raytracer CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(archimede src/main.cpp src/simulation.cpp)
target_compile_options(archimede PRIVATE -Wall -Wextra)

# OpenMP est optionnel : il accélère les balayages en utilisant tous les coeurs.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(archimede PRIVATE OpenMP::OpenMP_CXX)
endif()
