cmake_minimum_required (VERSION 3.6)
project (ProgParallel)
find_package(Threads REQUIRED)
find_package(SDL2 REQUIRED)
find_package(benchmark REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_RELEASE_FLAGS "${CMAKE_CXX_RELEASE_FLAGS} -O3 -march=native")

include_directories(lib/include)

add_library(renderer src/render.hpp src/render.cpp)
add_executable(view src/view.cpp)
target_include_directories(view PRIVATE ${SDL2_INCLUDE_DIRS})
target_link_libraries(view renderer SDL2::SDL2 Threads::Threads)


add_executable(bench src/bench.cpp)
target_link_libraries(bench renderer benchmark::benchmark)
