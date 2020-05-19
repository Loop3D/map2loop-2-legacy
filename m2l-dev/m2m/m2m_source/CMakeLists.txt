# Author: Vitaliy Ogarko, vogarko@gmail.com
cmake_minimum_required (VERSION 2.8)
project (Map2model-cpp)

option(DISABLE_LTO OFF)

SET (SRCS src/map2model.cpp
          src/parameters_utils.cpp
          src/parameters_reader.cpp
          src/converter_utils.cpp
          src/intersector.cpp
          src/data_reader.cpp
          src/converter.cpp
          src/graph.cpp
          src/graph_utils.cpp
          src/local_topology_analyzer.cpp
          src/topology_analyzer.cpp
          src/lib/clipper.cpp
          src/lib/svgbuilder.cpp)

SET (SRCS_INCLUDE src
                  src/lib)

add_executable(map2model ${SRCS})
target_include_directories(map2model PUBLIC ${SRCS_INCLUDE})

########################################################################################

if (CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=gnu++0x)
endif()

# Adding compiler flags.
if ( CMAKE_COMPILER_IS_GNUCC )
    SET (GCC_FLAGS -m64 -mfpmath=sse -O3 -ffast-math -march=native -funroll-loops -Wall -Wno-unknown-pragmas -Wextra)
    IF (NOT DISABLE_LTO)
        SET (GCC_FLAGS ${GCC_FLAGS} -flto)
    endif()

#    For printing vectorization info.
#    SET (GCC_FLAGS -m64 -mfpmath=sse -O3 -ffast-math -march=native -funroll-loops -ftree-vectorize -fopt-info-vec-optimized)
    target_compile_options(map2model PRIVATE ${GCC_FLAGS})
endif()
if ( MSVC )
    target_compile_options(map2model PRIVATE /W4)
endif()
