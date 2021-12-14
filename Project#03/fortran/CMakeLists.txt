# References: 
# https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/languages/fortran/ForFortranExample
# https://fortranwiki.org/fortran/show/Build+tools

cmake_minimum_required (VERSION 3.12)
project (els)
enable_language (Fortran)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set (dialect "-ffree-form -std=f2008 -fimplicit-none")
    set (extras "-fcheck=all -ggdb")
endif ()
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set (dialect "-stand f08 -free -implicitnone")
    set (extras "-check bounds")
endif ()

set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${extras}")
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

file (GLOB_RECURSE sources src/*.f90 src/*.F90)
add_executable(els.x ${sources})
target_link_libraries (els.x /opt/OpenBLAS/lib/libopenblas.a -lpthread)
