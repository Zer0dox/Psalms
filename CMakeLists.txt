cmake_minimum_required(VERSION 3.26)
project(Psalms C)

set(CMAKE_C_STANDARD 11)

add_executable(Psalms main.c
        durac.c
        durac.h
)
