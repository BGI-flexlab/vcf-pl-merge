cmake_minimum_required(VERSION 3.17)
project(vcfpaste C)

set(CMAKE_C_STANDARD 11)


include_directories(. /Users/huangzhibo/miniconda3/envs/lowpass2/include)
set(CMAKE_LIBRARY_PATH /Users/huangzhibo/miniconda3/envs/lowpass2/lib)

LINK_DIRECTORIES(/Users/huangzhibo/miniconda3/envs/lowpass2/lib)

find_library(hts_LIB libhts.a)

add_executable(vcfpaste main.c vcfpaste.c vcfpaste.h)

target_link_libraries(vcfpaste ${hts_LIB} curl z pthread bz2 lzma)

