cmake_minimum_required(VERSION 2.8)

project (main)


aux_source_directory(./ DIR_SRCS)
add_definitions("-Wall -g")

add_executable(${PROJECT_NAME} ${DIR_SRCS})
target_link_libraries(${PROJECT_NAME})