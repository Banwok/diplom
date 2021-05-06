cmake_minimum_required(VERSION 3.7)

#get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
#string(REPLACE " " "_" ProjectId ${ProjectId})
#project(${ProjectId})

#add_executable(${ProjectId} Dll1.cpp)

#find_package(OpenMP REQUIRED)
#target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

PYTHON_ADD_MODULE(Dll1 Dll1.cpp)
FILE(COPY bostTRy.py DESTINATION .)
ADD_TEST(NAME 01_sie COMMAND ${PYTHON_EXECUTABLE} bostTRy.py)


