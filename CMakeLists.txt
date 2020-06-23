cmake_minimum_required(VERSION 3.16)
project(connected_components_labeling)

set(CMAKE_CXX_STANDARD 17)

# Where to find CMake modules and OpenCV
#set(OpenCV_DIR "C:\\opencv\\mingw-build\\install")
#find_package(OpenCV REQUIRED)

if (MSVC)
    message("-- Configuring Microsoft Visual Studio Compiler for max performance")
    set(CMAKE_CXX_FLAGS "/Ob2 /O2 /GR /MDd /Zi /Qspectre /DWIN32 /D_WINDOWS /W3 /EHsc" CACHE INTERNAL "")
    macro(RemoveDebugCXXFlag flag)
        string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
        string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_DEBUG_INIT "${CMAKE_CXX_FLAGS_DEBUG_INIT}")
    endmacro()

    RemoveDebugCXXFlag("/RTC1")
    RemoveDebugCXXFlag("/Ob0")
    RemoveDebugCXXFlag("/Od")
endif()

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#message(${CMAKE_CXX_FLAGS})
#message(${CMAKE_EXE_LINKER_FLAGS})

add_executable(connected_components_labeling main.cpp include/stb_image.h include/image.h include/buffer.h include/cluster_perc.h)

if (GNU|CLANG)
    message("-- Configuring Clang or GNU Compiler for max performance")
    target_link_libraries(connected_components_labeling -ansi -pedantic -Wall -Wextra -O3)
endif()

# add libs you need
#set(OpenCV_LIBS opencv_core opencv_imgproc opencv_highgui opencv_imgcodecs)
#target_link_libraries(connected_components_labeling ${OpenCV_LIBS})
#include_directories(${OpenCV_INCLUDE_DIRS})

