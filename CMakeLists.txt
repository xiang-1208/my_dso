project(my_dso)
cmake_minimum_required(VERSION 2.8)
set(BUILD_TYPE DEBUG)

set(EXECUTABLE_OUTPUT_PATH bin)
set(LIBRARY_OUTPUT_PATH lib)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

message("Current Cmake version is : " ${CMAKE_VERSION})

find_package(LibZip QUIET)


########################
# set compile options  #
########################
#delete the cache to rebuild after changing the options
option(BUILD_WITH_OPENCV "Build with OpenCV3 for Image related Processing" ON)
option(BUILD_WITH_LibZip "Build with LibZip" ON)
option(BUILD_WITH_Eigen "Build with Eigen" ON)
option(BUILD_WITH_Pangolin "Build with Pangolin" ON)

include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/thirdparty/Sophus
  ${PROJECT_SOURCE_DIR}/thirdparty/sse2neon
) 

########################
# find libs & include  #
########################
#OpenCV3
find_package(LibZip QUIET)
if (BUILD_WITH_OPENCV)
    find_package(OpenCV 3 REQUIRED)
    include_directories(${OpenCV_INCLUDE_DIRS})
    message("OPENCV [OK]: " ${OpenCV_INCLUDE_DIRS})
    add_definitions(-DOPENCV_ON)
endif (BUILD_WITH_OPENCV)
#LibZip
# decide if we have ziplib.
if (LIBZIP_LIBRARY)
	message("--- found ziplib (${LIBZIP_VERSION}), compiling with zip capability.")
	add_definitions(-DHAS_ZIPLIB=1)
	include_directories( ${LIBZIP_INCLUDE_DIR_ZIP} ${LIBZIP_INCLUDE_DIR_ZIPCONF} ) 
else()
	message("--- not found ziplib (${LIBZIP_LIBRARY}), compiling without zip capability.")
	set(LIBZIP_LIBRARY "")
endif()
#Eigen
if (BUILD_WITH_Eigen)
    find_package(Eigen3 REQUIRED)
    include_directories(${EIGEN3_INCLUDE_DIR})
    message("Eigen [OK]: " ${OpenCV_INCLUDE_DIRS})
endif (BUILD_WITH_Eigen)
#Pangolin
if (BUILD_WITH_Pangolin)
    find_package(Pangolin 0.2 REQUIRED)
    include_directories(${Pangolin_INCLUDE_DIRS})
    message("Pangolin [OK]: " ${Pangolin_INCLUDE_DIRS})
endif (BUILD_WITH_Pangolin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")


########################
#      link libs       #
########################
if (BUILD_WITH_OPENCV)
    #link opencv (optional)
    set(DEP_LIBS ${DEP_LIBS} ${OpenCV_LIBS})
endif (BUILD_WITH_OPENCV)
if (BUILD_WITH_LibZip)
    #link LibZip (optional)
    set(DEP_LIBS ${DEP_LIBS} ${LIBZIP_LIBRARY})
endif (BUILD_WITH_LibZip)
if (BUILD_WITH_Eigen)
    #link Eigen (optional)
    set(DEP_LIBS ${DEP_LIBS} ${EIGEN3_LIBS})
endif(BUILD_WITH_Eigen)
if (BUILD_WITH_Pangolin)
    #link Eigen (optional)
    set(DEP_LIBS ${DEP_LIBS} ${Pangolin_LIBRARIES})
endif(BUILD_WITH_Pangolin)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin") # OSX
    set(BOOST_THREAD_LIBRARY boost_thread-mt)
else()
    set(BOOST_THREAD_LIBRARY boost_thread)
endif()
set(DEP_LIBS ${DEP_LIBS} boost_system)

########################
#   add executables    #
########################
#my_dso
message("--- compiling my_dso.")
add_library(my_lib src/util/DatasetReader.cpp /
                    src/util/Undistort.cpp
                    src/util/MinimalImage.cpp
                    src/util/globalutil.cpp
                    src/util/ImageAndExposure.cpp
                    src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp
                    src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp
                    src/IOWrapper/Pangolin/PangolinDSOViewer.cpp
                    src/IOWrapper/Pangolin/KeyFrameDisplay.cpp
                    src/FullSystem/HessianBlocks.cpp
                    src/FullSystem/CoarseInitializer.cpp
                    src/FullSystem/PixelSelector2.cpp
                    src/FullSystem/FullSystem.cpp)
add_executable(my_dso src/main.cpp )
target_link_libraries(my_dso my_lib)
target_link_libraries(my_dso boost_thread )
target_link_libraries(my_dso ${DEP_LIBS})

