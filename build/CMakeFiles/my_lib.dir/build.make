# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/xiang/Project/my_dso

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/xiang/Project/my_dso/build

# Include any dependencies generated for this target.
include CMakeFiles/my_lib.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/my_lib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/my_lib.dir/flags.make

CMakeFiles/my_lib.dir/src/util/DatasetReader.cpp.o: CMakeFiles/my_lib.dir/flags.make
CMakeFiles/my_lib.dir/src/util/DatasetReader.cpp.o: ../src/util/DatasetReader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiang/Project/my_dso/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/my_lib.dir/src/util/DatasetReader.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_lib.dir/src/util/DatasetReader.cpp.o -c /home/xiang/Project/my_dso/src/util/DatasetReader.cpp

CMakeFiles/my_lib.dir/src/util/DatasetReader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_lib.dir/src/util/DatasetReader.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiang/Project/my_dso/src/util/DatasetReader.cpp > CMakeFiles/my_lib.dir/src/util/DatasetReader.cpp.i

CMakeFiles/my_lib.dir/src/util/DatasetReader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_lib.dir/src/util/DatasetReader.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiang/Project/my_dso/src/util/DatasetReader.cpp -o CMakeFiles/my_lib.dir/src/util/DatasetReader.cpp.s

CMakeFiles/my_lib.dir/src/util/Undistort.cpp.o: CMakeFiles/my_lib.dir/flags.make
CMakeFiles/my_lib.dir/src/util/Undistort.cpp.o: ../src/util/Undistort.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiang/Project/my_dso/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/my_lib.dir/src/util/Undistort.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_lib.dir/src/util/Undistort.cpp.o -c /home/xiang/Project/my_dso/src/util/Undistort.cpp

CMakeFiles/my_lib.dir/src/util/Undistort.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_lib.dir/src/util/Undistort.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiang/Project/my_dso/src/util/Undistort.cpp > CMakeFiles/my_lib.dir/src/util/Undistort.cpp.i

CMakeFiles/my_lib.dir/src/util/Undistort.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_lib.dir/src/util/Undistort.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiang/Project/my_dso/src/util/Undistort.cpp -o CMakeFiles/my_lib.dir/src/util/Undistort.cpp.s

CMakeFiles/my_lib.dir/src/util/MinimalImage.cpp.o: CMakeFiles/my_lib.dir/flags.make
CMakeFiles/my_lib.dir/src/util/MinimalImage.cpp.o: ../src/util/MinimalImage.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiang/Project/my_dso/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/my_lib.dir/src/util/MinimalImage.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_lib.dir/src/util/MinimalImage.cpp.o -c /home/xiang/Project/my_dso/src/util/MinimalImage.cpp

CMakeFiles/my_lib.dir/src/util/MinimalImage.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_lib.dir/src/util/MinimalImage.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiang/Project/my_dso/src/util/MinimalImage.cpp > CMakeFiles/my_lib.dir/src/util/MinimalImage.cpp.i

CMakeFiles/my_lib.dir/src/util/MinimalImage.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_lib.dir/src/util/MinimalImage.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiang/Project/my_dso/src/util/MinimalImage.cpp -o CMakeFiles/my_lib.dir/src/util/MinimalImage.cpp.s

CMakeFiles/my_lib.dir/src/util/globalutil.cpp.o: CMakeFiles/my_lib.dir/flags.make
CMakeFiles/my_lib.dir/src/util/globalutil.cpp.o: ../src/util/globalutil.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiang/Project/my_dso/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/my_lib.dir/src/util/globalutil.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_lib.dir/src/util/globalutil.cpp.o -c /home/xiang/Project/my_dso/src/util/globalutil.cpp

CMakeFiles/my_lib.dir/src/util/globalutil.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_lib.dir/src/util/globalutil.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiang/Project/my_dso/src/util/globalutil.cpp > CMakeFiles/my_lib.dir/src/util/globalutil.cpp.i

CMakeFiles/my_lib.dir/src/util/globalutil.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_lib.dir/src/util/globalutil.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiang/Project/my_dso/src/util/globalutil.cpp -o CMakeFiles/my_lib.dir/src/util/globalutil.cpp.s

CMakeFiles/my_lib.dir/src/util/ImageAndExposure.cpp.o: CMakeFiles/my_lib.dir/flags.make
CMakeFiles/my_lib.dir/src/util/ImageAndExposure.cpp.o: ../src/util/ImageAndExposure.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiang/Project/my_dso/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/my_lib.dir/src/util/ImageAndExposure.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_lib.dir/src/util/ImageAndExposure.cpp.o -c /home/xiang/Project/my_dso/src/util/ImageAndExposure.cpp

CMakeFiles/my_lib.dir/src/util/ImageAndExposure.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_lib.dir/src/util/ImageAndExposure.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiang/Project/my_dso/src/util/ImageAndExposure.cpp > CMakeFiles/my_lib.dir/src/util/ImageAndExposure.cpp.i

CMakeFiles/my_lib.dir/src/util/ImageAndExposure.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_lib.dir/src/util/ImageAndExposure.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiang/Project/my_dso/src/util/ImageAndExposure.cpp -o CMakeFiles/my_lib.dir/src/util/ImageAndExposure.cpp.s

CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp.o: CMakeFiles/my_lib.dir/flags.make
CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp.o: ../src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiang/Project/my_dso/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp.o -c /home/xiang/Project/my_dso/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp

CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiang/Project/my_dso/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp > CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp.i

CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiang/Project/my_dso/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp -o CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp.s

CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp.o: CMakeFiles/my_lib.dir/flags.make
CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp.o: ../src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiang/Project/my_dso/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp.o -c /home/xiang/Project/my_dso/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp

CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiang/Project/my_dso/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp > CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp.i

CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiang/Project/my_dso/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp -o CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp.s

CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp.o: CMakeFiles/my_lib.dir/flags.make
CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp.o: ../src/IOWrapper/Pangolin/PangolinDSOViewer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiang/Project/my_dso/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp.o -c /home/xiang/Project/my_dso/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp

CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiang/Project/my_dso/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp > CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp.i

CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiang/Project/my_dso/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp -o CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp.s

CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp.o: CMakeFiles/my_lib.dir/flags.make
CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp.o: ../src/IOWrapper/Pangolin/KeyFrameDisplay.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiang/Project/my_dso/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp.o -c /home/xiang/Project/my_dso/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp

CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiang/Project/my_dso/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp > CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp.i

CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiang/Project/my_dso/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp -o CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp.s

CMakeFiles/my_lib.dir/src/FullSystem/HessianBlocks.cpp.o: CMakeFiles/my_lib.dir/flags.make
CMakeFiles/my_lib.dir/src/FullSystem/HessianBlocks.cpp.o: ../src/FullSystem/HessianBlocks.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiang/Project/my_dso/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/my_lib.dir/src/FullSystem/HessianBlocks.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_lib.dir/src/FullSystem/HessianBlocks.cpp.o -c /home/xiang/Project/my_dso/src/FullSystem/HessianBlocks.cpp

CMakeFiles/my_lib.dir/src/FullSystem/HessianBlocks.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_lib.dir/src/FullSystem/HessianBlocks.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiang/Project/my_dso/src/FullSystem/HessianBlocks.cpp > CMakeFiles/my_lib.dir/src/FullSystem/HessianBlocks.cpp.i

CMakeFiles/my_lib.dir/src/FullSystem/HessianBlocks.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_lib.dir/src/FullSystem/HessianBlocks.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiang/Project/my_dso/src/FullSystem/HessianBlocks.cpp -o CMakeFiles/my_lib.dir/src/FullSystem/HessianBlocks.cpp.s

CMakeFiles/my_lib.dir/src/FullSystem/CoarseInitializer.cpp.o: CMakeFiles/my_lib.dir/flags.make
CMakeFiles/my_lib.dir/src/FullSystem/CoarseInitializer.cpp.o: ../src/FullSystem/CoarseInitializer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiang/Project/my_dso/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/my_lib.dir/src/FullSystem/CoarseInitializer.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_lib.dir/src/FullSystem/CoarseInitializer.cpp.o -c /home/xiang/Project/my_dso/src/FullSystem/CoarseInitializer.cpp

CMakeFiles/my_lib.dir/src/FullSystem/CoarseInitializer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_lib.dir/src/FullSystem/CoarseInitializer.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiang/Project/my_dso/src/FullSystem/CoarseInitializer.cpp > CMakeFiles/my_lib.dir/src/FullSystem/CoarseInitializer.cpp.i

CMakeFiles/my_lib.dir/src/FullSystem/CoarseInitializer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_lib.dir/src/FullSystem/CoarseInitializer.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiang/Project/my_dso/src/FullSystem/CoarseInitializer.cpp -o CMakeFiles/my_lib.dir/src/FullSystem/CoarseInitializer.cpp.s

CMakeFiles/my_lib.dir/src/FullSystem/PixelSelector2.cpp.o: CMakeFiles/my_lib.dir/flags.make
CMakeFiles/my_lib.dir/src/FullSystem/PixelSelector2.cpp.o: ../src/FullSystem/PixelSelector2.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiang/Project/my_dso/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/my_lib.dir/src/FullSystem/PixelSelector2.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_lib.dir/src/FullSystem/PixelSelector2.cpp.o -c /home/xiang/Project/my_dso/src/FullSystem/PixelSelector2.cpp

CMakeFiles/my_lib.dir/src/FullSystem/PixelSelector2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_lib.dir/src/FullSystem/PixelSelector2.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiang/Project/my_dso/src/FullSystem/PixelSelector2.cpp > CMakeFiles/my_lib.dir/src/FullSystem/PixelSelector2.cpp.i

CMakeFiles/my_lib.dir/src/FullSystem/PixelSelector2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_lib.dir/src/FullSystem/PixelSelector2.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiang/Project/my_dso/src/FullSystem/PixelSelector2.cpp -o CMakeFiles/my_lib.dir/src/FullSystem/PixelSelector2.cpp.s

CMakeFiles/my_lib.dir/src/FullSystem/FullSystem.cpp.o: CMakeFiles/my_lib.dir/flags.make
CMakeFiles/my_lib.dir/src/FullSystem/FullSystem.cpp.o: ../src/FullSystem/FullSystem.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiang/Project/my_dso/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/my_lib.dir/src/FullSystem/FullSystem.cpp.o"
	/usr/bin/x86_64-linux-gnu-g++-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_lib.dir/src/FullSystem/FullSystem.cpp.o -c /home/xiang/Project/my_dso/src/FullSystem/FullSystem.cpp

CMakeFiles/my_lib.dir/src/FullSystem/FullSystem.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_lib.dir/src/FullSystem/FullSystem.cpp.i"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiang/Project/my_dso/src/FullSystem/FullSystem.cpp > CMakeFiles/my_lib.dir/src/FullSystem/FullSystem.cpp.i

CMakeFiles/my_lib.dir/src/FullSystem/FullSystem.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_lib.dir/src/FullSystem/FullSystem.cpp.s"
	/usr/bin/x86_64-linux-gnu-g++-5 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiang/Project/my_dso/src/FullSystem/FullSystem.cpp -o CMakeFiles/my_lib.dir/src/FullSystem/FullSystem.cpp.s

# Object files for target my_lib
my_lib_OBJECTS = \
"CMakeFiles/my_lib.dir/src/util/DatasetReader.cpp.o" \
"CMakeFiles/my_lib.dir/src/util/Undistort.cpp.o" \
"CMakeFiles/my_lib.dir/src/util/MinimalImage.cpp.o" \
"CMakeFiles/my_lib.dir/src/util/globalutil.cpp.o" \
"CMakeFiles/my_lib.dir/src/util/ImageAndExposure.cpp.o" \
"CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp.o" \
"CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp.o" \
"CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp.o" \
"CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp.o" \
"CMakeFiles/my_lib.dir/src/FullSystem/HessianBlocks.cpp.o" \
"CMakeFiles/my_lib.dir/src/FullSystem/CoarseInitializer.cpp.o" \
"CMakeFiles/my_lib.dir/src/FullSystem/PixelSelector2.cpp.o" \
"CMakeFiles/my_lib.dir/src/FullSystem/FullSystem.cpp.o"

# External object files for target my_lib
my_lib_EXTERNAL_OBJECTS =

lib/libmy_lib.a: CMakeFiles/my_lib.dir/src/util/DatasetReader.cpp.o
lib/libmy_lib.a: CMakeFiles/my_lib.dir/src/util/Undistort.cpp.o
lib/libmy_lib.a: CMakeFiles/my_lib.dir/src/util/MinimalImage.cpp.o
lib/libmy_lib.a: CMakeFiles/my_lib.dir/src/util/globalutil.cpp.o
lib/libmy_lib.a: CMakeFiles/my_lib.dir/src/util/ImageAndExposure.cpp.o
lib/libmy_lib.a: CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageRW_OpenCV.cpp.o
lib/libmy_lib.a: CMakeFiles/my_lib.dir/src/IOWrapper/OpenCV/ImageDisplay_Opencv.cpp.o
lib/libmy_lib.a: CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/PangolinDSOViewer.cpp.o
lib/libmy_lib.a: CMakeFiles/my_lib.dir/src/IOWrapper/Pangolin/KeyFrameDisplay.cpp.o
lib/libmy_lib.a: CMakeFiles/my_lib.dir/src/FullSystem/HessianBlocks.cpp.o
lib/libmy_lib.a: CMakeFiles/my_lib.dir/src/FullSystem/CoarseInitializer.cpp.o
lib/libmy_lib.a: CMakeFiles/my_lib.dir/src/FullSystem/PixelSelector2.cpp.o
lib/libmy_lib.a: CMakeFiles/my_lib.dir/src/FullSystem/FullSystem.cpp.o
lib/libmy_lib.a: CMakeFiles/my_lib.dir/build.make
lib/libmy_lib.a: CMakeFiles/my_lib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/xiang/Project/my_dso/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Linking CXX static library lib/libmy_lib.a"
	$(CMAKE_COMMAND) -P CMakeFiles/my_lib.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/my_lib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/my_lib.dir/build: lib/libmy_lib.a

.PHONY : CMakeFiles/my_lib.dir/build

CMakeFiles/my_lib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/my_lib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/my_lib.dir/clean

CMakeFiles/my_lib.dir/depend:
	cd /home/xiang/Project/my_dso/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/xiang/Project/my_dso /home/xiang/Project/my_dso /home/xiang/Project/my_dso/build /home/xiang/Project/my_dso/build /home/xiang/Project/my_dso/build/CMakeFiles/my_lib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/my_lib.dir/depend
