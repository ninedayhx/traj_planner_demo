cmake_minimum_required(VERSION 3.0.2)
project(traj_plan)

message(STATUS "a planner demo")
message(STATUS "authur:ninedayhx")
message(STATUS "email:1170535490@qq.com")

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)      
ENDIF()
MESSAGE("Build type: " ${CMAKE_BUILD_TYPE}) 

#  -fopenmp is important 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 ")

set(CATKIN_DEPENDS
	roscpp
	rospy
	nav_msgs
	sensor_msgs
	std_msgs
	tf
	# for grid map
	grid_map_core
	grid_map_msgs
	grid_map_ros
	pcl_ros
	pcl_conversions
)
set(GRID_MAP_FILES
	src/GridMapPclConverter.cpp
	src/GridMapPclLoader.cpp
	src/helpers.cpp
	src/PclLoaderParameters.cpp
	src/PointcloudProcessor.cpp
)


find_package(Eigen3 REQUIRED)

find_package(PCL REQUIRED)

find_package(OpenMP QUIET)
if (OpenMP_FOUND)
  add_compile_options("${OpenMP_CXX_FLAGS}")
  add_definitions(-DGRID_MAP_PCL_OPENMP_FOUND=${OpenMP_FOUND})
endif()

find_package(catkin REQUIRED
COMPONENTS
	${CATKIN_DEPENDS}
)

catkin_package(
	#头文件路径
	INCLUDE_DIRS
		include
	# 本项目中构建的库，构建可以在后面
	LIBRARIES
		${grid_map}
		yaml-cpp
	# catkin依赖项
	CATKIN_DEPENDS
		${CATKIN_DEPENDS}
	# 系统依赖项
	DEPENDS
		PCL
)


include_directories(
	include 
		${catkin_INCLUDE_DIRS}
		${EIGEN3_INCLUDE_DIR}
)

###########
## Build ##
###########
add_library(grid_map SHARED
	${GRID_MAP_FILES}
)
target_include_directories(grid_map PRIVATE
	include
)
target_include_directories(grid_map SYSTEM PUBLIC
	${catkin_INCLUDE_DIRS}
	${EIGEN3_INCLUDE_DIR}
	${PCL_INCLUDE_DIR}
)
target_link_libraries(grid_map
	${catkin_LIBRARIES}
	${PCL_LIBRARIES}
	yaml-cpp
)

add_library(minimum_traj
	src/minimum_traj.cpp
)
target_link_libraries(minimum_traj 
	${EIGEN3_LIBS}
)

add_executable(${PROJECT_NAME}_node 
	src/traj_plan_node.cpp
)
target_link_libraries(${PROJECT_NAME}_node
	minimum_traj
	${catkin_LIBRARIES}
	${EIGEN3_LIBRARIES}
	${PCL_LIBRARIES}
)

add_executable(map_server_node 
	src/map_server_node.cpp
)
target_link_libraries(map_server_node 
	grid_map
	${catkin_LIBRARIES}
	${PCL_LIBRARIES}
)
