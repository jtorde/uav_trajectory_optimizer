cmake_minimum_required(VERSION 2.8.3)
project(cvx)

## Add support for C++11
add_definitions(-std=c++11)

find_package( Eigen3 REQUIRED )
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(catkin REQUIRED)

include_directories(${catkin_INCLUDE_DIRS})

add_executable(${PROJECT_NAME}_exec src/main.cpp)
target_link_libraries(${PROJECT_NAME}_exec ${catkin_LIBRARIES} ${CMAKE_CURRENT_LIST_DIR}/src/solvers/cvxgen/libcvx_vel.so ${CMAKE_CURRENT_LIST_DIR}/src/solvers/cvxgen/libcvx_accel.so ${CMAKE_CURRENT_LIST_DIR}/src/solvers/cvxgen/libcvx_jerk.so)
add_dependencies(${PROJECT_NAME}_exec ${catkin_EXPORTED_TARGETS})


