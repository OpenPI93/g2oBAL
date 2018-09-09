cmake_minimum_required( VERSION 2.8 )
project( BAProblem )

set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "-std=c++11 -O3" )

# 添加cmake模块以使用g2o
list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules )

find_package( OpenCV 3.1 REQUIRED )
find_package(Sophus REQUIRED)
find_package( G2O REQUIRED )
find_package(Eigen3 REQUIRED)
find_package( CSparse REQUIRED )
Find_Package(Cholmod REQUIRED)
find_package(Pangolin REQUIRED)

include_directories(
    ${OpenCV_INCLUDE_DIRS}
    ${G2O_INCLUDE_DIRS}
    ${CSPARSE_INCLUDE_DIR}
    ${Sophus_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${Pangolin_INCLUDE_DIRS}

)

add_executable( BAProblem main.cpp balproblem.cpp)
target_link_libraries( BAProblem
   ${OpenCV_LIBS}
   ${CSPARSE_LIBRARY}
   ${Sophus_LIBS}
   ${EIGEN3_LIBS}
   ${CHOLMOD_LIBRARIES}
   ${Pangolin_LIBRARIES}
   g2o_core g2o_stuff g2o_types_sba g2o_csparse_extension
)
