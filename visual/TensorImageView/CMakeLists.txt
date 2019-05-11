
cmake_minimum_required(VERSION 2.8)

project(TensorImageView)
 


find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
 

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
add_executable(TensorImageView TensorImageView.cxx)
target_link_libraries(TensorImageView
  ${VTK_LIBRARIES} ${ITK_LIBRARIES})