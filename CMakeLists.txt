################################################################################
# General Information
################################################################################

cmake_minimum_required(VERSION 3.3)
project(RationalCCD)

################################################################################
set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

# Detects whether this is a top-level project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(ECCD_TOPLEVEL_PROJECT ON)
else()
    set(ECCD_TOPLEVEL_PROJECT OFF)
endif()

################################################################################
# Build static binaries
#set(BUILD_SHARED_LIBS OFF)
set(BUILD_STATIC_RUNTIME OFF)
set(BUILD_DYNAMIC_RUNTIME ON)

# Use folder in Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Export compile flags (used for autocompletion of the C++ code)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Generate position independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#ADD_DEFINITIONS("-D_DLL")

include(ECCDDownloadExternal)

# if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
# 	# If we use static library, we link with the static C++ runtime.
# 	foreach(config ${CMAKE_CONFIGURATION_TYPES})
# 		string(TOUPPER ${config} config)
# 		string(REPLACE /MD /MT CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
# 		string(REPLACE /MD /MT CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
# 		string(REPLACE /MDd /MTd CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
# 		string(REPLACE /MDd /MTd CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
# 	endforeach()
# endif()

################################################################################
# 3rd party libraries
################################################################################
find_package(GMP)

IF(NOT ${GMP_FOUND})
        MESSAGE(FATAL_ERROR "Cannot find GMP")
ENDIF()

# Eigen
if(NOT TARGET Eigen3::Eigen)
  eccd_download_eigen()
  add_library(${PROJECT_NAME}_eigen INTERFACE)
  target_include_directories(${PROJECT_NAME}_eigen SYSTEM INTERFACE
    $<BUILD_INTERFACE:${THIRD_PARTY_DIR}/eigen>
    $<INSTALL_INTERFACE:include>
  )
  set_property(TARGET ${PROJECT_NAME}_eigen PROPERTY EXPORT_NAME Eigen3::Eigen)
  add_library(Eigen3::Eigen ALIAS ${PROJECT_NAME}_eigen)
endif()

################################################################################

set(SOURCE_FILES
	src/ECCD.cpp
	src/ECCD.hpp
	src/Rational.hpp

	src/Utils.cpp
	src/Utils.hpp


	src/Plots.hpp
	src/Plots.cpp
	)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_library(${PROJECT_NAME}  ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

target_link_libraries(${PROJECT_NAME} PUBLIC gmp::gmp)







# Compile extras only if this is a top-level project
if(ECCD_TOPLEVEL_PROJECT)
	add_executable(${PROJECT_NAME}_bin src/main.cpp)
	target_link_libraries(${PROJECT_NAME}_bin ${PROJECT_NAME})


	eccd_download_catch2()
	list(APPEND CMAKE_MODULE_PATH ${THIRD_PARTY_DIR}/Catch2/contrib)

	add_library(catch INTERFACE)
	target_include_directories(catch SYSTEM INTERFACE ${THIRD_PARTY_DIR}/Catch2/single_include/catch2)

	include(CTest)
	enable_testing()
	add_subdirectory(tests)
endif()
