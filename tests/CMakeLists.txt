cmake_minimum_required(VERSION 3.1)
################################################################################


set(test_sources
    main.cpp

    ExactCCDTests.cpp
)

add_executable(${PROJECT_NAME}_tests ${test_sources})

foreach(source IN ITEMS ${test_sources})
   source_group("tests" FILES "${source}")
endforeach()

################################################################################
# Required Libraries
################################################################################

target_link_libraries(${PROJECT_NAME}_tests PUBLIC ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}_tests PUBLIC catch)

################################################################################
# Add Tests
################################################################################



# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(${PROJECT_NAME}_tests)
