include(FetchContent)

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        release-1.11.0
)

FetchContent_MakeAvailable(googletest)
add_library(GTest::GTest INTERFACE IMPORTED)
target_link_libraries(GTest::GTest INTERFACE gtest_main gmock_main)

add_executable(symmetry_test symmetry_test.cpp)

target_link_libraries(symmetry_test
        PRIVATE
        GTest::GTest
        symmetry
        )

add_test(symmetry_gtests symmetry_test)
