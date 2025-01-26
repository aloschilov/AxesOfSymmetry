add_library(symmetry symmetry.cpp matrix.cpp)
add_executable(symmetry_main main.cpp symmetry.cpp matrix.cpp)


target_include_directories(symmetry PUBLIC include)
target_include_directories(symmetry_main PUBLIC include)
