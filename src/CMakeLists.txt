include_directories (${PROJECT_SOURCE_DIR}/include)

add_executable (mm main.cpp)

set_property (TARGET mm PROPERTY CXX_STANDARD 17)

find_package (OpenMP REQUIRED)

if (OpenMP_CXX_FOUND)
    target_link_libraries (mm PRIVATE OpenMP::OpenMP_CXX)
endif ()
