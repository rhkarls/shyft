# CMake file for compiling the C++ api library
set(api_sources api.cpp api_serialization.cpp )
add_library(shyft_api STATIC ${api_sources})
