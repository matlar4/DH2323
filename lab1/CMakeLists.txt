# Top-Level CMakeList.txt

cmake_minimum_required (VERSION 2.6)
project ( FirstLab )

add_executable( FirstLab skeleton.cpp)
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
set(SDL_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/SDL-1.2.15/include)
#find_package (SDL REQUIRED)

add_custom_target(ExternalDlls
COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/SDL_1.2_x86_vs2015/SDL.dll          ${CMAKE_BINARY_DIR}/Release
COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/SDL_1.2_x86_vs2015/SDL.dll          ${CMAKE_BINARY_DIR}/Debug
)

include_directories(
	${SDL_INCLUDE_DIR}
	${PROJECT_SOURCE_DIR}/glm
)
#link_libraries(${SDL_LIBRARY})
#link_directories(${SDL_LIB_DIR})
target_link_libraries(FirstLab ${PROJECT_SOURCE_DIR}/SDL_1.2_x86_vs2015/SDL.lib ${PROJECT_SOURCE_DIR}/SDL_1.2_x86_vs2015/SDLmain.lib)
add_dependencies(FirstLab ExternalDlls)