# Top-Level CMakeList.txt

cmake_minimum_required (VERSION 2.6)
project ( ThirdLab )

add_executable( ThirdLab skeleton.cpp)

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
set(SDL_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/SDL-1.2.15/include)

add_custom_target(ExternalDlls
COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/SDL_1.2_x86_vs2015/SDL.dll          ${CMAKE_BINARY_DIR}/Release
COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/SDL_1.2_x86_vs2015/SDL.dll          ${CMAKE_BINARY_DIR}/Debug
)


include_directories(
	${SDL_INCLUDE_DIR}
	${PROJECT_SOURCE_DIR}/glm
)
#link_libraries(${SDL_LIBRARY})
target_link_libraries(ThirdLab ${PROJECT_SOURCE_DIR}/SDL_1.2_x86_vs2015/SDL.lib ${PROJECT_SOURCE_DIR}/SDL_1.2_x86_vs2015/SDLmain.lib)
add_dependencies(ThirdLab ExternalDlls)


