include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/src)

add_executable(demo1 demo1.c)
target_link_libraries(demo1 tabr ${MATH_LIBRARY})

if (THREADS_FOUND)
  add_executable(demothr1 demothr1.c)
  target_link_libraries(demothr1 tabr ${MATH_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
endif (THREADS_FOUND)
