
set(LIBTABR_SOURCE tabr.c)
set(LIBTABR_INCLUDE tabr.h)

if (THREADS_FOUND)
  set(LIBTABR_SOURCE ${LIBTABR_SOURCE} threadman.c)
  set(LIBTABR_INCLUDE ${LIBTABR_INCLUDE} tabr_threadman.h)
endif (THREADS_FOUND)

add_library(tabr ${LIBTABR_SOURCE})
install(TARGETS tabr DESTINATION lib)
install(FILES ${LIBTABR_INCLUDE} DESTINATION include)
