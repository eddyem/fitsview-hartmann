cmake_minimum_required(VERSION 2.6)
if(DEFINED DEBUG)
	add_definitions(-DEBUG)
endif()
set(SRC ${CMAKE_SOURCE_DIR}/src)
aux_source_directory(${SRC} SOURCES)
set(NOCUFILE ${SRC}/NOCUDA.c)
set(CUFILE ${SRC}/CUDA.cu)
set(CFLAGS -O3 -Wall -Werror -W -std=c99)
set(LCPATH ${SRC}/locale/ru)
set(PO_FILE ${LCPATH}/messages.po)
set(MO_FILE ${LCPATH}/LC_MESSAGES/${PROJ}.mo)
set(RU_FILE ${LCPATH}/ru.po)
find_package(PkgConfig REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GTK2 REQUIRED)
pkg_check_modules(${PROJ} REQUIRED
		#glib-2.0>=2.10
		gtkglext-1.0>=0.7.0
		gtkglext-x11-1.0>=0.7.0
		#gtk+-2.0>=2.6.0
		cfitsio>=3.0
		fftw3>=3.2.0
		)
include(FindOpenMP)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lfftw3_threads")
if(NOT DEFINED NO_LEPTONICA)
	pkg_check_modules(LIBLEPT liblept)
endif()
if(NOT DEFINED NO_GSL)
	pkg_check_modules(GSL gsl)
endif()
if(NOT DEFINED GSL_VERSION)
	message("GSL not found, some mathematics functions wouldn't be avialable")
else()
	add_definitions(-DGSL_FOUND)
endif()
if(NOT DEFINED LIBLEPT_VERSION)
	message("Leptonica library not found, some functions wouldn't be avialable")
else()
	add_definitions(-DLEPTONICA_FOUND)
endif()
if(CUDA_FOUND)
	list(REMOVE_ITEM SOURCES ${NOCUFILE})
	list(APPEND CUDA_NVCC_FLAGS --use_fast_math)
	cuda_include_directories(include)
	cuda_add_executable(${PROJ} ${SOURCES} ${CUFILE} ${PO_FILE} ${MO_FILE} ui.h)
	target_link_libraries(	${PROJ} ${${PROJ}_LIBRARIES}
							${GSL_LIBRARIES} ${LIBLEPT_LIBRARIES}
							${CUDA_CUFFT_LIBRARIES} -lcuda)
else(CUDA_FOUND)
	find_package(Threads)
	add_executable(${PROJ} ${SOURCES} ${PO_FILE} ${MO_FILE} ui.h)
	target_link_libraries(	${PROJ} ${${PROJ}_LIBRARIES}
							${GSL_LIBRARIES} ${LIBLEPT_LIBRARIES}
							${CMAKE_THREAD_LIBS_INIT}
							)
endif(CUDA_FOUND)
include_directories(${SRC}/include ${${PROJ}_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${LIBLEPT_INCLUDE_DIRS})
link_directories(${${PROJ}_LIBRARY_DIRS} ${GSL_LIBRARY_DIRS} ${LIBLEPT_LIBRARY_DIRS})
add_definitions(-DPACKAGE_VERSION=\"0.0.1\"  -DGETTEXT_PACKAGE=\"${PROJ}\"
				-DLOCALEDIR=\"~/.local/share/locale\" ${CFLAGS})

#if(DEFINED DEBUG)
	find_package(Gettext REQUIRED)
	find_program(GETTEXT_XGETTEXT_EXECUTABLE xgettext)
	if(NOT GETTEXT_XGETTEXT_EXECUTABLE )
		  message(FATAL_ERROR "xgettext not found")
	endif(NOT GETTEXT_XGETTEXT_EXECUTABLE )
add_custom_command(
	OUTPUT ${PO_FILE}
	COMMAND ${GETTEXT_XGETTEXT_EXECUTABLE} --from-code=koi8-r ${SOURCES} ${SRC}/${PROJ}.glade -k_ -kN_ -o ${PO_FILE}
	DEPENDS ${SOURCES})
add_custom_command(
	OUTPUT ${MO_FILE}
	COMMAND ${GETTEXT_MSGFMT_EXECUTABLE} ${RU_FILE} -o ${MO_FILE}
	DEPENDS ${RU_FILE})
add_custom_command(
	OUTPUT ui.h
	COMMAND gtk-builder-convert ${SRC}/${PROJ}.glade ${PROJ}.ui
	COMMAND ${SRC}/scripts/genh ${PROJ}.ui
	DEPENDS ${PROJ}.glade)
#endif(DEFINED DEBUG)
