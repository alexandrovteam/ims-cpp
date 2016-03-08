cmake_minimum_required(VERSION 2.8)
project(SM_standalone)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
    message(FATAL_ERROR "GCC version must be at least 4.9!")
  endif()
else()
  message(FATAL_ERROR "Please use recent GCC to compile this project")
endif()

include(ExternalProject)

if ( APPLE )
  # necessary for GCC to use system assembler
  # http://stackoverflow.com/questions/9840207/how-to-use-avx-pclmulqdq-on-mac-os-x-lion
  set(CMAKE_C_FLAGS "-Wa,-q -fvisibility-inlines-hidden")
  set(CMAKE_CXX_FLAGS "-Wa,-q -fvisibility-inlines-hidden")
elseif( WIN32 )
  # necessary for FFTW to work on 32-bit, as mentioned e.g. here: https://github.com/libgit2/libgit2/issues/3342
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mincoming-stack-boundary=2")
else()
  set(CMAKE_C_FLAGS "-fPIC")
  set(CMAKE_CXX_FLAGS "-fPIC")
endif()

# MessagePack - simple serialization library
externalproject_add(msgpack
  GIT_REPOSITORY https://github.com/msgpack/msgpack-c
  GIT_TAG cpp-1.4.0
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  ) # for C++ compilation is not required
ExternalProject_Get_Property(msgpack source_dir)
include_directories(${source_dir}/include)

find_package(ZLIB REQUIRED)
# blosc - fast compressor for numeric data
externalproject_add(blosc
  GIT_REPOSITORY https://github.com/Blosc/c-blosc
  GIT_TAG v1.7.1
  INSTALL_COMMAND ""
  BUILD_COMMAND make blosc_static
  CMAKE_ARGS -DBUILD_TESTS:BOOL=OFF -DBUILD_BENCHMARKS:BOOL=OFF -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
             -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
             -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
             -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  )
ExternalProject_Get_Property(blosc source_dir)
include_directories(${source_dir}/blosc/)
ExternalProject_Get_Property(blosc binary_dir)
set(BLOSC_LIBRARY ${binary_dir}/blosc/libblosc.a ${ZLIB_LIBRARY})

# cxxopts - lightweight command line options parser
externalproject_add(cxxopts
  GIT_REPOSITORY https://github.com/jarro2783/cxxopts
  GIT_TAG aec97a6f53c3486fc51e0d9857f10b683180d668
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  ) # header-only library, requires std::regex
ExternalProject_Get_Property(cxxopts source_dir)
include_directories(${source_dir}/src)

set(PROGRESSBAR_DIR ${PROJECT_SOURCE_DIR}/third-party/progressbar)
externalproject_add(progressbar
  GIT_REPOSITORY https://github.com/doches/progressbar
  GIT_TAG 760e20c245a72cf564687c1ecc18163799229c58
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
  DOWNLOAD_DIR ${PROGRESSBAR_DIR}
  SOURCE_DIR ${PROGRESSBAR_DIR}
  BINARY_DIR ${PROGRESSBAR_DIR}
  BUILD_COMMAND ${CMAKE_C_COMPILER} ${CMAKE_C_FLAGS} -c -Iinclude lib/progressbar.c  -o lib/libprogressbar.a
  )
ExternalProject_Get_Property(progressbar source_dir)
include_directories(${source_dir}/include)
set(PROGRESSBAR_LIB ${PROGRESSBAR_DIR}/lib/libprogressbar.a)

find_package(OpenMP REQUIRED)
find_package(LibXml2 REQUIRED)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)
include_directories(${LIBXML2_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11 -pthread")


find_library(FFTW_LIBRARY NAMES fftw3 fftw fftw-3.3)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are strongly discouraged. "
    "Please remove CMakeCache.txt and the CMakeFiles/ directory, then build out-of-source.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")
if (NOT APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++ -static-libgcc")
# on Mac OS X static linking leads to exceptions not being caught,
# so it's better to just use the 'delocate-wheel' tool (puts needed shared libs into the wheel)
endif()

# necessary for static compilation
add_definitions(-DIN_LIBXML)

set(SRC_COMMON
  imzb/reader.cpp
  imzb/writer.cpp
  ims/image_measures.cpp
  ms/isocalc.cpp
  ms/isotope_pattern.cpp
  ms/periodic_table.cpp
  )

add_library(libims STATIC ${SRC_COMMON})
add_dependencies(libims blosc)
set_target_properties(libims PROPERTIES OUTPUT_NAME ims)

if ( NOT APPLE )
target_link_libraries(libims ${BLOSC_LIBRARY})
endif()

add_executable(ims-bin
  tools/main.cpp
  tools/convert.cpp
  tools/isocalc.cpp
  tools/detect.cpp
  tools/fdr.cpp
  imzml/reader.cpp
  utils/isotope_pattern_db.cpp
  utils/metrics.cpp)
set_target_properties(ims-bin PROPERTIES OUTPUT_NAME ims)
if ( WIN32 )
  target_link_libraries(ims-bin ${PROGRESSBAR_LIB} -Wl,-Bstatic -lncurses -lxml2 -liconv libims -lfftw3 -lgomp -lwinpthread -lws2_32)
elseif ( APPLE )
  target_link_libraries(ims-bin ${BLOSC_LIBRARY} ${PROGRESSBAR_LIB} -lncurses -lxml2 libims /usr/local/opt/fftw/lib/libfftw3.a)
else()
  target_link_libraries(ims-bin ${PROGRESSBAR_LIB} -Wl,-Bstatic -ltermcap -lxml2 libims -lfftw3 -llzma -lstdc++ -Wl,-Bdynamic -lm)
endif()
add_dependencies(ims-bin blosc libims msgpack progressbar)

if ( SCILS_H5_SUPPORT)
target_compile_definitions(ims-bin PRIVATE -DSCILS_H5)
target_link_libraries(ims-bin -lhdf5 -lhdf5_hl)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/)

add_library(ms_cffi SHARED cffi/isocalc.cpp cffi/common.cpp ms/isocalc.cpp ms/isotope_pattern.cpp ms/periodic_table.cpp)
if ( WIN32 )
  target_link_libraries(ms_cffi -Wl,-Bstatic -lfftw3 -lwinpthread -lws2_32)
elseif ( APPLE )
  target_link_libraries(ms_cffi /usr/local/opt/fftw/lib/libfftw3.a)
else()
  target_link_libraries(ms_cffi -Wl,-Bstatic -lfftw3 -lstdc++ -Wl,-Bdynamic -lm)
endif()

add_library(imzb_cffi SHARED cffi/imzb.cpp cffi/common.cpp imzb/reader.cpp)
target_link_libraries(imzb_cffi ${BLOSC_LIBRARY} -Wl,-Bstatic -lstdc++ -Wl,-Bdynamic -lm)

# export compilation database for autocompletion
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
