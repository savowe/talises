PROJECT (atus2)

cmake_minimum_required(VERSION 3.1)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

message("**********************************************************************")
execute_process(COMMAND sh -c "module list")
message("**********************************************************************")

find_package(MPI)
find_package(Boost REQUIRED)
find_package(GSL REQUIRED)
find_package(FFTW REQUIRED)
find_package(MUPARSER REQUIRED)
message("**********************************************************************")


#SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS_RELEASE "-std=gnu++14 -march=native -O3 -funroll-loops -ftree-vectorize -fopenmp -w -s -Wall")
#set(CMAKE_CXX_FLAGS_DEBUG "-std=gnu++14 -g -Wall -Wextra -fopenmp -fsanitize=thread")
#set(CMAKE_CXX_FLAGS_DEBUG "-std=gnu++14 -g -Wall -Wextra -fopenmp -fsanitize=address")
set(CMAKE_CXX_FLAGS_DEBUG "-std=gnu++14 -g -Wall -Wextra -fopenmp -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug" )
endif()


set( HOME $ENV{HOME} CACHE STRING INTERNAL )
set( DIR_INC ${PROJECT_SOURCE_DIR}/include CACHE STRING INTERNAL )
set( DIR_MYLIB ${PROJECT_SOURCE_DIR}/source/libs/lib_myutils CACHE STRING INTERNAL )
set( DIR_MYLIB_MPI ${PROJECT_SOURCE_DIR}/source/libs/lib_myutils_mpi_2 CACHE STRING INTERNAL )


set( EXECUTABLE_OUTPUT_PATH  ${HOME}/bin )

include_directories( ${MPI_CXX_INCLUDE_PATH} 
                     ${Boost_INCLUDE_DIRS} 
                     ${DIR_INC} ${DIR_MYLIB} 
                     ${DIR_MYLIB_MPI}/src_mpi 
                     ${GSL_INCLUDE_DIR} 
                     ${FFTW_INCLUDE_DIR} 
                     ${MUPARSER_INCLUDE_DIR} 
)

# enable profiling
#set( CMAKE_EXE_LINKER_FLAGS -pg )

add_subdirectory( src )



# FILE(GLOB bash_sh "${PROJECT_SOURCE_DIR}/Bash/*")
#
# FOREACH( file_i ${bash_sh})
#     MESSAGE(STATUS ${file_i} )
#     INSTALL(FILES ${file_i} PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION "${HOME}/bin" )
# ENDFOREACH( file_i )