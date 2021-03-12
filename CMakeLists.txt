cmake_minimum_required(VERSION 3.0)

project(tubex-capd)



################################################################################
# Looking for Ibex, Tubex
################################################################################

  set(IBEX_ROOT "/home/neveub/tubex/ibex-lib/install") # Install path of ibex lib
  set(TUBEX_ROOT "/home/neveub/tubex/tubex-lib/")
  set(TUBEX_SOLVER_ROOT "/home/neveub/tubex/tubex-solve/")
 
  set(IBEX_LIB_DIR "${IBEX_ROOT}/lib"
                   "${IBEX_ROOT}/lib/ibex/3rd")
  
  set(IBEX_INCLUDE_DIR  "${IBEX_ROOT}/include"
		        "${IBEX_ROOT}/include/ibex"
		        "${IBEX_ROOT}/include/ibex/3rd"
			"${IBEX_ROOT}/include/lib" 
			"${IBEX_ROOT}/include/lib/ibex/3rd")

  set(TUBEX_INCLUDE_DIR "${TUBEX_ROOT}/build/include"
			"${TUBEX_ROOT}/build/src/core"
			"${TUBEX_ROOT}/include/tubex"
			"${TUBEX_ROOT}/include/lib")
  set(TUBEX_SOLVER_INCLUDE_DIR "${TUBEX_SOLVER_ROOT}/src" 
    "${TUBEX_SOLVER_ROOT}/build/src")

  set(TUBEX_LIB_DIR "${TUBEX_ROOT}/build/src/core")
    #"${TUBEX_ROOT}/build/src/robotics")
  set(TUBEX_SOLVER_LIB_DIR "${TUBEX_SOLVER_ROOT}/build/src")

    


  include_directories(${IBEX_INCLUDE_DIR})
  include_directories(${TUBEX_INCLUDE_DIR})
  include_directories(${TUBEX_SOLVER_INCLUDE_DIR})

  link_directories(${IBEX_LIB_DIR})
  link_directories(${TUBEX_LIB_DIR})
  link_directories(${TUBEX_SOLVER_LIB_DIR})
  
  

###########################################################"
set(CAPD_DIR "/home/neveub/capd/capd-capdDynSys-5.0.59")

set(CAPD_LIB_DIR "/home/neveub/capd/capd-capdDynSys-5.0.59/install/lib")

include_directories(${CAPD_DIR}/install/include/)

link_directories(${CAPD_LIB_DIR})

# CAPD LIBRARY BINARY
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math -D__USE_FILIB__")
set(CAPD_LIBS   tubex-solve tubex capd ibex prim)
include_directories(capd_tubex)
add_library(Capdtubex STATIC capd_tubex/tubex_CtcCapd.h capd_tubex/tubex_CtcCapd.cpp capd_tubex/tubex_TubeVectorODE.h capd_tubex/tubex_TubeVectorODE.cpp capd_tubex/tubex_capd2tubex.h capd_tubex/tubex_capd2tubex.cpp)
target_link_libraries(Capdtubex ${CAPD_LIBS})


add_executable(main examples/main.cpp)


target_link_libraries(main Capdtubex   ${CAPD_LIBS}  )
add_executable(ivp01 examples/ivp01.cpp)
target_link_libraries(ivp01 Capdtubex   ${CAPD_LIBS}  )
add_executable(ivp02 examples/ivp02.cpp)
target_link_libraries(ivp02 Capdtubex   ${CAPD_LIBS}  )
add_executable(ivp03 examples/ivp03.cpp)
target_link_libraries(ivp03 Capdtubex   ${CAPD_LIBS}  )

add_executable(ivp02-solver examples/ivp02-solver.cpp)
target_link_libraries(ivp02-solver Capdtubex   ${CAPD_LIBS}  )

add_executable(ivp12 examples/ivp12.cpp)
target_link_libraries(ivp12 Capdtubex   ${CAPD_LIBS}  )
add_executable(ivp13 examples/ivp13.cpp)
target_link_libraries(ivp13 Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp09 examples/bvp09.cpp)
target_link_libraries(bvp09 Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp04 examples/bvp04.cpp)
target_link_libraries(bvp04 Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp01a examples/bvpsolve/bvp01a.cpp)
target_link_libraries(bvp01a Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp01b examples/bvpsolve/bvp01b.cpp)
target_link_libraries(bvp01b Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp01c examples/bvpsolve/bvp01c.cpp)
target_link_libraries(bvp01c Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp01d examples/bvpsolve/bvp01d.cpp)
target_link_libraries(bvp01d Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp02a examples/bvpsolve/bvp02a.cpp)


target_link_libraries(bvp02a Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp02b examples/bvpsolve/bvp02b.cpp)
target_link_libraries(bvp02b Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp02c examples/bvpsolve/bvp02c.cpp)
target_link_libraries(bvp02c Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp02d examples/bvpsolve/bvp02d.cpp)
target_link_libraries(bvp02d Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp22a examples/bvpsolve/bvp22a.cpp)
target_link_libraries(bvp22a Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp22b examples/bvpsolve/bvp22b.cpp)
target_link_libraries(bvp22b Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp22c examples/bvpsolve/bvp22c.cpp)
target_link_libraries(bvp22c Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp25a examples/bvpsolve/bvp25a.cpp)
target_link_libraries(bvp25a Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp25b examples/bvpsolve/bvp25b.cpp)
target_link_libraries(bvp25b Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp25c examples/bvpsolve/bvp25c.cpp)
target_link_libraries(bvp25c Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp26a examples/bvpsolve/bvp26a.cpp)
target_link_libraries(bvp26a Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp26b examples/bvpsolve/bvp26b.cpp)
target_link_libraries(bvp26b Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp26c examples/bvpsolve/bvp26c.cpp)
target_link_libraries(bvp26c Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp27a examples/bvpsolve/bvp27a.cpp)
target_link_libraries(bvp27a Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp27b examples/bvpsolve/bvp27b.cpp)
target_link_libraries(bvp27b Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp32a examples/bvpsolve/bvp32a.cpp)
target_link_libraries(bvp32a Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp32b examples/bvpsolve/bvp32b.cpp)
target_link_libraries(bvp32b Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp34a examples/bvpsolve/bvp34a.cpp)
target_link_libraries(bvp34a Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp34b examples/bvpsolve/bvp34b.cpp)
target_link_libraries(bvp34b Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp34c examples/bvpsolve/bvp34c.cpp)
target_link_libraries(bvp34c Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp34d examples/bvpsolve/bvp34d.cpp)
target_link_libraries(bvp34d Capdtubex   ${CAPD_LIBS}  )

add_executable(bvp35a examples/bvpsolve/bvp35a.cpp)
target_link_libraries(bvp35a Capdtubex   ${CAPD_LIBS}  )

add_executable(bvp35b examples/bvpsolve/bvp35b.cpp)
target_link_libraries(bvp35b Capdtubex   ${CAPD_LIBS}  )
add_executable(bvp35c examples/bvpsolve/bvp35c.cpp)
target_link_libraries(bvp35c Capdtubex   ${CAPD_LIBS}  )

add_executable(bvp40 examples/bvpsolve/bvp40.cpp)
target_link_libraries(bvp40 Capdtubex   ${CAPD_LIBS}  )

###########################################################
