project(SIMRA-PostProc)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules)

if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")
endif()

set(SIMRA_SRCS ASMs3DSimra.C
               FSWallDistance.C
               SimraIntegrand.C
               SimraIO.C
               SIMFSWallDistance.C
               SIMSimraBase.C
               SIMSimraProject.C
               SIMSimraTransfer.C
               ASMs3DSimra.h
               FSWallDistance.h
               SIMFSWallDistance.h
               SimraIntegrand.h
               SimraIO.h
               SIMSimraBase.h
               SIMSimraProject.h
               SIMSimraTransfer.h)

include_directories(${IFEM_INCLUDES} ${PROJECT_SOURCE_DIR} ../Common)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

add_library(SIMRA ${SIMRA_SRCS})

add_executable(SIMRA-PostProc main_PostProc.C)
add_executable(SIMRA-SolutionTransfer main_SolutionTransfer.C)
add_executable(SIMRA-CompareSetup main_CompareSetup.C)
target_link_libraries(SIMRA-PostProc SIMRA IFEMAppCommon ${IFEM_LIBRARIES})
target_link_libraries(SIMRA-SolutionTransfer SIMRA IFEMAppCommon ${IFEM_LIBRARIES})
target_link_libraries(SIMRA-CompareSetup SIMRA IFEMAppCommon ${IFEM_LIBRARIES})

# Installation
install(TARGETS SIMRA-PostProc SIMRA-SolutionTransfer DESTINATION bin COMPONENT bin)

# For generating the doxy
set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}")
add_doc_target(SIMRA-PostProc SIMRA-PostProc)

# Regression tests
enable_testing()
include(IFEMTesting)
if(MPI_FOUND)
endif()
if(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  ifem_add_test(Box.reg SIMRA-PostProc)
  ifem_add_vtf_test(Box.vreg SIMRA-PostProc)
  ifem_add_hdf5_test(Box.hreg SIMRA-PostProc)
  add_test(Box_transfer ${PROJECT_SOURCE_DIR}/scripts/runCompareTest.sh
                        ${PROJECT_BINARY_DIR}/bin/SIMRA-SolutionTransfer
                        ${PROJECT_BINARY_DIR}/bin/SIMRA-CompareSetup
                        Box-transfer.xinp
                        ${PROJECT_SOURCE_DIR}/Test
                        1e-5)
  add_test(Box_classify ${PROJECT_SOURCE_DIR}/scripts/runCompareTest.sh
                        ${PROJECT_BINARY_DIR}/bin/SIMRA-SolutionTransfer
                        ${PROJECT_BINARY_DIR}/bin/SIMRA-CompareSetup
                        Box-classify.xinp
                        ${PROJECT_SOURCE_DIR}/Test
                        1e-5)
endif()
list(APPEND TEST_APPS SIMRA-PostProc SIMRA-SolutionTransfer SIMRA-CompareSetup)

# Unit tests
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  SIMRA-PostProc 
                  IFEMAppCommon ${IFEM_LIBRARIES})

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
