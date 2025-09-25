cmake_minimum_required(VERSION 3.22)

project(SIMRA-PostProc LANGUAGES C CXX)

if(NOT TARGET IFEM)
  if(IFEM_AS_SUBMODULE)
    add_subdirectory(../.. IFEM)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR}/../..)
  else()
    get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
    set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../../${BUILD_DIR})
    find_package(IFEM REQUIRED)
  endif()
endif()

ifem_add_library(
  NAME
    SIMRA
  SOURCES
    ASMs3DSimra.C
    FSWallDistance.C
    SimraIntegrand.C
    SimraIO.C
    SIMFSWallDistance.C
    SIMSimraBase.C
    SIMSimraProject.C
    SIMSimraTransfer.C
  HEADERS
    ASMs3DSimra.h
    FSWallDistance.h
    SIMFSWallDistance.h
    SimraIntegrand.h
    SimraIO.h
    SIMSimraBase.h
    SIMSimraProject.h
    SIMSimraTransfer.h
  LIBRARIES
    IFEM
)

ifem_add_application(
  NAME
    SIMRA-PostProc
  SOURCES
    main_PostProc.C
  LIBRARIES
    SIMRA
)

ifem_add_application(
  NAME
    SIMRA-SolutionTransfer
  SOURCES
    main_SolutionTransfer.C
  LIBRARIES
    SIMRA
)

ifem_add_application(
  NAME
    SIMRA-CompareSetup
  SOURCES
    main_CompareSetup.C
  LIBRARIES
    SIMRA
)

# Regression tests
enable_testing()

ifem_add_hdf5_test(
  TARGET
    SIMRA-PostProc
  TEST_FILES
    Box.hreg
)

ifem_add_vtf_test(
  TARGET
    SIMRA-PostProc
  TEST_FILES
    Box.vreg
)

ifem_add_regression_test(
  TARGET
    SIMRA-PostProc
  TEST_FILES
    Box.reg
)

if(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  add_test(
    NAME
      SIMRA-SolutionTransfer+Box_transfer
    COMMAND
      ${PROJECT_SOURCE_DIR}/scripts/runCompareTest.sh
      $<TARGET_FILE:SIMRA-SolutionTransfer>
      $<TARGET_FILE:SIMRA-CompareSetup>
      Box-transfer.xinp
      ${PROJECT_SOURCE_DIR}/Test
      1e-5
  )
  add_test(
    NAME
      SIMRA-SolutionTransfer+Box_classify
    COMMAND
      ${PROJECT_SOURCE_DIR}/scripts/runCompareTest.sh
      $<TARGET_FILE:SIMRA-SolutionTransfer>
      $<TARGET_FILE:SIMRA-CompareSetup>
      Box-classify.xinp
      ${PROJECT_SOURCE_DIR}/Test
      1e-5
    )
endif()

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  ifem_add_check_target()
endif()

# For generating the doxy
ifem_add_doc_target(TARGET SIMRA-PostProc)

# Installation
include(GNUInstallDirs)
install(
  TARGETS
    SIMRA-PostProc SIMRA-SolutionTransfer
  DESTINATION
    ${CMAKE_INSTALL_BINDIR}
)
