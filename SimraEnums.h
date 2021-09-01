//==============================================================================
//!
//! \file SimraEnums.h
//!
//! \date Jul 6 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various enumerations in used SIMRA applications.
//!
//==============================================================================

#ifndef SIMRA_ENUMS_H_
#define SIMRA_ENUMS_H_


//! \brief Enumeration of solution components.
enum SolutionEntry {
  U_X = 0, //!< X velocity
  U_Y, //!< Y velocity
  U_Z, //!< Z velocity
  PS, //!< Hydrostatic pressure
  TK, //!< Turbulent kinetic energy
  TD, //!< Turbulent dissipation
  VTEF, //!< Effective viscosity
  PT, //!< Potential temperature
  PTS, //!< Hydrostatic temperature
  RHO, //!< Densitity
  RHOS, //!< Stratified densities
  STRAT //!< Stratification
};


//! \brief Enumeration of result file types.
enum ResultsType {
  BOUNDARY_FILE, //!< Boundary condition file
  RESTART_FILE, //!< Restart file, holds one time step
  HISTORY_FILE, //!< History file, holds multiple time steps
  INIT_FILE     //!< Init file, holds one time step
};

#endif
