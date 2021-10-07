//==============================================================================
//!
//! \file SimraFieldGenerator.h
//!
//! \date Sep 1 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Generators for SIMRA solution fields.
//!
//==============================================================================

#ifndef SIMRA_FIELD_GENERATOR_H_
#define SIMRA_FIELD_GENERATOR_H_

#include <array>

#include "MatVec.h"


class ASMs3DSimra;
class TiXmlElement;

/*!
 * \brief Class for generating idealized initial conditions for SIMRA models.
*/

class SimraFieldGenerator {
public:
  //! \brief Constructor.
  //! \param solution Reference to simulators solution vectors
  SimraFieldGenerator(Vectors& solution);

  //! \brief Set patch for mesh.
  //! param patch Pointer to patch to use
  void setPatch(const ASMs3DSimra* patch);

  //! \brief Read generator settings from an XML element.
  //! \param elem Element to read
  void read(const TiXmlElement* elem);

  //! \brief Print generator settings to terminal.
  void print();

  //! \brief Generate hydrostatic conditions.
  void hydrostaticConditions();

  //! \brief Returns true if we should generate a uniform velocity field.
  bool generateVelocity() const { return uni_velocity; }

  //! \brief Generate a uniform velocity field.
  void uniformVelocity();

  //! \brief Generate compatible turbulence fields.
  void turbulenceFields();

  //! \brief Returns configured constant surface roughness.
  double surfaceRoughness() const { return z0; }

private:
  //! \brief Generates a 1D temperature profile.
  std::vector<double> potentialTemperature() const;

  //! \brief Generates 1D hydrostatic pressure and density profiles.
  std::array<std::vector<double>,3>
  hydrostaticPressureDensity (const std::vector<double>& T_prof) const;

  const ASMs3DSimra* pch = nullptr; //!< Pointer to mesh
  std::array<size_t,3> n; //!< Cartesian dimensions of mesh
  Vectors& solution; //!< Reference to simulators solution vectors

  bool uni_velocity = false;  //!< True to generate uniform velocity
  double bl_height = 2000.0; //!< Height of boundary layer
  double wind_speed = 20.0; //!< Wind speed at boundary layer
  double wind_dir = 330.0; //!< Wind direction in degrees
  double z0 = 0.3; //!< Surface roughness
};

#endif
