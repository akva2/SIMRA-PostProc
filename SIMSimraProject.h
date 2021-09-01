//==============================================================================
//!
//! \file SIMSimraProject.h
//!
//! \date Jul 6 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for projection of SIMRA results.
//!
//==============================================================================

#ifndef SIMSIMRAPROJECT_H
#define SIMSIMRAPROJECT_H

#include "SimraEnums.h"
#include "SimraIntegrand.h"
#include "SimraIO.h"
#include "SIMSimraBase.h"

#include <fstream>


class DataExporter;

/*!
  \brief Simulation driver for projection of SIMRA results.
*/

class SIMSimraProject : public SIMSimraBase
{
public:
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

  //! \brief Default constructor.
  //! \param context The base xml tag to parse
  SIMSimraProject(const std::string& context = "simra");

  //! \brief Empty destructor.
  virtual ~SIMSimraProject() { myInts.clear(); myProblem = nullptr; }

  //! \brief Read results from the result file.
  bool readResults();

  //! \brief Returns solution time.
  double getSolutionTime() const { return solTime; }

  //! \brief Returns primary solution vector.
  Vector getSolution() const;

  //! \brief Set the solution vector.
  void setSolutions(const Vector& sol);

  //! \brief Returns a reference to wall distance.
  Vector& getDistance();

  //! \brief Returns a reference to the element-wise pressures.
  const Vector& getElmPressures() const { return itg.elmPressure; }

  //! \brief Returns a mutable reference to the element-wise pressures.
  Vector& getElmPressures() { return itg.elmPressure; }

  //! \brief Write all solution vectors to VTF.
  //! \param nBlock Running VTF block counter
  //! \param sol Primary solution vectors
  bool writeSolutionVectors(int& nBlock, const Vector& sol);

  //! \brief Print solution norms to terminal.
  //! \param gNorm Global norms
  void printSolutionNorms(const Vectors& gNorm) const;

  //! \brief Apply post-processing tasks on norms.
  bool postProcessNorms(Vectors& gNorm, Matrix* eNormp) override;

  //! \brief Register fields for exporter output.
  void registerFields(DataExporter& exporter, const Vector& sol,
                      const Vectors& projs, const Matrix& eNorm) const;

  //! \brief Name for data exporter.
  std::string getName() const override { return "Simra"; }

  //! \brief Returns whether to calculate terrain distance based on orthogonal mesh.
  //! \details If enabled, this function also calculates the distance
  bool orthogonalDistance();

  //! \brief Returns a single solution vector.
  //! \param idx 0-based index for solution to obtain
  const Vector& getSol(size_t idx) const;

  //! \brief Returns the file name of the result file.
  const std::string& getResultFile() const { return resultFile; }

  //! \brief Returns current step.
  int currentStep() const { return iStep; }

protected:
  //! \brief Prints a norm group to the log stream.
  void printNormGroup(const Vector& rNorm,
                      const Vector& fNorm,
                      const std::string& name) const override;

  //! \brief Print exact solution and error norms.
  void printExactNorms(const Vector& gNorm, size_t w) const;

  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  bool parse(const TiXmlElement *elem) override;

  ResultsType rType = RESTART_FILE; //!< Type for result file
  std::string resultFile; //!< File with result vectors
  Vectors solution; //!< Solution vector
  bool stratified = true; //!< Include temperature
  bool calcYp = true; //!< True to calculate y+ term
  bool useOrthogonalMesh = false; //!< True to assume terrain-orthogonal mesh lines.
  double solTime; //!< Solution time
  double uRef = 1.0; //!< Reference velocity value
  double lRef = 1.0; //!< Reference length scale

  SimraIntegrand itg; //!< Integrand to use

  int iStep = 0; //!< Current time step to read
  int initStep = 0; //!< Initial step to read from history file
  std::ifstream ifs; //!< File stream for reading
  std::ifstream::pos_type fileSize = 0; //!< Size of file
  std::string inputContext; //!< Input context
};


#endif // SIMSIMRAPROJECT_H
