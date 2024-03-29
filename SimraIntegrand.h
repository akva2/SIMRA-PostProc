// $Id$
//==============================================================================
//!
//! \file SimraIntegrand.h
//!
//! \date Jul 6 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations simra result projection.
//!
//==============================================================================

#ifndef _SIMRA_INTEGRAND_H_
#define _SIMRA_INTEGRAND_H_

#include "IntegrandBase.h"
#include "MatVec.h"
#include "Vec3.h"

#include <cstddef>
#include <string>
#include <vector>


class AnaSol;
class FiniteElement;
class LocalIntegral;
class Tensor;


/*!
  \brief Class representing the integrand of the Darcy problem.
*/

class SimraIntegrand : public IntegrandBase
{
public:
  //! \brief The constructor initializes all pointers to zero.
  explicit SimraIntegrand();

  //! \brief Empty destructor.
  virtual ~SimraIntegrand() {}

  using IntegrandBase::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  //!
  //! \details This method is invoked once before starting the numerical
  //! integration loop over the Gaussian quadrature points over an element.
  //! It is supposed to perform all the necessary internal initializations
  //! needed before the numerical integration is started for current element.
  bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt) override;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  bool evalSol(Vector& s, const FiniteElement& fe, const Vec3& X,
               const std::vector<int>& MNPC) const override;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution (optional)
  NormBase* getNormIntegrand(AnaSol* asol) const override;

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  size_t getNoFields(int group = 0) const override;
  //! \brief Returns the name of a primary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  std::string getField1Name(size_t i, const char* prefix = 0) const override;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  std::string getField2Name(size_t i, const char* prefix = 0) const override;

  //! \brief Calculates the velocity at current point.
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] vec Primary solution vector for current element
  Vec3 velocity(const FiniteElement& fe, const Vectors& vec) const;

  //! \brief Calculates the velocity gradient at current point.
  //! \param[in] fe Finite element data of current integration point
  //! \param[out] grad Velocity gradient tensor at current point
  //! \param[in] vec Primary solution vector for current element
  void velocityGradient(const FiniteElement& fe,
                        Tensor& grad, const Vectors& vec) const;

  //! \brief Symmetrize the supplied velocity gradient if necessary.
  //! \param[out] eps Strain tensor at current point
  void strain(Tensor& eps) const;

  //! \brief Calculates the (Cauchy) stress tensor at current point.
  //! \param[in] fe Finite element data of current integration point
  //! \param[out] sigma Strain tensor at current point
  //! \param[in] vec Primary solution vector for current element
  void stress(const FiniteElement& fe,
              Tensor& sigma, const Vectors& vec) const;

  //! \brief Calculates the pressure at current point.
  //! \param[in] fe Finite element data of current integration point
  double pressure(const FiniteElement& fe) const;

  //! \brief Calculates the viscosity at current point.
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] vec Primary solution vector for current element
  double viscosity(const FiniteElement& fe, const Vectors& vec) const;

  //! \brief Calculates the turbulent kinetic energy at current point.
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] vec Primary solution vector for current element
  double TKE(const FiniteElement& fe, const Vectors& vec) const;

  //! \brief Calculates the turbulent dissipation rate at current point.
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] vec Primary solution vector for current element
  double TD(const FiniteElement& fe, const Vectors& vec) const;

  //! \brief Calculates the temperature at current point.
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] vec Primary solution vector for current element
  double temperature(const FiniteElement& fe, const Vectors& vec) const;

  //! \brief Calculates the temperature gradient at current point.
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] vec Primary solution vector for current element
  Vec3 temperatureGradient(const FiniteElement& fe,
                           const Vectors& vec) const;

  Vector elmPressure; //!< Element-wise pressures
  Vector dist; //!< Wall distance
  double nu = 1.0; //!< Fluid viscosity
};


/*!
  \brief Class representing the integrand of Darcy energy norms.
*/

class SimraNorm : public NormBase
{
public:
  //! \brief Enumeration of primary norms.
  enum PrimaryNormEntries {
    L2_Uh     = 0,
    H1_Uh,
    L2_DIV_Uh,
    L2_Ph,
    L2_SIGMAh,
    L2_pTh,
    H1_pTh,
    L2_TK,
    L2_TD,
    L2_U,
    L2_E_U,
    H1_U,
    H1_E_U,
    L2_P,
    L2_E_P,
    L2_SIGMA,
    L2_E_SIGMA,
    L2_pT,
    L2_E_pT,
    H1_pT,
    H1_E_pT,
    TOTAL_ERROR
  };

  //! \brief Enumeration of recovery norms.
  enum RecoveryNormEntries {
    H1_Ur_Uh = 0,
    TOTAL_NORM_REC,
    L2_Pr_Ph,
    L2_DIV_Ur,
    L2_SIGMA_REC,
    H1_pTr_pTh,
    H1_Ur_U,
    TOTAL_E_REC,
    L2_Pr_P,
    L2_SIGMA_E_REC,
    H1_pTr_pT,
    EFF_REC_VEL,
    EFF_REC_PRESS,
    EFF_REC_STRESS,
    EFF_REC_TEMP,
    EFF_REC_TOTAL
  };

  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Poisson problem to evaluate norms for
  //! \param[in] aSol The analytical solution (optional)
  SimraNorm(SimraIntegrand& p, AnaSol* aSol);
  //! \brief Empty destructor.
  virtual ~SimraNorm() {}

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X) const override;

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  size_t getNoFields(int group) const override;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  std::string getName(size_t i, size_t j, const char* prefix) const override;

protected:
  AnaSol* aSol; //!< Analytical solution
};

#endif
