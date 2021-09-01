//==============================================================================
//!
//! \file SimraFieldGenerator.C
//!
//! \date Sep 1 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Generators for SIMRA solution fields.
//!
//==============================================================================

#include "SimraFieldGenerator.h"

#include "ASMs3DSimra.h"
#include "SimraEnums.h"

#include "IFEM.h"
#include "Utilities.h"

#include "tinyxml.h"


SimraFieldGenerator::SimraFieldGenerator (Vectors& solutions) :
  solution(solutions)
{
}


void SimraFieldGenerator::setPatch (const ASMs3DSimra *patch)
{
  pch = patch;
  pch->getNoStructNodes(n[0],n[1],n[2]);
}


void SimraFieldGenerator::read (const TiXmlElement* elem)
{
  utl::getAttribute(elem, "bl_height", bl_height);
  utl::getAttribute(elem, "wind_speed", wind_speed);
  utl::getAttribute(elem, "wind_direction", wind_dir);
  utl::getAttribute(elem, "surface_roughness", z0);
  utl::getAttribute(elem, "velocity", uni_velocity);
}


void SimraFieldGenerator::print ()
{
  IFEM::cout << "Field generator settings:\n"
             << "\t Boundary layer height = " << bl_height
             << "\n\t Wind speed at boundary layer = " << wind_speed
             << "\n\t Wind direction = " << wind_dir
             << "\n\t Surface roughness = " << z0
             << "\n\t Generate velocity = " << (uni_velocity ? "true" : "false")
             << std::endl;
}


void SimraFieldGenerator::hydrostaticConditions ()
{
  auto T_prof = potentialTemperature();
  auto psrho = hydrostaticPressureDensity(T_prof);

  // initialize field with potential temperature
  size_t idx = 0;
  for (size_t k = 0; k < n[2]; ++k)
    for (size_t j = 0; j < n[1]; ++j)
      for (size_t i = 0; i < n[0]; ++i, ++idx)
        solution[PT][idx] = T_prof[k];

  // set hydrostatic pressure and density
  idx = 0;
  const auto& ps1 = psrho[0];
  const auto& rhos1 = psrho[1];
  const auto& z_ref = psrho[2];

  for (size_t k = 0; k < n[2]; ++k)
    for (size_t j = 0; j < n[1]; ++j)
      for (size_t i = 0; i < n[0]; ++i, ++idx) {
        Vec3 coord = pch->getCoord(idx+1);
        auto it = std::lower_bound(z_ref.begin(), z_ref.end(), coord[2]);
        if (it == z_ref.end()) {
          solution[PS][idx] = ps1.back();
          solution[RHOS][idx] = rhos1.back();
        } else if (it == z_ref.begin()) {
          solution[PS][idx] = ps1.front();
          solution[RHOS][idx] = rhos1.front();
        } else {
          size_t kidx = it - z_ref.begin();
          double den = z_ref[kidx] - z_ref[kidx-1];
          double alpha = coord[2] - z_ref[kidx-1];
          solution[PS][idx] = ps1[kidx-1] + (ps1[kidx] - ps1[kidx-1]) * alpha / den;
          solution[RHOS][idx] = rhos1[kidx-1] + (rhos1[kidx] - rhos1[kidx-1]) * alpha / den;
        }
      }
  solution[RHO] = solution[RHOS];
  solution[PTS] = solution[PT];
}


void SimraFieldGenerator::uniformVelocity ()
{
  double cappa=0.42;
  double u_a = wind_speed;//  ! wind at delta
  double alpha = wind_dir / 180 * M_PI;

  double delta = bl_height;
  delta += z0;

  double eta = 1.0;
  double wake = 3.0*eta - 2.0*eta*eta;
  double ustar = u_a*cappa / (std::log(delta/z0)+wake);

  std::vector<double> pro_u(n[2]), pro_v(n[2]), pro_w(n[2]);

  Vec3 z_bakke = pch->getCoord(1);
  for (size_t k = 0; k < n[2]; ++k) {
    Vec3 coord = pch->getCoord(1 + k*n[0]*n[1]);
    double z = coord[2]-z_bakke[2]+z0;
    eta = std::min(z/delta, 1.0);
    wake = 3.0*eta - 2.0*eta*eta;
    double profu = std::min((ustar/cappa)*(std::log(z/z0)+wake), u_a);
    pro_u[k] = profu*cos(alpha);
    pro_v[k] = profu*sin(alpha);
    pro_w[k] = 0.0;
  }

  size_t idx = 0;
  for (size_t k = 0; k < n[2]; ++k)
    for (size_t j = 0; j < n[1]; ++j)
      for (size_t i = 0; i < n[0]; ++i, ++idx) {
        solution[U_X][idx] = pro_u[k];
        solution[U_Y][idx] = pro_v[k];
        solution[U_Z][idx] = pro_w[k];
      }
}


void SimraFieldGenerator::turbulenceFields ()
{
  static constexpr double av0 = 1.0e-5;
  static constexpr double rho0 = 1.3;
  static constexpr double cm34 = 0.164316773;

  double delta = bl_height;
  delta += z0;

  std::vector<double> z_ref(n[2]), rho(n[2]),
                      pro_u(n[2]), pro_v(n[2]), pro_w(n[2]);
  for (size_t i = 0; i < n[0]; ++i)
    for (size_t j = 0; j < n[1]; ++j) {
      auto&& idx = [i,j,this](size_t k) { return i + (j + k*n[1])*n[0]; };
      for (size_t k = 0; k < n[2]; ++k) {
        Vec3 coord = pch->getCoord(1 + k*n[0]*n[1]);
        z_ref[k] = coord[2];
        pro_u[k] = solution[U_X][idx(k)];
        pro_v[k] = solution[U_Y][idx(k)];
        rho[k] = solution[RHO][idx(k)];
      }

    std::vector<double> av(n[2]), td(n[2]), tk(n[2]);
    for (size_t k = 0; k < n[2]; ++k) {
      size_t kp = k < n[2]-1 ? k+1 : n[2]-1;
      size_t km = k > 0 ? k-1 : 0;
      double z1 = z_ref[k]-z_ref[0] + z0;
      double al0 = 0.42*z1/(1.0+4*z1/delta);
      double dz = z_ref[kp]-z_ref[km];
      double drdz = (rho[kp]-rho[km]) / dz;
      double any2 = -(9.81/rho0)*drdz;
      double dudz = (pro_u[kp]-pro_u[km]) / dz;
      double dvdz = (pro_v[kp]-pro_v[km]) / dz;
      double prod = dudz*dudz + dvdz*dvdz + 1.0e-8;
      double Ri = any2 / prod;
      double fact = std::max(1.0-0.5*Ri, 0.00001);
      av[k] = std::max(al0*al0*sqrt(prod)*sqrt(fact), av0);
      double u_top = hypot(pro_u.back(), pro_v.back());
      double tkmin = 0.001*u_top*u_top;
      tk[k] = std::max(av[k]*av[k]/(0.3*al0*al0), tkmin);
      td[k] = cm34*pow(tk[k],1.5)/al0;
    }

    for (size_t k = 0; k < n[2]; ++k) {
      solution[VTEF][idx(k)] = av[k];
      solution[TK][idx(k)] = tk[k];
      solution[TD][idx(k)] = td[k];
    }
  }
}


std::vector<double> SimraFieldGenerator::potentialTemperature () const
{
  auto&& teta_profile = [](double z)
  {
    static constexpr double T0 = 272.3;
    static constexpr double z1 = 800.0;
    static constexpr double a1 = 0.0007;
    static constexpr double a12 = 0.009;
    static constexpr double dz12 = 1100.0;
    static constexpr double dz23 = 700.0;

    double T = T0;
    if (z <= z1)
      T = T0 + a1*z;
    double T1 = T0 + a1*z1;
    double z2 = z1 + dz12;
    if (z > z1 && z <= z2)
      T1 = T1 + a12*(z - z1);
    double T2 = T1 + a12*dz12;
    double a23 = a12;
    double z3 = z2 + dz23;
    if (z > z2 && z <= z3)
      T = T2 + a23*(z - z2);
    double T3 = T2 + a23*dz23;
    double a3 = a1;
    if (z > z3)
      T = T3 + a3*(z - z3);

    return T;
  };

  std::vector<double> T_prof(n[2]);

  Vec3 z_bakke = pch->getCoord(1);
  for (size_t k = 0; k < n[2]; ++k) {
    Vec3 coord = pch->getCoord(1 + k*n[0]*n[1]);
    T_prof[k] = teta_profile(coord[2]-z_bakke[2]);
  }

  return T_prof;
}


std::array<std::vector<double>,3>
SimraFieldGenerator::hydrostaticPressureDensity (const std::vector<double>& T_prof) const
{
  static constexpr double grav = 9.81;
  static constexpr double p_ref = 100000.0;  // static pressure at sea level [Pa]
  static constexpr double R_gas = 287.0; // Gas constant
  static constexpr double gamma = 2.0 / 7.0;
  std::array<std::vector<double>,3> result;
  result[0].resize(n[2]);
  result[1].resize(n[2]);
  result[2].resize(n[2]);

  auto& ps1 = result[0];
  auto& rhos1 = result[1];
  auto& z_ref = result[2];
  for (size_t k = 1; k < n[2]; ++k) {
    Vec3 coord = pch->getCoord(1 + k*n[0]*n[1]);
    z_ref[k] = coord[2];
  }

  z_ref[0] = 0.0;
  result[0][0] = p_ref;
  result[1][0] = p_ref / (R_gas * T_prof[0]);

  for (size_t it = 0; it < 2; ++it)
    for (size_t k = 1; k < n[2]; ++k) {
      double dz = z_ref[k] - z_ref[k-1];
      double rhos_avr = it == 0 ? rhos1[k-1] : 0.5*(rhos1[k-1] + rhos1[k]);
      ps1[k] =  ps1[k-1] - grav*rhos_avr*dz;
      rhos1[k] = ps1[k] / (R_gas*T_prof[k]) * pow(p_ref/ps1[k], gamma);
    }

  return result;
}
