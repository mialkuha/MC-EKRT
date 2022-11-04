//Copyright (c) 2021 Mikko Kuha

#ifndef TYPEDEFS_HPP
#define TYPEDEFS_HPP

#ifndef FMGEV
#define FMGEV 5.068
#endif

#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <tuple>

typedef double momentum;
typedef double rapidity;
typedef double xsectval;
typedef double spatial;
typedef int particle_id;

struct dijet_specs
{
    momentum kt{0};
    rapidity y1{0};
    rapidity y2{0};
    particle_id init1{0};
    particle_id init2{0};
    particle_id final1{0};
    particle_id final2{0};
};

struct xsection_id
{
  xsectval sigma;
  particle_id init1;
  particle_id init2;
  particle_id final1;
  particle_id final2;
};

enum class B2_normalization_mode
{
    total,
    inelastic,
    nothing
};

inline std::ostream& operator<<(std::ostream& os, const B2_normalization_mode& mode)
{
    switch (mode)
    {
    case B2_normalization_mode::total:
      os << "total";
      break;
    case B2_normalization_mode::inelastic:
      os << "inelastic";
      break;
    case B2_normalization_mode::nothing:
      os << "nothing";
      break;
    default:
      break;
    }
    return os;
}

struct AA_collision_params
{
  bool mc_glauber_mode;
  bool pp_scattering;
  bool pA_scattering;
  bool spatial_pdfs;
  bool spatial_averaging;
  bool calculate_end_state;
  bool reduce_nucleon_energies;
  xsectval sigma_inel_for_glauber;
  const std::function<spatial(const spatial&)> Tpp;
  B2_normalization_mode normalize_to;
  momentum sqrt_s;
  momentum energy_threshold;
};

struct coords {
 public:
  double x{0};
  double y{0};
  double z{0};
  double mag2() const
  {
    return pow(this->x,2)+pow(this->y,2)+pow(this->z,2);
  }
  double magt2() const
  {
    return pow(this->x,2)+pow(this->y,2);
  }
  coords& operator+=(const coords& rhs)
  {
    this->x += rhs.x;
    this->y += rhs.y;
    this->z += rhs.z;
    return *this;
  }
  
  friend coords operator+(coords lhs, const coords& rhs)
  {
    lhs += rhs;
    return lhs;
  }
  coords& operator-=(const coords& rhs)
  {
    this->x -= rhs.x;
    this->y -= rhs.y;
    this->z -= rhs.z;
    return *this;
  }
  
  friend coords operator-(coords lhs, const coords& rhs)
  {
    lhs -= rhs;
    return lhs;
  }
  coords& operator/=(const double& rhs)
  {
    this->x /= rhs;
    this->y /= rhs;
    this->z /= rhs;
    return *this;
  }
  
  friend coords operator/(coords lhs, const double& rhs)
  {
    lhs /= rhs;
    return lhs;
  }
  coords& operator*=(const double& rhs)
  {
    this->x *= rhs;
    this->y *= rhs;
    this->z *= rhs;
    return *this;
  }
  
  friend coords operator*(coords lhs, const double& rhs)
  {
    lhs *= rhs;
    return lhs;
  }

  friend std::ostream& operator<<(std::ostream& os, const coords& co)
  {
    os << '{' << co.x << ", " << co.y << ", " << co.z<<'}';
    return os;
  }
};

struct dijet_with_coords
{
    dijet_specs dijet;
    coords co;
    double t01;
    double t02;
    double tata;
};

#endif // TYPEDEFS_HPP
