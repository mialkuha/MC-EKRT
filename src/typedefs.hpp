//Copyright (c) 2023 Mikko Kuha

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

typedef int_fast8_t particle_id;

struct dijet_specs
{
    double kt{0};
    double y1{0};
    double y2{0};
    particle_id init1{0};
    particle_id init2{0};
    particle_id final1{0};
    particle_id final2{0};
    double pro_pdf{0};
    double tar_pdf{0};
};

struct envelope_func
{
    double min_kt{0};
    double norm1{0};
    double norm2{0};
    double power{0};
    double switch_kt{0};
    double prim_integ_constant{0};
    double prim_switch_y{0};
    const std::function<double(const double&)> func;
    const std::function<double(const double&)> prim;
    const std::function<double(const double&)> prim_inv;
};

struct xsection_id
{
    double sigma;
    particle_id init1;
    particle_id init2;
    particle_id final1;
    particle_id final2;
    double pro_pdf{0};
    double tar_pdf{0};
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

inline bool almostEqual(double x, double y)    
{    
    double maxXYOne = std::max( { 1.0, std::abs(x) , std::abs(y) } ) ;
    return std::abs(x - y) <= std::numeric_limits<double>::epsilon()*maxXYOne ;
}

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
    os << co.x << "," << co.y << "," << co.z;
    return os;
  }

  friend std::ostream& operator==(std::ostream& os, const coords& co)
  {
    os << co.x << "," << co.y << "," << co.z;
    return os;
  }

  friend bool operator==(const coords &lhs, const coords &rhs)
  {
    return almostEqual(lhs.x, rhs.x)&&
           almostEqual(lhs.y, rhs.y)&&
           almostEqual(lhs.z, rhs.z);
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

struct hotspot_info
{
    coords co{0,0,0};
};

#endif // TYPEDEFS_HPP
