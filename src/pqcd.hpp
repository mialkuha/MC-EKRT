// Copyright (c) 2025 Mikko Kuha (University of Jyväskylä).
// This program is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details. You should have received a copy
// of the GNU General Public License along with this program. If not, see
// <http://www.gnu.org/licenses/>.

#ifndef PQCD_HPP
#define PQCD_HPP

#include <cmath>
#include <functional>
#include <memory>
#include <mutex>
#include <random>
#include <tuple>

#include "cubature.h"
#include "linear_interpolator.hpp"
#include "nn_coll.hpp"
#include "pdf_builder.hpp"
#include "typedefs.hpp"

class nn_coll; // Declaration to get away from the circular references
class pqcd
{
public:
  constexpr static double g_error_tolerance = 1e-4;
  class diff_sigma
  {
  public:
    static auto full_partonic_bookkeeping(
        const std::array<double, 7> &f_i_x1,
        const std::array<double, 7> &f_i_x2,
        const std::array<double, 7> &f_ai_x1,
        const std::array<double, 7> &f_ai_x2,
        const uint_fast8_t &num_flavors,
        const double &s_hat,
        const double &t_hat,
        const double &u_hat

        ) noexcept -> double;

    static auto full_partonic_bookkeeping_1jet(
        const std::array<double, 7> &f_i_x1,
        const std::array<double, 7> &f_i_x2,
        const std::array<double, 7> &f_ai_x1,
        const std::array<double, 7> &f_ai_x2,
        const uint_fast8_t &num_flavors,
        const double &s_hat,
        const double &t_hat,
        const double &u_hat

        ) noexcept -> double;

    static auto sigma_qiqj_qiqj(
        const double &s,
        const double &t,
        const double &u) noexcept -> double;

    static auto sigma_qiqi_qiqi(
        const double &s,
        const double &t,
        const double &u) noexcept -> double;

    static auto sigma_qiaqi_qjaqj(
        const double &s,
        const double &t,
        const double &u) noexcept -> double;

    static auto sigma_qiaqi_qiaqi(
        const double &s,
        const double &t,
        const double &u) noexcept -> double;

    static auto sigma_qiaqi_gg(
        const double &s,
        const double &t,
        const double &u) noexcept -> double;

    static auto sigma_gg_qaq(
        const double &s,
        const double &t,
        const double &u) noexcept -> double;

    static auto sigma_gq_gq(
        const double &s,
        const double &t,
        const double &u) noexcept -> double;

    static auto sigma_gg_gg(
        const double &s,
        const double &t,
        const double &u) noexcept -> double;
  };

  enum scale_choice
  {
    scaled_from_kt,
    constant
  };

  struct nn_coll_params
  {
    double sum_tppa;
    double sum_tppb;
    bool target_neutron;
    bool projectile_neutron;

    explicit nn_coll_params(
        double sum_tppa_ = 0.0,
        double sum_tppb_ = 0.0,
        bool target_neutron_ = false,
        bool projectile_neutron_ = false) noexcept
        : sum_tppa(sum_tppa_),
          sum_tppb(sum_tppb_),
          target_neutron(target_neutron_),
          projectile_neutron(projectile_neutron_)
    {
    }
    nn_coll_params &operator=(const nn_coll_params &rhs)
    {
      this->sum_tppa = rhs.sum_tppa;
      this->sum_tppb = rhs.sum_tppb;
      this->target_neutron = rhs.target_neutron;
      this->projectile_neutron = rhs.projectile_neutron;
      return *this;
    }
  };

  static auto sigma_jet(
      const double &x1,
      const double &x2,
      const double &q2,
      std::shared_ptr<pdf_builder> pdf,
      const double &s_hat,
      const double &t_hat,
      const double &u_hat,
      const double &K_factor,
      pqcd::nn_coll_params &nn_params,
      const bool &average) noexcept -> double;

  struct sigma_jet_params
  {
    scale_choice scale_c;
    double scalar;
    double K_factor;
    bool use_ses;

    explicit sigma_jet_params(
        scale_choice scale_c_ = scaled_from_kt,
        double scalar_ = 1.0,
        double K_factor_ = 1.0,
        bool use_ses_ = false) noexcept
        : scale_c(scale_c_),
          scalar(scalar_),
          K_factor(K_factor_),
          use_ses(use_ses_)
    {
    }
    sigma_jet_params &operator=(const sigma_jet_params &rhs)
    {
      this->scale_c = rhs.scale_c;
      this->scalar = rhs.scalar;
      this->K_factor = rhs.K_factor;
      this->use_ses = rhs.use_ses;
      return *this;
    }
  };

  static auto diff_cross_section_2jet(
      const double &sqrt_s,
      const double &kt,
      const double &y1,
      const double &y2,
      std::shared_ptr<pdf_builder> p_p_pdf,
      pqcd::sigma_jet_params sigma_params,
      pqcd::nn_coll_params &nn_params,
      bool max,
      bool debug // true prints the calculated processes
      ) noexcept -> std::vector<xsection_id>;

  static auto generate_2_to_2_scatt(
      const double &sqrt_s,
      const double &kt_min,
      const double &kt_max,
      std::uniform_real_distribution<double> unirand,
      std::shared_ptr<std::mt19937> eng,
      std::shared_ptr<pdf_builder> p_p_pdf,
      pqcd::sigma_jet_params params,
      const double &power_law,
      envelope_func &env_func,
      pqcd::nn_coll_params &nn_params) noexcept -> dijet_specs;

  static auto sigma_jet_integrand(
      unsigned ndim,
      const double *p_x,
      void *p_fdata,
      unsigned fdim,
      double *p_fval) noexcept -> int;

  static auto sigma_1jet_integrand_binned(
      unsigned ndim,
      const double *p_x,
      void *p_fdata,
      unsigned fdim,
      double *p_fval) noexcept -> int;

  static auto sigma_jet_integrand_binned(
      unsigned ndim,
      const double *p_x,
      void *p_fdata,
      unsigned fdim,
      double *p_fval) noexcept -> int;

  static auto sigma_dijet_integrand_binned(
      unsigned ndim,
      const double *p_x,
      void *p_fdata,
      unsigned fdim,
      double *p_fval) noexcept -> int;

  static auto scale_limits_from_0_1(
      const double &z1,               // variables between 0 and 1
      const double &z2,               // variables between 0 and 1
      const double &z3,               // variables between 0 and 1
      const double &kt2_lower_cutoff, // parameters
      const double &mand_s,           // parameters
      double &kt2,                    // output
      double &y1,                     // output
      double &y2,                     // output
      double &jacobian                // output
      ) noexcept -> void;

  static auto f_ses(
      const double &x,
      const double &q2,
      std::shared_ptr<pdf_builder> p_pdf) noexcept -> double;

  static auto s_hat_from_ys(
      const double &y1,
      const double &y2,
      const double &kt2) noexcept -> double;

  static auto t_hat_from_ys(
      const double &y1,
      const double &y2,
      const double &kt2) noexcept -> double;

  static auto u_hat_from_ys(
      const double &y1,
      const double &y2,
      const double &kt2) noexcept -> double;

  static auto calculate_sigma_jet(
      std::shared_ptr<pdf_builder> p_pdf,
      const double *const p_mand_s,
      const double *const p_kt2_lower_cutoff,
      pqcd::sigma_jet_params params,
      pqcd::nn_coll_params nn_params,
      const bool &average) noexcept -> double;

  static auto calculate_sigma_1jet_binned(
      std::shared_ptr<pdf_builder> p_pdf,
      const double *const p_mand_s,
      const std::tuple<double, double, double, double> *const p_bin,
      pqcd::sigma_jet_params params,
      pqcd::nn_coll_params nn_params) noexcept -> double;

  static auto calculate_sigma_jet_binned(
      std::shared_ptr<pdf_builder> p_pdf,
      const double *const p_mand_s,
      const std::tuple<double, double, double, double> *const p_bin,
      pqcd::sigma_jet_params params,
      pqcd::nn_coll_params nn_params) noexcept -> double;

  static auto calculate_sigma_dijet_binned(
      std::shared_ptr<pdf_builder> p_pdf,
      const double *const p_mand_s,
      const std::tuple<double, double, double, double> *const p_bin,
      pqcd::sigma_jet_params params,
      pqcd::nn_coll_params nn_params) noexcept -> double;

protected:
private:
  std::mutex envelope_maximum_mutex;
};

#endif // PQCD_HPP
