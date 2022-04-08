//Copyright (c) 2021 Mikko Kuha

#ifndef PQCD_HPP
#define PQCD_HPP

#include <cmath>
#include <memory>
#include <mutex>
#include <tuple>

#include <gsl/gsl>

#include "cubature.h"
#include "LHAPDF/GridPDF.h"
#include "nn_coll.hpp"
#include "typedefs.hpp"

class nn_coll;//Declaration to get away from the circular references
class pqcd
{
public:
    constexpr static double g_error_tolerance = 1e-4;
    class diff_sigma
    {
    public:
        struct params
        {
            bool projectile_with_npdfs;
            bool target_with_npdfs;
            bool isoscalar_projectile;
            bool isoscalar_target;
            int npdf_setnumber;

            explicit params
            (
              bool projectile_with_npdfs_ = false, 
              bool target_with_npdfs_     = false, 
              bool isoscalar_projectile_  = false, 
              bool isoscalar_target_      = false, 
              int npdf_setnumber_         = 1
            ) noexcept
            : projectile_with_npdfs(projectile_with_npdfs_), 
              target_with_npdfs(target_with_npdfs_), 
              isoscalar_projectile(isoscalar_projectile_), 
              isoscalar_target(isoscalar_target_), 
              npdf_setnumber(npdf_setnumber_) 
            {}
        };
        
        static auto sigma_qiqj_qiqj
        (
          const momentum &s, 
          const momentum &t, 
          const momentum &u, 
          const momentum &alpha_s
        ) noexcept -> xsectval;

        static auto sigma_qiqi_qiqi
        (
          const momentum &s, 
          const momentum &t, 
          const momentum &u, 
          const momentum &alpha_s
        ) noexcept -> xsectval;

        static auto sigma_qiaqi_qjaqj
        (
          const momentum &s, 
          const momentum &t, 
          const momentum &u, 
          const momentum &alpha_s
        ) noexcept -> xsectval;

        static auto sigma_qiaqi_qiaqi
        (
          const momentum &s, 
          const momentum &t, 
          const momentum &u, 
          const momentum &alpha_s
        ) noexcept -> xsectval;

        static auto sigma_qiaqi_gg
        (
          const momentum &s, 
          const momentum &t, 
          const momentum &u, 
          const momentum &alpha_s
        ) noexcept -> xsectval;

        static auto sigma_gg_qaq
        (
          const momentum &s, 
          const momentum &t, 
          const momentum &u, 
          const momentum &alpha_s
        ) noexcept -> xsectval;

        static auto sigma_gq_gq
        (
          const momentum &s, 
          const momentum &t, 
          const momentum &u, 
          const momentum &alpha_s
        ) noexcept -> xsectval;

        static auto sigma_gg_gg
        (
          const momentum &s, 
          const momentum &t, 
          const momentum &u, 
          const momentum &alpha_s
        ) noexcept -> xsectval;

        static auto sigma_jet
        (
          const rapidity &x1,    
          const rapidity &x2, 
          const momentum &q2,    
          std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
          const momentum &s_hat, 
          const momentum &t_hat, 
          const momentum &u_hat, 
          const params *const p_params,
          std::shared_ptr<LHAPDF::GridPDF> p_n_pdf
        ) noexcept -> xsectval;

        static auto spatial_sigma_jet_mf
        (
          const rapidity &x1,    
          const rapidity &x2, 
          const momentum &q2,    
          std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
          const momentum &s_hat, 
          const momentum &t_hat, 
          const momentum &u_hat, 
          const pqcd::diff_sigma::params *const p_params, 
          std::shared_ptr<LHAPDF::GridPDF> p_n_pdf, 
          std::function<double(double const&)> rA_spatial,
          std::function<double(double const&)> rB_spatial
        ) noexcept -> xsectval;

        static auto spatial_sigma_jet_full
        (
          const rapidity &x1,    
          const rapidity &x2, 
          const momentum &q2,    
          std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
          const momentum &s_hat, 
          const momentum &t_hat, 
          const momentum &u_hat, 
          const pqcd::diff_sigma::params *const p_params, 
          std::shared_ptr<LHAPDF::GridPDF> p_n_pdf, 
          std::function<double(double const&)> rA_spatial,
          std::function<double(double const&)> rB_spatial,
          const std::array<const double, 3> &T_sums
        ) noexcept -> xsectval;
    };

    enum scale_choice
    {
      scaled_from_kt,
      constant
    };

    struct sigma_jet_params
    {
        scale_choice scale_c;
        momentum scalar;
        diff_sigma::params *p_d_params;
        bool use_ses;

        explicit sigma_jet_params
          (
            scale_choice scale_c_       = scaled_from_kt, 
            momentum scalar_            = 1.0, 
            bool projectile_with_npdfs_ = false, 
            bool target_with_npdfs_     = false, 
            bool isoscalar_projectile_  = false, 
            bool isoscalar_target_      = false, 
            int npdf_setnumber_         = 1, 
            bool use_ses_               = false
          ) noexcept
          : scale_c(scale_c_), 
            scalar(scalar_), 
            p_d_params(new diff_sigma::params
                        (
                          projectile_with_npdfs_, 
                          target_with_npdfs_, 
                          isoscalar_projectile_, 
                          isoscalar_target_, 
                          npdf_setnumber_
                        )), 
            use_ses(use_ses_)
          {}

        explicit sigma_jet_params
          (
            diff_sigma::params *p_d_params_, 
            scale_choice scale_c_             = scaled_from_kt, 
            momentum scalar_                  = 1.0, 
            bool use_ses_                     = false
          ) noexcept
          : scale_c(scale_c_), 
            scalar(scalar_), 
            p_d_params(p_d_params_), 
            use_ses(use_ses_)
          {}

    };

    static void generate_bin_NN_coll(nn_coll * coll) noexcept;

    static auto sigma_jet_integrand
      (
        unsigned ndim, 
        const double *p_x, 
        void *p_fdata, 
        unsigned fdim, 
        double *p_fval
      ) noexcept -> int;

    static void scale_limits_from_0_1
      (
        const rapidity &z1, //variables between 0 and 1
        const rapidity &z2, //variables between 0 and 1
        const rapidity &z3, //variables between 0 and 1
        const momentum &kt2_lower_cutoff, //parameters
        const momentum &mand_s, //parameters
        momentum &kt2, //output
        rapidity &y1, //output
        rapidity &y2, //output
        xsectval &jacobian //output
      ) noexcept;

    static auto spatial_sigma_jet_integrand_mf
      (
        unsigned ndim, 
        const double *p_x, 
        void *p_fdata, 
        unsigned fdim, 
        double *p_fval
      ) noexcept -> int;

    static auto spatial_sigma_jet_integrand_full
      (
        unsigned ndim, 
        const double *p_x, 
        void *p_fdata, 
        unsigned fdim, 
        double *p_fval
      ) noexcept -> int;

    static auto spatial_sigma_jet_integrand_factored
      (
        unsigned ndim, 
        const double *p_x, 
        void *p_fdata, 
        unsigned fdim, 
        double *p_fval
      ) noexcept -> int;

    static auto f_ses
      (
        const rapidity &x, 
        const momentum &q2, 
        std::shared_ptr<LHAPDF::GridPDF> p_pdf
      ) noexcept -> xsectval;

    static auto s_hat_from_ys
      (
        const rapidity &y1,
        const rapidity &y2, 
        const momentum &kt2
      ) noexcept -> momentum;

    static auto t_hat_from_ys
      (
        const rapidity &y1,
        const rapidity &y2, 
        const momentum &kt2
      ) noexcept -> momentum;

    static auto u_hat_from_ys
      (
        const rapidity &y1,
        const rapidity &y2, 
        const momentum &kt2
      ) noexcept -> momentum;


    static auto calculate_sigma_jet
      (
        std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
        const momentum *const p_mand_s, 
        const momentum *const p_kt2_lower_cutoff, 
        const pqcd::sigma_jet_params *const p_params
      ) noexcept -> xsectval;

    static auto calculate_spatial_sigma_jet_mf
      (
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
        std::shared_ptr<LHAPDF::GridPDF> p_n_pdf, 
        const momentum *const p_mand_s, 
        const momentum *const p_kt2_lower_cutoff, 
        const pqcd::sigma_jet_params *const p_params, 
        const spatial &sum_tppa, 
        const spatial &sum_tppb, 
        const spatial &tAA_0, 
        const spatial &tBB_0
      ) noexcept -> xsectval;

    static auto calculate_spatial_sigma_jet_full
      (
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
        std::shared_ptr<LHAPDF::GridPDF> p_n_pdf, 
        const momentum *const p_mand_s, 
        const momentum *const p_kt2_lower_cutoff, 
        const pqcd::sigma_jet_params *const p_params, 
        const std::array<const double, 3> &T_sums,
        const spatial &tAA_0, 
        const spatial &tBB_0
      ) noexcept -> xsectval;

    static auto calculate_spatial_sigma_jet_factored
      (
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
        std::shared_ptr<LHAPDF::GridPDF> p_n_pdf, 
        const momentum *const p_mand_s, 
        const momentum *const p_kt2_lower_cutoff, 
        const pqcd::sigma_jet_params *const p_params
      ) noexcept -> std::array<xsectval,4>;

protected:
private:
};

#endif // PQCD_HPP
