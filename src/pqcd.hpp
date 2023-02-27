//Copyright (c) 2022 Mikko Kuha

#ifndef PQCD_HPP
#define PQCD_HPP

#include <cmath>
#include <functional>
#include <memory>
#include <mutex>
#include <random>
#include <tuple>


#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshorten-64-to-32"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "LHAPDF/GridPDF.h"
#pragma GCC diagnostic pop

#include "cubature.h"
#include "linear_interpolator.hpp"
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
            bool projectile_with_npdfs{false};
            bool target_with_npdfs{false};
            bool isoscalar_projectile{false};
            bool isoscalar_target{false};
            bool npdfs_spatial{false};
            int npdf_setnumber{1};
            double spatial_cutoff{0.0};
            double K_factor{1};
            uint_fast16_t A{1};
            uint_fast16_t B{1};
            uint_fast16_t ZA{1};
            uint_fast16_t ZB{1};
            std::shared_ptr<LHAPDF::GridPDF> p_n_pdf{nullptr};
            std::function<double(double const&)> rA_spatial{nullptr};
            std::function<double(double const&)> rB_spatial{nullptr};

            explicit params
            (
              auto projectile_with_npdfs_, 
              auto target_with_npdfs_,     
              auto isoscalar_projectile_,  
              auto isoscalar_target_,
              auto npdfs_spatial_,
              auto npdf_setnumber_,        
              auto spatial_cutoff_,
              auto K_factor_,
              auto A_,                     
              auto B_,                     
              auto ZA_,                    
              auto ZB_,                    
              auto p_n_pdf_,               
              auto rA_spatial_,            
              auto rB_spatial_            
            ) noexcept
            : projectile_with_npdfs(projectile_with_npdfs_), 
              target_with_npdfs(target_with_npdfs_), 
              isoscalar_projectile(isoscalar_projectile_), 
              isoscalar_target(isoscalar_target_), 
              npdfs_spatial(npdfs_spatial_), 
              npdf_setnumber(npdf_setnumber_),
              spatial_cutoff(spatial_cutoff_),
              K_factor(K_factor_), 
              A(A_),
              B(B_),
              ZA(ZA_),
              ZB(ZB_),
              p_n_pdf(p_n_pdf_),
              rA_spatial(std::move(rA_spatial_)),
              rB_spatial(std::move(rB_spatial_))
            {}

            explicit params
            (
              auto projectile_with_npdfs_, 
              auto target_with_npdfs_,     
              auto isoscalar_projectile_,  
              auto isoscalar_target_,  
              auto npdfs_spatial_, 
              auto npdf_setnumber_,  
              auto spatial_cutoff_,     
              auto K_factor_,      
              auto A_,                     
              auto B_,                     
              auto ZA_,                    
              auto ZB_
            ) noexcept
            : projectile_with_npdfs(projectile_with_npdfs_), 
              target_with_npdfs(target_with_npdfs_), 
              isoscalar_projectile(isoscalar_projectile_), 
              isoscalar_target(isoscalar_target_), 
              npdfs_spatial(npdfs_spatial_), 
              npdf_setnumber(npdf_setnumber_),
              spatial_cutoff(spatial_cutoff_),
              K_factor(K_factor_), 
              A(A_),
              B(B_),
              ZA(ZA_),
              ZB(ZB_)
            {}
            params& operator=(const params& rhs)
            {
              this->projectile_with_npdfs = rhs.projectile_with_npdfs;
              this->target_with_npdfs = rhs.target_with_npdfs;
              this->isoscalar_projectile = rhs.isoscalar_projectile;
              this->isoscalar_target = rhs.isoscalar_target;
              this->npdfs_spatial = rhs.npdfs_spatial;
              this->npdf_setnumber = rhs.npdf_setnumber;
              this->spatial_cutoff = rhs.spatial_cutoff;
              this->K_factor = rhs.K_factor;
              this->A = rhs.A;
              this->B = rhs.B;
              this->ZA = rhs.ZA;
              this->ZB = rhs.ZB;
              this->p_n_pdf = rhs.p_n_pdf;
              this->rA_spatial = rhs.rA_spatial;
              this->rB_spatial = rhs.rB_spatial;
              return *this;
            }
        };

        static auto make_pdfs
        (
          const double &x1, 
          const double &x2, 
          const double &q2, 
          std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
          /*std::shared_ptr<LHAPDF::GridPDF> p_n_pdf, */
          const bool &projectile_with_npdfs,
          const bool &target_with_npdfs,
          const bool &isoscalar_projectile,
          const bool &isoscalar_target,
          const bool &npdfs_spatial,
          const int &npdf_setnumber,
          const uint_fast16_t &A,
          const uint_fast16_t &B,
          const uint_fast16_t &ZA,
          const uint_fast16_t &ZB,
          const std::function<double (const double &)> &rA_spatial,
          const std::function<double (const double &)> &rB_spatial
        ) noexcept -> std::tuple
          <
            std::array<double, 7>,
            std::array<double, 7>,
            std::array<double, 7>,
            std::array<double, 7>
          >;

        static auto full_partonic_bookkeeping
        (
            const std::array<double, 7> &f_i_x1, 
            const std::array<double, 7> &f_i_x2,
            const std::array<double, 7> &f_ai_x1,
            const std::array<double, 7> &f_ai_x2,
            const uint_fast8_t &num_flavors,
            const double &s_hat,
            const double &t_hat, 
            const double &u_hat

        ) noexcept -> double;

        static auto full_partonic_bookkeeping_1jet
        (
            const std::array<double, 7> &f_i_x1, 
            const std::array<double, 7> &f_i_x2,
            const std::array<double, 7> &f_ai_x1,
            const std::array<double, 7> &f_ai_x2,
            const uint_fast8_t &num_flavors,
            const double &s_hat,
            const double &t_hat, 
            const double &u_hat

        ) noexcept -> double;
        
        static auto sigma_qiqj_qiqj
        (
          const double &s, 
          const double &t, 
          const double &u
        ) noexcept -> double;

        static auto sigma_qiqi_qiqi
        (
          const double &s, 
          const double &t, 
          const double &u
        ) noexcept -> double;

        static auto sigma_qiaqi_qjaqj
        (
          const double &s, 
          const double &t, 
          const double &u
        ) noexcept -> double;

        static auto sigma_qiaqi_qiaqi
        (
          const double &s, 
          const double &t, 
          const double &u
        ) noexcept -> double;

        static auto sigma_qiaqi_gg
        (
          const double &s, 
          const double &t, 
          const double &u
        ) noexcept -> double;

        static auto sigma_gg_qaq
        (
          const double &s, 
          const double &t, 
          const double &u
        ) noexcept -> double;

        static auto sigma_gq_gq
        (
          const double &s, 
          const double &t, 
          const double &u
        ) noexcept -> double;

        static auto sigma_gg_gg
        (
          const double &s, 
          const double &t, 
          const double &u
        ) noexcept -> double;

        static auto sigma_1jet
        (
          const double &x1,    
          const double &x2, 
          const double &q2,    
          std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
          const double &s_hat, 
          const double &t_hat, 
          const double &u_hat, 
          const params p_params/*,
          std::shared_ptr<LHAPDF::GridPDF> p_n_pdf*/
        ) noexcept -> double;

        static auto sigma_jet
        (
          const double &x1,    
          const double &x2, 
          const double &q2,    
          std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
          const double &s_hat, 
          const double &t_hat, 
          const double &u_hat, 
          const params p_params/*,
          std::shared_ptr<LHAPDF::GridPDF> p_n_pdf*/
        ) noexcept -> double;
    };

    enum scale_choice
    {
      scaled_from_kt,
      constant
    };

    struct sigma_jet_params
    {
        scale_choice scale_c;
        double scalar;
        diff_sigma::params d_params;
        bool use_ses;
        bool snPDFs_linear;
        linear_interpolator c_A_func;

        explicit sigma_jet_params
          (
            diff_sigma::params d_params_, 
            scale_choice scale_c_    = scaled_from_kt, 
            double scalar_           = 1.0,
            bool use_ses_            = false,
            bool snPDFs_linear_      = true,
            linear_interpolator c_A_func_ = linear_interpolator(std::vector<double>{0.0},std::vector<double>{0.0})
          ) noexcept
          : scale_c(scale_c_), 
            scalar(scalar_), 
            d_params(d_params_), 
            use_ses(use_ses_), 
            snPDFs_linear(snPDFs_linear_), 
            c_A_func(c_A_func_)
          {}
        sigma_jet_params& operator=(const sigma_jet_params& rhs)
        {
          this->scale_c       = rhs.scale_c;
          this->scalar        = rhs.scalar;
          this->d_params      = rhs.d_params;
          this->use_ses       = rhs.use_ses;
          this->snPDFs_linear = rhs.snPDFs_linear;
          this->c_A_func      = rhs.c_A_func;
          return *this;
        }
    };

    static auto diff_cross_section_2jet
    (
      const double &sqrt_s,
      const double &kt, 
      const double &y1, 
      const double &y2,
      std::shared_ptr<LHAPDF::GridPDF> p_p_pdf,
      pqcd::sigma_jet_params sigma_params,
      bool debug = false //true prints the calculated processes
    ) noexcept -> std::vector<xsection_id>;

    static auto throw_0_truncated_poissonian
      (
        const double &lambda, 
        std::uniform_real_distribution<double> unirand, 
        std::shared_ptr<std::mt19937> eng
      ) noexcept -> uint_fast8_t;

    static auto generate_2_to_2_scatt
      (
        const double &sqrt_s,
        const double &kt_min,
        const double &kt_max,
        std::uniform_real_distribution<double> unirand, 
        std::shared_ptr<std::mt19937> eng,
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf,
        pqcd::sigma_jet_params params,
        const double &power_law,
        double &envelope_maximum
      ) noexcept -> dijet_specs;

    static auto generate_2_to_2_scatt
      (
        const double &sqrt_s,
        const double &kt_min,
        const double &kt_max,
        std::uniform_real_distribution<double> unirand, 
        std::shared_ptr<std::mt19937> eng,
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf,
        pqcd::sigma_jet_params params,
        const double &power_law,
        envelope_func &env_func
      ) noexcept -> dijet_specs;

    static auto generate_bin_NN_coll
      (
        nn_coll &coll,
        const double &sigma_jet,
        const double &Tpp_b, 
        const double &kt0,
        std::uniform_real_distribution<double> unirand, 
        std::shared_ptr<std::mt19937> eng,
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf,
        pqcd::sigma_jet_params params,
        const double &power_law,
        double &envelope_maximum
      ) noexcept -> void;

    static auto generate_bin_NN_coll
      (
        nn_coll &coll,
        const double &sigma_jet,
        const double &Tpp_b, 
        const double &kt0,
        std::uniform_real_distribution<double> unirand, 
        std::shared_ptr<std::mt19937> eng,
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf,
        pqcd::sigma_jet_params params,
        const double &power_law,
        envelope_func &env_func
      ) noexcept -> void;

    static auto sigma_jet_integrand
      (
        unsigned ndim, 
        const double *p_x, 
        void *p_fdata, 
        unsigned fdim, 
        double *p_fval
      ) noexcept -> int;

    static auto sigma_1jet_integrand_binned
      (
        unsigned ndim, 
        const double *p_x, 
        void *p_fdata, 
        unsigned fdim, 
        double *p_fval
      ) noexcept -> int;

    static auto sigma_jet_integrand_binned
      (
        unsigned ndim, 
        const double *p_x, 
        void *p_fdata, 
        unsigned fdim, 
        double *p_fval
      ) noexcept -> int;

    static auto sigma_dijet_integrand_binned
      (
        unsigned ndim, 
        const double *p_x, 
        void *p_fdata, 
        unsigned fdim, 
        double *p_fval
      ) noexcept -> int;

    static auto scale_limits_from_0_1
      (
        const double &z1, //variables between 0 and 1
        const double &z2, //variables between 0 and 1
        const double &z3, //variables between 0 and 1
        const double &kt2_lower_cutoff, //parameters
        const double &mand_s, //parameters
        double &kt2, //output
        double &y1, //output
        double &y2, //output
        double &jacobian //output
      ) noexcept -> void;

    static auto f_ses
      (
        const double &x, 
        const double &q2, 
        std::shared_ptr<LHAPDF::GridPDF> p_pdf
      ) noexcept -> double;

    static auto s_hat_from_ys
      (
        const double &y1,
        const double &y2, 
        const double &kt2
      ) noexcept -> double;

    static auto t_hat_from_ys
      (
        const double &y1,
        const double &y2, 
        const double &kt2
      ) noexcept -> double;

    static auto u_hat_from_ys
      (
        const double &y1,
        const double &y2, 
        const double &kt2
      ) noexcept -> double;


    static auto calculate_sigma_jet
      (
        std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
        const double *const p_mand_s, 
        const double *const p_kt2_lower_cutoff, 
        pqcd::sigma_jet_params params
      ) noexcept -> double;

    static auto calculate_sigma_1jet_binned
      (
        std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
        const double *const p_mand_s, 
        const std::tuple<double, double, double, double> *const p_bin, 
        pqcd::sigma_jet_params params
      ) noexcept -> double;

    static auto calculate_sigma_jet_binned
      (
        std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
        const double *const p_mand_s, 
        const std::tuple<double, double, double, double> *const p_bin, 
        pqcd::sigma_jet_params params
      ) noexcept -> double;

    static auto calculate_sigma_dijet_binned
      (
        std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
        const double *const p_mand_s, 
        const std::tuple<double, double, double, double> *const p_bin, 
        pqcd::sigma_jet_params params
      ) noexcept -> double;

    static auto calculate_spatial_sigma_jet
      (
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
        /*std::shared_ptr<LHAPDF::GridPDF> p_n_pdf,*/ 
        const double *const p_mand_s, 
        const double *const p_kt2_lower_cutoff, 
        pqcd::sigma_jet_params params, 
        const double &sum_tppa, 
        const double &sum_tppb, 
        const double &tAA_0, 
        const double &tBB_0
      ) noexcept -> double;

protected:
private:
  std::mutex envelope_maximum_mutex;
};

#endif // PQCD_HPP
