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

            explicit params(bool projectile_with_npdfs_ = false, bool target_with_npdfs_ = false, bool isoscalar_projectile_ = false, bool isoscalar_target_ = false, int npdf_setnumber_ = 1) noexcept
                : projectile_with_npdfs(projectile_with_npdfs_), target_with_npdfs(target_with_npdfs_), isoscalar_projectile(isoscalar_projectile_), isoscalar_target(isoscalar_target_), npdf_setnumber(npdf_setnumber_) {}
        };
        ///sigma_qiqj_qiqj
        ///Calculates and returns the cross section of the subprocess qiqj->qiqj. Return value should be real.
        ///
        ///\param p_s_hat = pointer to subprocess mandelstam variable s
        ///\param p_t_hat = pointer to subprocess mandelstam variable t
        ///\param p_u_hat = pointer to subprocess mandelstam variable u
        ///\param p_alpha_s = pointer to strong interaction constant
        ///
        static xsectval sigma_qiqj_qiqj(const momentum * p_s_hat, const momentum * p_t_hat, const momentum * p_u_hat, const momentum alpha_s) noexcept;

        ///sigma_qiqi_qiqi
        ///Calculates and returns the cross section of the subprocess qiqi->qiqi. Return value should be real.
        ///
        ///\param p_s_hat = pointer to subprocess mandelstam variable s
        ///\param p_t_hat = pointer to subprocess mandelstam variable t
        ///\param p_u_hat = pointer to subprocess mandelstam variable u
        ///\param p_alpha_s = pointer to strong interaction constant
        ///
        static xsectval sigma_qiqi_qiqi(const momentum * p_s_hat, const momentum * p_t_hat, const momentum * p_u_hat, const momentum alpha_s) noexcept;

        ///sigma_qiaqi_qjaqj
        ///Calculates and returns the cross section of the subprocess qiaqi->qjaqj. Return value should be real.
        ///
        ///\param p_s_hat = pointer to subprocess mandelstam variable s
        ///\param p_t_hat = pointer to subprocess mandelstam variable t
        ///\param p_u_hat = pointer to subprocess mandelstam variable u
        ///\param p_alpha_s = pointer to strong interaction constant
        ///
        static xsectval sigma_qiaqi_qjaqj( const momentum * p_s_hat,  const momentum * p_t_hat,  const momentum * p_u_hat,  const momentum alpha_s) noexcept;

        ///sigma_qiaqi_qiaqi
        ///Calculates and returns the cross section of the subprocess qiaqi->qiaqi. Return value should be real.
        ///
        ///\param p_s_hat = pointer to subprocess mandelstam variable s
        ///\param p_t_hat = pointer to subprocess mandelstam variable t
        ///\param p_u_hat = pointer to subprocess mandelstam variable u
        ///\param p_alpha_s = pointer to strong interaction constant
        ///
        static xsectval sigma_qiaqi_qiaqi( const momentum * p_s_hat,  const momentum * p_t_hat,  const momentum * p_u_hat,  const momentum alpha_s) noexcept;

        ///sigma_qiaqi_gg
        ///Calculates and returns the cross section of the subprocess qiaqi->gg. Return value should be real.
        ///
        ///\param p_s_hat = pointer to subprocess mandelstam variable s
        ///\param p_t_hat = pointer to subprocess mandelstam variable t
        ///\param p_u_hat = pointer to subprocess mandelstam variable u
        ///\param p_alpha_s = pointer to strong interaction constant
        ///
        static xsectval sigma_qiaqi_gg( const momentum * p_s_hat,  const momentum * p_t_hat,  const momentum * p_u_hat,  const momentum alpha_s) noexcept;

        ///sigma_gg_qaq
        ///Calculates and returns the cross section of the subprocess gg->qaq. Return value should be real.
        ///
        ///\param p_s_hat = pointer to subprocess mandelstam variable s
        ///\param p_t_hat = pointer to subprocess mandelstam variable t
        ///\param p_u_hat = pointer to subprocess mandelstam variable u
        ///\param p_alpha_s = pointer to strong interaction constant
        ///
        static xsectval sigma_gg_qaq( const momentum * p_s_hat,  const momentum * p_t_hat,  const momentum * p_u_hat,  const momentum alpha_s) noexcept;

        ///sigma_gq_gq
        ///Calculates and returns the cross section of the subprocess gg->gg. Return value should be real.
        ///
        ///\param p_s_hat = pointer to subprocess mandelstam variable s
        ///\param p_t_hat = pointer to subprocess mandelstam variable t
        ///\param p_u_hat = pointer to subprocess mandelstam variable u
        ///\param p_alpha_s = pointer to strong interaction constant
        ///
        static xsectval sigma_gq_gq( const momentum * p_s_hat,  const momentum * p_t_hat,  const momentum * p_u_hat,  const momentum alpha_s) noexcept;

        ///sigma_gg_gg
        ///Calculates and returns the cross section of the subprocess gg->gg. Return value should be real.
        ///
        ///\param p_s_hat = pointer to subprocess mandelstam variable s
        ///\param p_t_hat = pointer to subprocess mandelstam variable t
        ///\param p_u_hat = pointer to subprocess mandelstam variable u
        ///\param p_alpha_s = pointer to strong interaction constant
        ///
        static xsectval sigma_gg_gg( const momentum * p_s_hat,  const momentum * p_t_hat,  const momentum * p_u_hat,  const momentum alpha_s) noexcept;

        ///diff_sigma_jet
        ///Calculates and returns the differential cross section.
        ///
        ///\param p_x1 = pointer to const momentum fraction of the first particle in question
        ///\param p_x2 = pointer to const momentum fraction of the first particle in question
        ///\param p_q2 = pointer to const momentum exchange squared
        ///\param p_pdf = pointer to a member of the PDF class from LHAPDF
        ///\param p_s_hat = pointer to subprocess mandelstam variable s
        ///\param p_t_hat = pointer to subprocess mandelstam variable t
        ///\param p_u_hat = pointer to subprocess mandelstam variable u
        ///
        static xsectval sigma_jet(const rapidity * p_x1,    const rapidity * p_x2, 
                                  const momentum * p_q2,    std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
                                  const momentum * p_s_hat, const momentum * p_t_hat, 
                                  const momentum * p_u_hat, const params * p_params) noexcept;
        static xsectval spatial_sigma_jet(const rapidity *const p_x1,    const rapidity *const p_x2, 
                                          const momentum *const p_q2,    std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
                                          const momentum *const p_s_hat, const momentum *const p_t_hat, 
                                          const momentum *const p_u_hat, const pqcd::diff_sigma::params *const p_params, 
                                          std::shared_ptr<LHAPDF::GridPDF> p_n_pdf, 
                                          std::function<double(double const&)> rA_spatial,
                                          std::function<double(double const&)> rB_spatial) noexcept;
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

        explicit sigma_jet_params(scale_choice scale_c_ = scaled_from_kt, momentum scalar_ = 1.0, bool projectile_with_npdfs_ = false, bool target_with_npdfs_ = false, bool isoscalar_projectile_ = false, bool isoscalar_target_ = false, int npdf_setnumber_ = 1, bool use_ses_ = false) noexcept
            : scale_c(scale_c_), scalar(scalar_), p_d_params(new diff_sigma::params(projectile_with_npdfs_, target_with_npdfs_, isoscalar_projectile_, isoscalar_target_, npdf_setnumber_)), use_ses(use_ses_)
        { /*p_d_params = new diff_sigma::params(projectile_with_npdfs_, target_with_npdfs_, isoscalar_projectile_, isoscalar_target_, npdf_setnumber_);*/
        }
        explicit sigma_jet_params(diff_sigma::params *p_d_params_, scale_choice scale_c_ = scaled_from_kt, momentum scalar_ = 1.0, bool use_ses_ = false) noexcept
            : scale_c(scale_c_), scalar(scalar_), p_d_params(p_d_params_), use_ses(use_ses_) {}
    };
    static void generate_bin_NN_coll(nn_coll * coll) noexcept;
    static int sigma_jet_integrand(unsigned ndim, const double *p_x, void *p_fdata, unsigned fdim, double *p_fval) noexcept;
    static void scale_limits_from_0_1(const rapidity & z1, const rapidity & z2, const rapidity & z3,                      //variables between 0 and 1
                                      const momentum *const kt2_lower_cutoff, const momentum *const mand_s,               //parameters
                                      momentum * p_kt2, rapidity * p_y1, rapidity * p_y2, xsectval * p_jacobian) noexcept;//output
    static int spatial_sigma_jet_integrand(unsigned ndim, const double *p_x, void *p_fdata, unsigned fdim, double *p_fval) noexcept;
    static xsectval f_ses(const rapidity * p_x, const momentum * p_q2, std::shared_ptr<LHAPDF::GridPDF> p_pdf) noexcept;
    static momentum s_hat_from_ys(const rapidity * p_y1, const rapidity * p_y2, const momentum * p_kt2) noexcept;
    static momentum t_hat_from_ys(const rapidity * p_y1, const rapidity * p_y2, const momentum * p_kt2) noexcept;
    static momentum u_hat_from_ys(const rapidity * p_y1, const rapidity * p_y2, const momentum * p_kt2) noexcept;
    static xsectval calculate_sigma_jet(std::shared_ptr<LHAPDF::GridPDF> p_pdf, const momentum *const p_mand_s, const momentum *const p_kt2_lower_cutoff, const pqcd::sigma_jet_params *const p_params) noexcept;
    static xsectval calculate_spatial_sigma_jet_mf(std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, std::shared_ptr<LHAPDF::GridPDF> p_n_pdf, const momentum *const p_mand_s, 
                                                    const momentum *const p_kt2_lower_cutoff, const pqcd::sigma_jet_params *const p_params, const spatial *const p_sum_tppa, 
                                                    const spatial *const p_sum_tppb, const spatial *const p_tAA_0, const spatial *const p_tBB_0) noexcept;
    static xsectval calculate_spatial_sigma_jet_full(std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, std::shared_ptr<LHAPDF::GridPDF> p_n_pdf, const momentum *const p_mand_s, 
                                                    const momentum *const p_kt2_lower_cutoff, const pqcd::sigma_jet_params *const p_params, const spatial *const p_sum_tppa, 
                                                    const spatial *const p_sum_tppb, const spatial *const p_tAA_0, const spatial *const p_tBB_0) noexcept;
protected:
private:
};

#endif // PQCD_HPP
