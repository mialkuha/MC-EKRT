// Copyright (c) 2025 Mikko Kuha (University of Jyväskylä).
// This program is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details. You should have received a copy
// of the GNU General Public License along with this program. If not, see
// <http://www.gnu.org/licenses/>.

#include "pqcd.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846 /* pi */
#endif

auto pqcd::diff_sigma::full_partonic_bookkeeping(
    const std::array<double, 7> &f_i_x1,
    const std::array<double, 7> &f_i_x2,
    const std::array<double, 7> &f_ai_x1,
    const std::array<double, 7> &f_ai_x2,
    const uint_fast8_t &num_flavors,
    const double &s_hat,
    const double &t_hat,
    const double &u_hat) noexcept -> double
{
    double sum = 0;

    ///* GG->XX
    sum += 0.5 * f_i_x1[0] * f_i_x2[0] * pqcd::diff_sigma::sigma_gg_gg(s_hat, t_hat, u_hat);
    sum += num_flavors * f_i_x1[0] * f_i_x2[0] * pqcd::diff_sigma::sigma_gg_qaq(s_hat, t_hat, u_hat);

    //*/
    ///* GQ->XX
    for (uint_fast8_t flavor = 1; flavor <= num_flavors; ++flavor)
    {
        sum += f_i_x1[0] * f_i_x2[flavor] * pqcd::diff_sigma::sigma_gq_gq(s_hat, t_hat, u_hat);
        sum += f_i_x1[0] * f_ai_x2[flavor] * pqcd::diff_sigma::sigma_gq_gq(s_hat, t_hat, u_hat);
        sum += f_i_x1[flavor] * f_i_x2[0] * pqcd::diff_sigma::sigma_gq_gq(s_hat, t_hat, u_hat);
        sum += f_ai_x1[flavor] * f_i_x2[0] * pqcd::diff_sigma::sigma_gq_gq(s_hat, t_hat, u_hat);
    }
    //*/
    ///* QQ->XX
    for (uint_fast8_t flavor = 1; flavor <= num_flavors; ++flavor)
    {
        sum += 0.5 * f_i_x1[flavor] * f_i_x2[flavor] * pqcd::diff_sigma::sigma_qiqi_qiqi(s_hat, t_hat, u_hat);
        sum += 0.5 * f_ai_x1[flavor] * f_ai_x2[flavor] * pqcd::diff_sigma::sigma_qiqi_qiqi(s_hat, t_hat, u_hat);
    }

    for (uint_fast8_t flavor1 = 1; flavor1 <= num_flavors; ++flavor1)
    {
        for (uint_fast8_t flavor2 = 1; flavor2 <= num_flavors; ++flavor2)
        {
            if (flavor1 != flavor2)
            {
                sum += f_i_x1[flavor1] * f_i_x2[flavor2] * pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, t_hat, u_hat);
                sum += f_ai_x1[flavor1] * f_ai_x2[flavor2] * pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, t_hat, u_hat);
                sum += f_i_x1[flavor1] * f_ai_x2[flavor2] * pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, t_hat, u_hat);
                sum += f_ai_x1[flavor1] * f_i_x2[flavor2] * pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, t_hat, u_hat);
            }
        }
    }

    for (uint_fast8_t flavor = 1; flavor <= num_flavors; ++flavor)
    {
        sum += f_i_x1[flavor] * f_ai_x2[flavor] * (pqcd::diff_sigma::sigma_qiaqi_qiaqi(s_hat, t_hat, u_hat) + 0.5 * pqcd::diff_sigma::sigma_qiaqi_gg(s_hat, t_hat, u_hat) + (num_flavors - 1) * pqcd::diff_sigma::sigma_qiaqi_qjaqj(s_hat, t_hat, u_hat));
        sum += f_ai_x1[flavor] * f_i_x2[flavor] * (pqcd::diff_sigma::sigma_qiaqi_qiaqi(s_hat, t_hat, u_hat) + 0.5 * pqcd::diff_sigma::sigma_qiaqi_gg(s_hat, t_hat, u_hat) + (num_flavors - 1) * pqcd::diff_sigma::sigma_qiaqi_qjaqj(s_hat, t_hat, u_hat));
    }
    //*/

    return sum;
}

auto pqcd::diff_sigma::full_partonic_bookkeeping_1jet(
    const std::array<double, 7> &f_i_x1,
    const std::array<double, 7> &f_i_x2,
    const std::array<double, 7> &f_ai_x1,
    const std::array<double, 7> &f_ai_x2,
    const uint_fast8_t &num_flavours,
    const double &s_hat,
    const double &t_hat,
    const double &u_hat

    ) noexcept -> double
{
    double sum = 0;

    ///* GG->XX
    sum += f_i_x1[0] * f_i_x2[0] * pqcd::diff_sigma::sigma_gg_gg(s_hat, t_hat, u_hat);
    sum += num_flavours * f_i_x1[0] * f_i_x2[0] * (pqcd::diff_sigma::sigma_gg_qaq(s_hat, t_hat, u_hat) + pqcd::diff_sigma::sigma_gg_qaq(s_hat, u_hat, t_hat));
    //*/
    ///* GQ->XX
    for (uint_fast8_t flavor = 1; flavor <= num_flavours; ++flavor)
    {
        sum += f_i_x1[0] * f_i_x2[flavor] * (pqcd::diff_sigma::sigma_gq_gq(s_hat, t_hat, u_hat) + pqcd::diff_sigma::sigma_gq_gq(s_hat, u_hat, t_hat));
        sum += f_i_x1[0] * f_ai_x2[flavor] * (pqcd::diff_sigma::sigma_gq_gq(s_hat, t_hat, u_hat) + pqcd::diff_sigma::sigma_gq_gq(s_hat, u_hat, t_hat));
        sum += f_i_x1[flavor] * f_i_x2[0] * (pqcd::diff_sigma::sigma_gq_gq(s_hat, t_hat, u_hat) + pqcd::diff_sigma::sigma_gq_gq(s_hat, u_hat, t_hat));
        sum += f_ai_x1[flavor] * f_i_x2[0] * (pqcd::diff_sigma::sigma_gq_gq(s_hat, t_hat, u_hat) + pqcd::diff_sigma::sigma_gq_gq(s_hat, u_hat, t_hat));
    }
    //*/
    ///* QQ->XX
    for (uint_fast8_t flavor = 1; flavor <= num_flavours; ++flavor)
    {
        sum += f_i_x1[flavor] * f_i_x2[flavor] * pqcd::diff_sigma::sigma_qiqi_qiqi(s_hat, t_hat, u_hat);
        sum += f_ai_x1[flavor] * f_ai_x2[flavor] * pqcd::diff_sigma::sigma_qiqi_qiqi(s_hat, t_hat, u_hat);
    }

    for (uint_fast8_t flavor1 = 1; flavor1 <= num_flavours; ++flavor1)
    {
        for (uint_fast8_t flavor2 = 1; flavor2 <= num_flavours; ++flavor2)
        {
            if (flavor1 != flavor2)
            {
                sum += f_i_x1[flavor1] * f_i_x2[flavor2] * (pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, t_hat, u_hat) + pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, u_hat, t_hat));
                sum += f_ai_x1[flavor1] * f_ai_x2[flavor2] * (pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, t_hat, u_hat) + pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, u_hat, t_hat));
                sum += f_i_x1[flavor1] * f_ai_x2[flavor2] * (pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, t_hat, u_hat) + pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, u_hat, t_hat));
                sum += f_ai_x1[flavor1] * f_i_x2[flavor2] * (pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, t_hat, u_hat) + pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, u_hat, t_hat));
            }
        }
    }

    for (uint_fast8_t flavor = 1; flavor <= num_flavours; ++flavor)
    {
        sum += f_i_x1[flavor] * f_ai_x2[flavor] *
               ((
                    pqcd::diff_sigma::sigma_qiaqi_qiaqi(s_hat, t_hat, u_hat) + pqcd::diff_sigma::sigma_qiaqi_qiaqi(s_hat, u_hat, t_hat)) +
                pqcd::diff_sigma::sigma_qiaqi_gg(s_hat, t_hat, u_hat) + (num_flavours - 1) * (pqcd::diff_sigma::sigma_qiaqi_qjaqj(s_hat, t_hat, u_hat) + pqcd::diff_sigma::sigma_qiaqi_qjaqj(s_hat, u_hat, t_hat)));
        sum += f_ai_x1[flavor] * f_i_x2[flavor] *
               ((
                    pqcd::diff_sigma::sigma_qiaqi_qiaqi(s_hat, t_hat, u_hat) + pqcd::diff_sigma::sigma_qiaqi_qiaqi(s_hat, u_hat, t_hat)) +
                pqcd::diff_sigma::sigma_qiaqi_gg(s_hat, t_hat, u_hat) + (num_flavours - 1) * (pqcd::diff_sigma::sigma_qiaqi_qjaqj(s_hat, t_hat, u_hat) + pqcd::diff_sigma::sigma_qiaqi_qjaqj(s_hat, u_hat, t_hat)));
    }
    //*/
    return sum;
}

auto pqcd::generate_2_to_2_scatt(
    const double &sqrt_s,
    const double &kt_min,
    const double &kt_max,
    std::uniform_real_distribution<double> unirand,
    std::shared_ptr<std::mt19937> eng,
    std::shared_ptr<pdf_builder> pdf,
    pqcd::sigma_jet_params params,
    const double &power_law,
    envelope_func &env_func,
    pqcd::nn_coll_params &nn_params) noexcept -> dijet_specs
{
    dijet_specs event;

    double y1_min, y1_max, y2_min, y2_max, y1, y2;
    double kt;
    double rand[4];
    double ratio;
    bool is_below = false;

    y1_min = static_cast<double>(-log(sqrt_s / kt_min));
    y1_max = static_cast<double>(log(sqrt_s / kt_min));
    y2_min = static_cast<double>(-log(sqrt_s / kt_min));
    y2_max = static_cast<double>(log(sqrt_s / kt_min));

    while (!is_below)
    {
        for (auto &r : rand)
        {
            r = unirand(*eng);
        }

        // kT with importance sampling:
        // If kT follows the distribution f(x), then take an envelope g(x) s.t. f(x)<=g(x) everywhere.
        // g(x) has a primitive G(x) and that has an inverse Ginv(x).
        // Now take x = Ginv( G(kt_min) + R1*(G(kt_max) - G(kt_min)) ),
        // and accept kt=x if f(x)>R2*g(x), where R1 and R2 are random numbers
        kt = env_func.prim_inv(env_func.prim(kt_min) + rand[0] * (env_func.prim(kt_max) - env_func.prim(kt_min)));

        // ys from uniform distribution
        y1 = y1_min + rand[1] * (y1_max - y1_min);
        y2 = y2_min + rand[2] * (y2_max - y2_min);

        auto xsection = pqcd::diff_cross_section_2jet(sqrt_s, kt, y1, y2, pdf, params, nn_params, false, false);

        double total_xsection = 0;

        for (auto xsect : xsection)
        {
            total_xsection += xsect.sigma;
        }

        ratio = total_xsection / env_func.func(kt);

        if (ratio > 1)
        {
            std::cout << std::endl
                      << "Limiting function smaller than cross section!!" << std::endl;
            std::cout << "kT = " << kt << std::endl;
            std::cout << "y1 = " << y1 << std::endl;
            std::cout << "y2 = " << y2 << std::endl;
            std::cout << "total_xsection = " << total_xsection << std::endl;
            std::cout << "limit = " << env_func.func(kt) << std::endl
                      << std::endl;
            std::cout << "Maybe increase envelope_marginal." << std::endl;
            continue;
        }

        if (ratio > rand[3])
        {
            is_below = true;
            double rand_xsect = total_xsection * unirand(*eng);
            double sum = 0;

            for (auto xsect : xsection)
            {
                sum += xsect.sigma;
                if (sum > rand_xsect)
                {
                    event.init1 = xsect.init1;
                    event.init2 = xsect.init2;
                    event.final1 = xsect.final1;
                    event.final2 = xsect.final2;
                    event.pro_pdf = xsect.pro_pdf;
                    event.tar_pdf = xsect.tar_pdf;
                    break;
                }
            }
        }
    }

    event.kt = kt;
    event.y1 = y1;
    event.y2 = y2;

    if (event.init1 == 1 || event.init1 == 2)
    {
        double x1 = (kt / sqrt_s) * (exp(y1) + exp(y2));
        double rando = unirand(*eng);
        event.valence1 = pdf->check_valence(x1, kt * kt, event.init1, nn_params.projectile_neutron, rando);
    }

    if (event.init2 == 1 || event.init2 == 2)
    {
        double x2 = (kt / sqrt_s) * (exp(-y1) + exp(-y2));
        double rando = unirand(*eng);
        event.valence2 = pdf->check_valence(x2, kt * kt, event.init2, nn_params.target_neutron, rando);
    }

    return event;
}

auto pqcd::calculate_sigma_jet(
    std::shared_ptr<pdf_builder> p_pdf,
    const double *const p_mand_s,
    const double *const p_kt2_lower_cutoff,
    pqcd::sigma_jet_params params,
    pqcd::nn_coll_params nn_params,
    const bool &average) noexcept -> double
{
    double sigma_jet, error;
    const double upper_limits[3] = {1, 1, 1};
    const double lower_limits[3] = {0, 0, 0};
    const unsigned fdim = 1;
    std::tuple<
        std::shared_ptr<pdf_builder>,
        const double *const,
        const double *const,
        pqcd::sigma_jet_params,
        pqcd::nn_coll_params,
        const bool>
        fdata = {p_pdf, p_mand_s, p_kt2_lower_cutoff, params, nn_params, average};

    int not_success;

    not_success = hcubature(
        fdim,                      // Integrand dimension
        pqcd::sigma_jet_integrand, // Integrand function
        &fdata,                    // Pointer to additional arguments
        3,                         // Variable dimension
        lower_limits,              // Variables minimum
        upper_limits,              // Variables maximum
        0,                         // Max n:o of function evaluations
        0,                         // Required absolute error
        pqcd::g_error_tolerance,   // Required relative error
        ERROR_INDIVIDUAL,          // Enumerate of which norm is used on errors
        &sigma_jet,                // Pointer to output
        &error                     // Pointer to error output
    );

    if (not_success != 0)
    {
        std::cout << "Problem with integration" << std::endl;
        return -1;
    }

    return sigma_jet;
}

auto pqcd::calculate_sigma_1jet_binned(
    std::shared_ptr<pdf_builder> p_pdf,
    const double *const p_mand_s,
    const std::tuple<double, double, double, double> *const p_bin,
    pqcd::sigma_jet_params params,
    pqcd::nn_coll_params nn_coll) noexcept -> double
{
    double sigma_jet, error;
    const double upper_limits[3] = {1, 1, 1};
    const double lower_limits[3] = {0, 0, 0};
    const unsigned fdim = 1;
    std::tuple<
        std::shared_ptr<pdf_builder>,
        const double *const,
        const std::tuple<double, double, double, double> *const,
        pqcd::sigma_jet_params,
        pqcd::nn_coll_params>
        fdata = {p_pdf, p_mand_s, p_bin, params, nn_coll};

    int not_success;

    not_success = hcubature(
        fdim,                              // Integrand dimension
        pqcd::sigma_1jet_integrand_binned, // Integrand function
        &fdata,                            // Pointer to additional arguments
        3,                                 // Variable dimension
        lower_limits,                      // Variables minimum
        upper_limits,                      // Variables maximum
        0,                                 // Max n:o of function evaluations
        0,                                 // Required absolute error
        pqcd::g_error_tolerance,           // Required relative error
        ERROR_INDIVIDUAL,                  // Enumerate of which norm is used on errors
        &sigma_jet,                        // Pointer to output
        &error                             // Pointer to error output
    );

    if (not_success != 0)
    {
        std::cout << "Problem with integration" << std::endl;
        return -1;
    }

    return sigma_jet;
}

auto pqcd::calculate_sigma_jet_binned(
    std::shared_ptr<pdf_builder> p_pdf,
    const double *const p_mand_s,
    const std::tuple<double, double, double, double> *const p_bin,
    pqcd::sigma_jet_params params,
    pqcd::nn_coll_params nn_params) noexcept -> double
{
    double sigma_jet, error;
    const double upper_limits[3] = {1, 1, 1};
    const double lower_limits[3] = {0, 0, 0};
    const unsigned fdim = 1;
    std::tuple<
        std::shared_ptr<pdf_builder>,
        const double *const,
        const std::tuple<double, double, double, double> *const,
        pqcd::sigma_jet_params,
        pqcd::nn_coll_params>
        fdata = {p_pdf, p_mand_s, p_bin, params, nn_params};

    int not_success;

    not_success = hcubature(
        fdim,                             // Integrand dimension
        pqcd::sigma_jet_integrand_binned, // Integrand function
        &fdata,                           // Pointer to additional arguments
        3,                                // Variable dimension
        lower_limits,                     // Variables minimum
        upper_limits,                     // Variables maximum
        0,                                // Max n:o of function evaluations
        0,                                // Required absolute error
        pqcd::g_error_tolerance,          // Required relative error
        ERROR_INDIVIDUAL,                 // Enumerate of which norm is used on errors
        &sigma_jet,                       // Pointer to output
        &error                            // Pointer to error output
    );

    if (not_success != 0)
    {
        std::cout << "Problem with integration" << std::endl;
        return -1;
    }

    return sigma_jet;
}

auto pqcd::calculate_sigma_dijet_binned(
    std::shared_ptr<pdf_builder> p_pdf,
    const double *const p_mand_s,
    const std::tuple<double, double, double, double> *const p_bin,
    pqcd::sigma_jet_params params,
    pqcd::nn_coll_params nn_params) noexcept -> double
{
    double sigma_jet, error;
    const double upper_limits[3] = {1, 1, 1};
    const double lower_limits[3] = {0, 0, 0};
    const unsigned fdim = 1;
    std::tuple<
        std::shared_ptr<pdf_builder>,
        const double *const,
        const std::tuple<double, double, double, double> *const,
        pqcd::sigma_jet_params,
        pqcd::nn_coll_params>
        fdata = {p_pdf, p_mand_s, p_bin, params, nn_params};

    int not_success;

    not_success = hcubature(
        fdim,                               // Integrand dimension
        pqcd::sigma_dijet_integrand_binned, // Integrand function
        &fdata,                             // Pointer to additional arguments
        3,                                  // Variable dimension
        lower_limits,                       // Variables minimum
        upper_limits,                       // Variables maximum
        0,                                  // Max n:o of function evaluations
        0,                                  // Required absolute error
        pqcd::g_error_tolerance,            // Required relative error
        ERROR_INDIVIDUAL,                   // Enumerate of which norm is used on errors
        &sigma_jet,                         // Pointer to output
        &error                              // Pointer to error output
    );

    if (not_success != 0)
    {
        std::cout << "Problem with integration" << std::endl;
        return -1;
    }

    return sigma_jet;
}

auto pqcd::diff_sigma::sigma_qiqj_qiqj(
    const double &s,
    const double &t,
    const double &u) noexcept -> double
{
    return 4 * ((s * s + u * u) / (t * t)) / 9;
}

auto pqcd::diff_sigma::sigma_qiqi_qiqi(
    const double &s,
    const double &t,
    const double &u) noexcept -> double
{
    return 4 * ((s * s + u * u) / (t * t) + (s * s + t * t) / (u * u) - (2 * s * s) / (3 * t * u)) / 9;
}

auto pqcd::diff_sigma::sigma_qiaqi_qjaqj(
    const double &s,
    const double &t,
    const double &u) noexcept -> double
{
    return 4 * ((t * t + u * u) / (s * s)) / 9;
}

auto pqcd::diff_sigma::sigma_qiaqi_qiaqi(
    const double &s,
    const double &t,
    const double &u) noexcept -> double
{
    return 4 * ((s * s + u * u) / (t * t) + (t * t + u * u) / (s * s) - (2 * u * u) / (3 * s * t)) / 9;
}

auto pqcd::diff_sigma::sigma_qiaqi_gg(
    const double &s,
    const double &t,
    const double &u) noexcept -> double
{
    return 8 * (t * t + u * u) * (4 / (9 * t * u) - 1 / (s * s)) / 3;
}

auto pqcd::diff_sigma::sigma_gg_qaq(
    const double &s,
    const double &t,
    const double &u) noexcept -> double
{
    return 3 * (t * t + u * u) * (4 / (9 * t * u) - 1 / (s * s)) / 8;
}

auto pqcd::diff_sigma::sigma_gq_gq(
    const double &s,
    const double &t,
    const double &u) noexcept -> double
{
    return (s * s + u * u) * (1 / (t * t) - 4 / (9 * s * u));
}

auto pqcd::diff_sigma::sigma_gg_gg(
    const double &s,
    const double &t,
    const double &u) noexcept -> double
{
    return 4.5 * (3.0 - (u * t) / (s * s) - (u * s) / (t * t) - (s * t) / (u * u));
}

auto pqcd::diff_cross_section_2jet(
    const double &sqrt_s,
    const double &kt,
    const double &y1,
    const double &y2,
    std::shared_ptr<pdf_builder> pdf,
    pqcd::sigma_jet_params sigma_params,
    pqcd::nn_coll_params &nn_params,
    bool max,
    bool debug // true prints the calculated processes
    ) noexcept -> std::vector<xsection_id>
{
    std::vector<xsection_id> xsection;

    const particle_id num_flavors = pdf->num_flavors();
    double xsect;
    xsection_id process;

    // Mandelstam variables for the subprocesses
    const double kt2 = kt * kt;
    const auto s_hat = pqcd::s_hat_from_ys(y1, y2, kt2);
    const auto t_hat = pqcd::t_hat_from_ys(y1, y2, kt2);
    const auto u_hat = pqcd::u_hat_from_ys(y1, y2, kt2);

    // x_1 and x_2 values
    const double x1 = (kt / sqrt_s) * (exp(y1) + exp(y2));
    const double x2 = (kt / sqrt_s) * (exp(-y1) + exp(-y2));

    // Factorization / renormalization scale
    double q2;
    switch (sigma_params.scale_c)
    {
    case scaled_from_kt:
        q2 = pow(sigma_params.scalar, 2) * kt2;
        break;
    case constant:
        q2 = pow(sigma_params.scalar, 2);
        break;
    default:
        q2 = kt2;
        break;
    }

    // Units for dsigma/dktdy1dy2 as [mb/GeV]
    const double alpha_s = pdf->alphasQ2(q2);
    const double units = sigma_params.K_factor * 2.0 * kt * (M_PI * alpha_s * alpha_s / (s_hat * s_hat)) * 10 / pow(FMGEV, 2);

    // Returns zero when outside the physical boundaries
    if ((x1 < 0) || (x1 > 1) || (x2 > 1) || (x2 < 0))
    {
        process.sigma = 0;
        process.init1 = 0;
        process.init2 = 0;
        process.final1 = 0;
        process.final2 = 0;
        xsection.push_back(process);
        return xsection;
    }

    // Nuclear modifications and PDF manipulations
    auto [f_i_x1, f_i_x2, f_ai_x1, f_ai_x2] = pdf->make_pdfs(
        x1,
        x2,
        q2,
        nn_params.target_neutron,
        nn_params.projectile_neutron,
        nn_params.sum_tppa,
        nn_params.sum_tppb,
        false,
        max);

    if (debug)
    {
        std::cout << "x1=" << x1 << " x2=" << x2 << " q2=" << q2 << std::endl;
        std::cout << "f_i_x1 : " << std::endl;
        for (auto f : f_i_x1)
        {
            std::cout << f << std::endl;
        }
        std::cout << "f_i_x2 : " << std::endl;
        for (auto f : f_i_x2)
        {
            std::cout << f << std::endl;
        }
        std::cout << "f_ai_x1 : " << std::endl;
        for (auto f : f_ai_x1)
        {
            std::cout << f << std::endl;
        }
        std::cout << "f_ai_x2 : " << std::endl;
        for (auto f : f_ai_x2)
        {
            std::cout << f << std::endl;
        }
    }

    xsection.reserve(176); // 176 different processes if n:o quark flavors = 5

    ///* GG -> XX
    // gg -> gg
    {
        xsect = units * 0.5 * f_i_x1[0] * f_i_x2[0] * pqcd::diff_sigma::sigma_gg_gg(s_hat, t_hat, u_hat);
        process.sigma = xsect;
        process.init1 = 0;
        process.init2 = 0;
        process.final1 = 0;
        process.final2 = 0;
        process.pro_pdf = f_i_x1[0];
        process.tar_pdf = f_i_x2[0];
        xsection.push_back(process);
        if (debug)
        {
            std::cout << "gg -> gg: 0+0 -> 0+0" << " : " << xsect << std::endl;
        }
    }

    // gg -> qaq
    for (uint_fast8_t flavor = 1; flavor <= num_flavors; ++flavor)
    {
        xsect = units * f_i_x1[0] * f_i_x2[0] * pqcd::diff_sigma::sigma_gg_qaq(s_hat, t_hat, u_hat);
        process.sigma = xsect;
        process.init1 = 0;
        process.init2 = 0;
        process.final1 = static_cast<int_fast8_t>(flavor);
        process.final2 = -static_cast<int_fast8_t>(flavor);
        process.pro_pdf = f_i_x1[0];
        process.tar_pdf = f_i_x2[0];
        xsection.push_back(process);
        if (debug)
        {
            std::cout << "gg -> qaq: 0+0 -> " << static_cast<int>(flavor) << '+' << -static_cast<int>(flavor) << " : " << xsect << std::endl;
        }
    }
    //*/

    ///* GQ->XX
    // gq -> gq
    for (uint_fast8_t flavor = 1; flavor <= num_flavors; ++flavor)
    {
        xsect = units * f_i_x1[0] * f_i_x2[flavor] * pqcd::diff_sigma::sigma_gq_gq(s_hat, t_hat, u_hat);
        process.sigma = xsect;
        process.init1 = 0;
        process.init2 = static_cast<int_fast8_t>(flavor);
        process.final1 = 0;
        process.final2 = static_cast<int_fast8_t>(flavor);
        process.pro_pdf = f_i_x1[0];
        process.tar_pdf = f_i_x2[flavor];
        xsection.push_back(process);
        if (debug)
        {
            std::cout << "gq -> gq: 0+" << static_cast<int>(flavor) << " -> 0+" << static_cast<int>(flavor) << " : " << xsect << std::endl;
        }

        xsect = units * f_i_x1[0] * f_ai_x2[flavor] * pqcd::diff_sigma::sigma_gq_gq(s_hat, t_hat, u_hat);
        process.sigma = xsect;
        process.init1 = 0;
        process.init2 = -static_cast<int_fast8_t>(flavor);
        process.final1 = 0;
        process.final2 = -static_cast<int_fast8_t>(flavor);
        process.pro_pdf = f_i_x1[0];
        process.tar_pdf = f_ai_x2[flavor];
        xsection.push_back(process);
        if (debug)
        {
            std::cout << "gq -> gq: 0+" << -static_cast<int>(flavor) << " -> 0+" << -static_cast<int>(flavor) << " : " << xsect << std::endl;
        }

        xsect = units * f_i_x1[flavor] * f_i_x2[0] * pqcd::diff_sigma::sigma_gq_gq(s_hat, t_hat, u_hat);
        process.sigma = xsect;
        process.init1 = static_cast<int_fast8_t>(flavor);
        process.init2 = 0;
        process.final1 = static_cast<int_fast8_t>(flavor);
        process.final2 = 0;
        process.pro_pdf = f_i_x1[flavor];
        process.tar_pdf = f_i_x2[0];
        xsection.push_back(process);
        if (debug)
        {
            std::cout << "gq -> gq: " << static_cast<int>(flavor) << "+0 -> " << static_cast<int>(flavor) << "+0" << " : " << xsect << std::endl;
        }

        xsect = units * f_ai_x1[flavor] * f_i_x2[0] * pqcd::diff_sigma::sigma_gq_gq(s_hat, t_hat, u_hat);
        process.sigma = xsect;
        process.init1 = -static_cast<int_fast8_t>(flavor);
        process.init2 = 0;
        process.final1 = -static_cast<int_fast8_t>(flavor);
        process.final2 = 0;
        process.pro_pdf = f_ai_x1[flavor];
        process.tar_pdf = f_i_x2[0];
        xsection.push_back(process);
        if (debug)
        {
            std::cout << "gq -> gq: " << -static_cast<int>(flavor) << "+0 -> " << -static_cast<int>(flavor) << "+0" << " : " << xsect << std::endl;
        }
    }

    //*/
    ///* QQ->XX
    // qiqi -> qiqi
    for (uint_fast8_t flavor = 1; flavor <= num_flavors; ++flavor)
    {
        xsect = units * 0.5 * f_i_x1[flavor] * f_i_x2[flavor] * pqcd::diff_sigma::sigma_qiqi_qiqi(s_hat, t_hat, u_hat);
        process.sigma = xsect;
        process.init1 = static_cast<int_fast8_t>(flavor);
        process.init2 = static_cast<int_fast8_t>(flavor);
        process.final1 = static_cast<int_fast8_t>(flavor);
        process.final2 = static_cast<int_fast8_t>(flavor);
        process.pro_pdf = f_i_x1[flavor];
        process.tar_pdf = f_i_x2[flavor];
        xsection.push_back(process);
        if (debug)
        {
            std::cout << "qiqi -> qiqi: " << static_cast<int>(flavor) << '+' << static_cast<int>(flavor) << " -> " << static_cast<int>(flavor) << '+' << static_cast<int>(flavor) << " : " << xsect << std::endl;
        }

        xsect = units * 0.5 * f_ai_x1[flavor] * f_ai_x2[flavor] * pqcd::diff_sigma::sigma_qiqi_qiqi(s_hat, t_hat, u_hat);
        process.sigma = xsect;
        process.init1 = -static_cast<int_fast8_t>(flavor);
        process.init2 = -static_cast<int_fast8_t>(flavor);
        process.final1 = -static_cast<int_fast8_t>(flavor);
        process.final2 = -static_cast<int_fast8_t>(flavor);
        process.pro_pdf = f_ai_x1[flavor];
        process.tar_pdf = f_ai_x2[flavor];
        xsection.push_back(process);
        if (debug)
        {
            std::cout << "qiqi -> qiqi: " << -static_cast<int>(flavor) << '+' << -static_cast<int>(flavor) << " -> " << -static_cast<int>(flavor) << '+' << -static_cast<int>(flavor) << " : " << xsect << std::endl;
        }
    }

    // qiqj -> qiqj
    for (uint_fast8_t flavor1 = 1; flavor1 <= num_flavors; ++flavor1)
    {
        for (uint_fast8_t flavor2 = 1; flavor2 <= num_flavors; ++flavor2)
        {
            if (flavor1 != flavor2)
            {
                xsect = units * f_i_x1[flavor1] * f_i_x2[flavor2] * pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, t_hat, u_hat);
                process.sigma = xsect;
                process.init1 = static_cast<int_fast8_t>(flavor1);
                process.init2 = static_cast<int_fast8_t>(flavor2);
                process.final1 = static_cast<int_fast8_t>(flavor1);
                process.final2 = static_cast<int_fast8_t>(flavor2);
                process.pro_pdf = f_i_x1[flavor1];
                process.tar_pdf = f_i_x2[flavor2];
                xsection.push_back(process);
                if (debug)
                {
                    std::cout << "qiqj -> qiqj: " << static_cast<int_fast8_t>(flavor1) << '+' << static_cast<int_fast8_t>(flavor2) << " -> " << static_cast<int_fast8_t>(flavor1) << '+' << static_cast<int_fast8_t>(flavor2) << " : " << xsect << std::endl;
                }

                xsect = units * f_ai_x1[flavor1] * f_ai_x2[flavor2] * pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, t_hat, u_hat);
                process.sigma = xsect;
                process.init1 = -static_cast<int_fast8_t>(flavor1);
                process.init2 = -static_cast<int_fast8_t>(flavor2);
                process.final1 = -static_cast<int_fast8_t>(flavor1);
                process.final2 = -static_cast<int_fast8_t>(flavor2);
                process.pro_pdf = f_ai_x1[flavor1];
                process.tar_pdf = f_ai_x2[flavor2];
                xsection.push_back(process);
                if (debug)
                {
                    std::cout << "qiqj -> qiqj: " << -static_cast<int_fast8_t>(flavor1) << '+' << -static_cast<int_fast8_t>(flavor2) << " -> " << -static_cast<int_fast8_t>(flavor1) << '+' << -static_cast<int_fast8_t>(flavor2) << " : " << xsect << std::endl;
                }

                xsect = units * f_i_x1[flavor1] * f_ai_x2[flavor2] * pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, t_hat, u_hat);
                process.sigma = xsect;
                process.init1 = static_cast<int_fast8_t>(flavor1);
                process.init2 = -static_cast<int_fast8_t>(flavor2);
                process.final1 = static_cast<int_fast8_t>(flavor1);
                process.final2 = -static_cast<int_fast8_t>(flavor2);
                process.pro_pdf = f_i_x1[flavor1];
                process.tar_pdf = f_ai_x2[flavor2];
                xsection.push_back(process);
                if (debug)
                {
                    std::cout << "qiqj -> qiqj: " << static_cast<int_fast8_t>(flavor1) << '+' << -static_cast<int_fast8_t>(flavor2) << " -> " << static_cast<int_fast8_t>(flavor1) << '+' << -static_cast<int_fast8_t>(flavor2) << " : " << xsect << std::endl;
                }

                xsect = units * f_ai_x1[flavor1] * f_i_x2[flavor2] * pqcd::diff_sigma::sigma_qiqj_qiqj(s_hat, t_hat, u_hat);
                process.sigma = xsect;
                process.init1 = -static_cast<int_fast8_t>(flavor1);
                process.init2 = static_cast<int_fast8_t>(flavor2);
                process.final1 = -static_cast<int_fast8_t>(flavor1);
                process.final2 = static_cast<int_fast8_t>(flavor2);
                process.pro_pdf = f_ai_x1[flavor1];
                process.tar_pdf = f_i_x2[flavor2];
                xsection.push_back(process);
                if (debug)
                {
                    std::cout << "qiqj -> qiqj: " << -static_cast<int_fast8_t>(flavor1) << '+' << static_cast<int_fast8_t>(flavor2) << " -> " << -static_cast<int_fast8_t>(flavor1) << '+' << static_cast<int_fast8_t>(flavor2) << " : " << xsect << std::endl;
                }
            }
        }
    }

    // qiaqi -> qiaqi
    for (uint_fast8_t flavor = 1; flavor <= num_flavors; ++flavor)
    {
        xsect = units * f_i_x1[flavor] * f_ai_x2[flavor] * pqcd::diff_sigma::sigma_qiaqi_qiaqi(s_hat, t_hat, u_hat);
        process.sigma = xsect;
        process.init1 = static_cast<int_fast8_t>(flavor);
        process.init2 = -static_cast<int_fast8_t>(flavor);
        process.final1 = static_cast<int_fast8_t>(flavor);
        process.final2 = -static_cast<int_fast8_t>(flavor);
        process.pro_pdf = f_i_x1[flavor];
        process.tar_pdf = f_ai_x2[flavor];
        xsection.push_back(process);
        if (debug)
        {
            std::cout << "qiaqi -> qiaqi: " << static_cast<int>(flavor) << '+' << -static_cast<int>(flavor) << " -> " << static_cast<int>(flavor) << '+' << -static_cast<int>(flavor) << " : " << xsect << std::endl;
        }

        xsect = units * f_ai_x1[flavor] * f_i_x2[flavor] * pqcd::diff_sigma::sigma_qiaqi_qiaqi(s_hat, t_hat, u_hat);
        process.sigma = xsect;
        process.init1 = -static_cast<int_fast8_t>(flavor);
        process.init2 = static_cast<int_fast8_t>(flavor);
        process.final1 = -static_cast<int_fast8_t>(flavor);
        process.final2 = static_cast<int_fast8_t>(flavor);
        process.pro_pdf = f_ai_x1[flavor];
        process.tar_pdf = f_i_x2[flavor];
        xsection.push_back(process);
        if (debug)
        {
            std::cout << "qiaqi -> qiaqi: " << -static_cast<int>(flavor) << '+' << static_cast<int>(flavor) << " -> " << -static_cast<int>(flavor) << '+' << static_cast<int>(flavor) << " : " << xsect << std::endl;
        }
    }

    // qiaqi -> gg
    for (uint_fast8_t flavor = 1; flavor <= num_flavors; ++flavor)
    {
        xsect = units * 0.5 * f_i_x1[flavor] * f_ai_x2[flavor] * pqcd::diff_sigma::sigma_qiaqi_gg(s_hat, t_hat, u_hat);
        process.sigma = xsect;
        process.init1 = static_cast<int_fast8_t>(flavor);
        process.init2 = -static_cast<int_fast8_t>(flavor);
        process.final1 = 0;
        process.final2 = 0;
        process.pro_pdf = f_i_x1[flavor];
        process.tar_pdf = f_ai_x2[flavor];
        xsection.push_back(process);
        if (debug)
        {
            std::cout << "qiaqi -> gg: " << static_cast<int>(flavor) << '+' << -static_cast<int>(flavor) << " -> " << 0 << '+' << 0 << " : " << xsect << std::endl;
        }

        xsect = units * 0.5 * f_ai_x1[flavor] * f_i_x2[flavor] * pqcd::diff_sigma::sigma_qiaqi_gg(s_hat, t_hat, u_hat);
        process.sigma = xsect;
        process.init1 = -static_cast<int_fast8_t>(flavor);
        process.init2 = static_cast<int_fast8_t>(flavor);
        process.final1 = 0;
        process.final2 = 0;
        process.pro_pdf = f_ai_x1[flavor];
        process.tar_pdf = f_i_x2[flavor];
        xsection.push_back(process);
        if (debug)
        {
            std::cout << "qiaqi -> gg: " << -static_cast<int>(flavor) << '+' << static_cast<int>(flavor) << " -> " << 0 << '+' << 0 << " : " << xsect << std::endl;
        }
    }

    // qiaqi -> qjaqj
    for (uint_fast8_t flavor1 = 1; flavor1 <= num_flavors; ++flavor1)
    {
        for (int_fast8_t flavor2 = 1; flavor2 <= num_flavors; ++flavor2)
        {
            if (flavor1 != flavor2)
            {
                xsect = units * f_i_x1[flavor1] * f_ai_x2[flavor1] * pqcd::diff_sigma::sigma_qiaqi_qjaqj(s_hat, t_hat, u_hat);
                process.sigma = xsect;
                process.init1 = static_cast<int_fast8_t>(flavor1);
                process.init2 = -static_cast<int_fast8_t>(flavor1);
                process.final1 = static_cast<int_fast8_t>(flavor2);
                process.final2 = -static_cast<int_fast8_t>(flavor2);
                process.pro_pdf = f_i_x1[flavor1];
                process.tar_pdf = f_ai_x2[flavor1];
                xsection.push_back(process);
                if (debug)
                {
                    std::cout << "qiaqi -> qjaqj: " << static_cast<int_fast8_t>(flavor1) << '+' << -static_cast<int_fast8_t>(flavor1) << " -> " << static_cast<int_fast8_t>(flavor2) << '+' << -static_cast<int_fast8_t>(flavor2) << " : " << xsect << std::endl;
                }

                xsect = units * f_ai_x1[flavor1] * f_i_x2[flavor1] * pqcd::diff_sigma::sigma_qiaqi_qjaqj(s_hat, t_hat, u_hat);
                process.sigma = xsect;
                process.init1 = -static_cast<int_fast8_t>(flavor1);
                process.init2 = static_cast<int_fast8_t>(flavor1);
                process.final1 = -static_cast<int_fast8_t>(flavor2);
                process.final2 = static_cast<int_fast8_t>(flavor2);
                process.pro_pdf = f_ai_x1[flavor1];
                process.tar_pdf = f_i_x2[flavor1];
                xsection.push_back(process);
                if (debug)
                {
                    std::cout << "qiaqi -> qjaqj: " << -static_cast<int_fast8_t>(flavor1) << '+' << static_cast<int_fast8_t>(flavor1) << " -> " << -static_cast<int_fast8_t>(flavor2) << '+' << static_cast<int_fast8_t>(flavor2) << " : " << xsect << std::endl;
                }
            }
        }
    }
    //*/

    return xsection;
}

auto pqcd::sigma_jet(
    const double &x1,
    const double &x2,
    const double &q2,
    std::shared_ptr<pdf_builder> pdf,
    const double &s_hat,
    const double &t_hat,
    const double &u_hat,
    const double &K_factor,
    pqcd::nn_coll_params &nn_params,
    const bool &average) noexcept -> double
{
    const double alpha_s = pdf->alphasQ2(q2);
    const uint_fast8_t num_flavors = pdf->num_flavors();

    auto [f_i_x1, f_i_x2, f_ai_x1, f_ai_x2] = pdf->make_pdfs(
        x1,
        x2,
        q2,
        nn_params.target_neutron,
        nn_params.projectile_neutron,
        nn_params.sum_tppa,
        nn_params.sum_tppb,
        average,
        false);

    auto d_sigma = pqcd::diff_sigma::full_partonic_bookkeeping(
        f_i_x1,
        f_i_x2,
        f_ai_x1,
        f_ai_x2,
        num_flavors,
        s_hat,
        t_hat,
        u_hat);

    return K_factor * (M_PI * alpha_s * alpha_s / (s_hat * s_hat)) * d_sigma;
}

auto pqcd::sigma_1jet_integrand_binned(
    unsigned ndim,
    const double *p_x,
    void *p_fdata,
    unsigned fdim,
    double *p_fval) noexcept -> int
{
    (void)ndim;
    (void)fdim; // To silence "unused" warnings

    auto [p_pdf, p_mand_s, p_bin, params, nn_params] =
        *(static_cast<std::tuple<
              std::shared_ptr<pdf_builder>,
              const double *const,
              const std::tuple<double, double, double, double> *const,
              pqcd::sigma_jet_params,
              pqcd::nn_coll_params> *>(p_fdata));

    auto [kt_low, kt_upp, y_low, y_upp] = *p_bin;

    // const double kt2_upp = fmin(pow(kt_upp,2) , *p_mand_s/4);
    kt_upp = fmin(kt_upp, sqrt(*p_mand_s) / 2.0);
    // const double kt2_low = pow(kt_low,2);
    // const double kt2 = kt2_low + p_x[0] * (kt2_upp - kt2_low);
    const double kt = kt_low + p_x[0] * (kt_upp - kt_low);
    // const auto sqrt_s_per_kt = sqrt(*p_mand_s / kt2);
    const auto sqrt_s_per_kt = sqrt(*p_mand_s) / kt;
    const auto kt2 = kt * kt;

    y_low = fmax(y_low, -acosh(sqrt_s_per_kt / 2));
    y_upp = fmin(y_upp, acosh(sqrt_s_per_kt / 2));
    const double y1 = y_low + p_x[1] * (y_upp - y_low);

    const auto y2_upp = log(sqrt_s_per_kt - exp(y1));
    const auto y2_low = -log(sqrt_s_per_kt - exp(-y1));
    const double y2 = y2_low + p_x[2] * (y2_upp - y2_low);

    // const double jacobian = (kt2_upp - kt2_low) * (y_upp - y_low) * (y2_upp - y2_low);
    const double jacobian = 2.0 * kt * (kt_upp - kt_low) * (y_upp - y_low) * (y2_upp - y2_low);

    double fac_scale;
    switch (params.scale_c)
    {
    case scaled_from_kt:
        fac_scale = pow(params.scalar, 2) * kt2;
        break;
    case constant:
        fac_scale = pow(params.scalar, 2);
        break;
    default:
        fac_scale = kt2;
        break;
    }

    const auto x1 = (exp(y1) + exp(y2)) / sqrt_s_per_kt;
    const auto x2 = (exp(-y1) + exp(-y2)) / sqrt_s_per_kt;

    if (std::isnan(y1) || std::isnan(y2) || std::isnan(x1) || std::isnan(x2) || x1 > 1 || x2 > 1)
    {
        p_fval[0] = 0;
        return 0;
    }

    const auto s_hat = pqcd::s_hat_from_ys(y1, y2, kt2);
    const auto t_hat = pqcd::t_hat_from_ys(y1, y2, kt2);
    const auto u_hat = pqcd::u_hat_from_ys(y1, y2, kt2);

    p_fval[0] = pqcd::sigma_jet(x1, x2, fac_scale, p_pdf, s_hat, t_hat, u_hat, params.K_factor, nn_params, false) * jacobian * 10 / pow(FMGEV, 2); // UNITS: mb

    return 0; // success
}

auto pqcd::sigma_jet_integrand_binned(
    unsigned ndim,
    const double *p_x,
    void *p_fdata,
    unsigned fdim,
    double *p_fval) noexcept -> int
{
    (void)ndim;
    (void)fdim; // To silence "unused" warnings

    auto [p_pdf, p_mand_s, p_bin, params, nn_params] =
        *(static_cast<std::tuple<
              std::shared_ptr<pdf_builder>,
              const double *const,
              const std::tuple<double, double, double, double> *const,
              pqcd::sigma_jet_params,
              pqcd::nn_coll_params> *>(p_fdata));

    auto [kt_low, kt_upp, y_low, y_upp] = *p_bin;

    // const double kt2_upp = fmin(pow(kt_upp,2) , *p_mand_s/4);
    kt_upp = fmin(kt_upp, sqrt(*p_mand_s) / 2.0);
    // const double kt2_low = pow(kt_low,2);
    // const double kt2 = kt2_low + p_x[0] * (kt2_upp - kt2_low);
    const double kt = kt_low + p_x[0] * (kt_upp - kt_low);
    // const auto sqrt_s_per_kt = sqrt(*p_mand_s / kt2);
    const auto sqrt_s_per_kt = sqrt(*p_mand_s) / kt;
    const auto kt2 = kt * kt;

    y_low = fmax(y_low, -acosh(sqrt_s_per_kt / 2));
    y_upp = fmin(y_upp, acosh(sqrt_s_per_kt / 2));
    const double y1 = y_low + p_x[1] * (y_upp - y_low);

    const auto y2_upp = log(sqrt_s_per_kt - exp(y1));
    const auto y2_low = -log(sqrt_s_per_kt - exp(-y1));
    const double y2 = y2_low + p_x[2] * (y2_upp - y2_low);

    // const double jacobian = (kt2_upp - kt2_low) * (y_upp - y_low) * (y2_upp - y2_low);
    const double jacobian = 2.0 * kt * (kt_upp - kt_low) * (y_upp - y_low) * (y2_upp - y2_low);

    double fac_scale;
    switch (params.scale_c)
    {
    case scaled_from_kt:
        fac_scale = pow(params.scalar, 2) * kt2;
        break;
    case constant:
        fac_scale = pow(params.scalar, 2);
        break;
    default:
        fac_scale = kt2;
        break;
    }

    const auto x1 = (exp(y1) + exp(y2)) / sqrt_s_per_kt;
    const auto x2 = (exp(-y1) + exp(-y2)) / sqrt_s_per_kt;

    if (std::isnan(y1) || std::isnan(y2) || std::isnan(x1) || std::isnan(x2) || x1 > 1 || x2 > 1)
    {
        p_fval[0] = 0;
        return 0;
    }

    const auto s_hat = pqcd::s_hat_from_ys(y1, y2, kt2);
    const auto t_hat = pqcd::t_hat_from_ys(y1, y2, kt2);
    const auto u_hat = pqcd::u_hat_from_ys(y1, y2, kt2);

    p_fval[0] = pqcd::sigma_jet(x1, x2, fac_scale, p_pdf, s_hat, t_hat, u_hat, params.K_factor, nn_params, false) * jacobian * 10 / pow(FMGEV, 2); // UNITS: mb

    return 0; // success
}

auto pqcd::sigma_dijet_integrand_binned(
    unsigned ndim,
    const double *p_x,
    void *p_fdata,
    unsigned fdim,
    double *p_fval) noexcept -> int
{
    (void)ndim;
    (void)fdim; // To silence "unused" warnings

    auto [p_pdf, p_mand_s, p_bin, params, nn_params] =
        *(static_cast<std::tuple<
              std::shared_ptr<pdf_builder>,
              const double *const,
              const std::tuple<double, double, double, double> *const,
              pqcd::sigma_jet_params,
              pqcd::nn_coll_params> *>(p_fdata));

    auto [kt_low, kt_upp, eta_low, eta_upp] = *p_bin;

    const double kt2_upp = fmin(pow(kt_upp, 2), *p_mand_s / 4);
    const double kt2_low = pow(kt_low, 2);
    const double kt2 = kt2_low + p_x[0] * (kt2_upp - kt2_low);
    const auto sqrt_s_per_kt = sqrt(*p_mand_s / kt2);

    const double eta = eta_low + p_x[1] * (eta_upp - eta_low);

    const auto ystar_upp = acosh(sqrt_s_per_kt * exp(-abs(eta)) / 2.0);
    const auto ystar_low = 0.0;
    const double ystar = ystar_low + p_x[2] * (ystar_upp - ystar_low);

    const double jacobian = (kt2_upp - kt2_low) * (eta_upp - eta_low) * (ystar_upp - ystar_low);

    double fac_scale;
    switch (params.scale_c)
    {
    case scaled_from_kt:
        fac_scale = pow(params.scalar, 2) * kt2;
        break;
    case constant:
        fac_scale = pow(params.scalar, 2);
        break;
    default:
        fac_scale = kt2;
        break;
    }

    const auto x1 = 2.0 * exp(eta) * cosh(ystar) / sqrt_s_per_kt;
    const auto x2 = 2.0 * exp(-eta) * cosh(ystar) / sqrt_s_per_kt;

    if (std::isnan(eta) || std::isnan(ystar) || std::isnan(x1) || std::isnan(x2) || x1 > 1 || x2 > 1)
    {
        p_fval[0] = 0;
        return 0;
    }

    const auto s_hat = 2.0 * kt2 * (1 + cosh(2.0 * ystar));
    const auto t_hat = -kt2 * (1 + exp(-2.0 * ystar));
    const auto u_hat = -kt2 * (1 + exp(2.0 * ystar));

    p_fval[0] = pqcd::sigma_jet(x1, x2, fac_scale, p_pdf, s_hat, t_hat, u_hat, params.K_factor, nn_params, false) * 2.0 * jacobian * 10 / pow(FMGEV, 2); // UNITS: mb

    return 0; // success
}

auto pqcd::sigma_jet_integrand(
    unsigned ndim,
    const double *p_x,
    void *p_fdata,
    unsigned fdim,
    double *p_fval) noexcept -> int
{
    (void)ndim;
    (void)fdim; // To silence "unused" warnings
    auto [p_pdf, p_mand_s, p_p02, params, nn_params, average] = *(static_cast<std::tuple<std::shared_ptr<pdf_builder>,
                                                                                         const double *const,
                                                                                         const double *const,
                                                                                         pqcd::sigma_jet_params,
                                                                                         pqcd::nn_coll_params,
                                                                                         const bool> *>(p_fdata));

    double kt2;
    double y1, y2;
    double jacobian;
    pqcd::scale_limits_from_0_1(p_x[0], p_x[1], p_x[2], *p_p02, *p_mand_s, kt2, y1, y2, jacobian);
    const auto sqrt_s_per_kt = sqrt(*p_mand_s / kt2);

    double fac_scale;

    switch (params.scale_c)
    {
    case scaled_from_kt:
        fac_scale = pow(params.scalar, 2) * kt2;
        break;
    case constant:
        fac_scale = pow(params.scalar, 2);
        break;
    default:
        fac_scale = kt2;
        break;
    }

    const auto x1 = (exp(y1) + exp(y2)) / sqrt_s_per_kt;
    const auto x2 = (exp(-y1) + exp(-y2)) / sqrt_s_per_kt;

    if (std::isnan(y1) || std::isnan(y2) || std::isnan(x1) || std::isnan(x2) || x1 > 1 || x2 > 1)
    {
        p_fval[0] = 0;
        return 0;
    }

    const auto s_hat = pqcd::s_hat_from_ys(y1, y2, kt2);
    const auto t_hat = pqcd::t_hat_from_ys(y1, y2, kt2);
    const auto u_hat = pqcd::u_hat_from_ys(y1, y2, kt2);

    // SES
    if (params.use_ses)
    {
        // TODO
        // const auto alpha_s = p_pdf->alphasQ2(fac_scale);
        // const auto subprocess_cs = (M_PI * alpha_s * alpha_s / (s_hat * s_hat)) * pqcd::diff_sigma::sigma_gg_gg(s_hat, t_hat, u_hat);
        // p_fval[0] = 0.5 * pqcd::f_ses(x1, fac_scale, p_pdf) * pqcd::f_ses(x2, fac_scale, p_pdf) * subprocess_cs * jacobian * 10 / pow(FMGEV, 2); //UNITS: mb
    }
    else // FULL SUMMATION
    {
        p_fval[0] = pqcd::sigma_jet(x1, x2, fac_scale, p_pdf, s_hat, t_hat, u_hat, params.K_factor, nn_params, average) * jacobian * 10 / pow(FMGEV, 2); // UNITS: mb
    }

    return 0; // success
}

// auto pqcd::f_ses
// (
//     const double &x,
//     const double &q2,
//     std::shared_ptr<pdf_builder> p_pdf
// ) noexcept -> double
// {
//     double sum = 0;
//     for (int_fast8_t flavor = -5; flavor <= 5; ++flavor)
//     {
//         sum += (flavor == 0) ? p_pdf->xfxQ2(flavor, x, q2) : 4 * p_pdf->xfxQ2(flavor, x, q2) / 9;
//     }
//     //cout<<"x="<<*p_x<<", q2="<<*p_q2<<" SES="<<sum<<endl;
//     return sum;
// }

auto pqcd::scale_limits_from_0_1(
    const double &z1,
    const double &z2,
    const double &z3, // variables between 0 and 1
    const double &kt2_lower_cutoff,
    const double &mand_s, // parameters
    double &kt2,
    double &y1,
    double &y2,
    double &jacobian) noexcept -> void
{
    // kt2 = kt2_lower_cutoff + z1 * ((mand_s / 4) - kt2_lower_cutoff);
    // const auto sqrt_s_per_kt = sqrt(mand_s / kt2);
    const auto kt_upp = sqrt(mand_s) / 2.0;
    const auto kt_low = sqrt(kt2_lower_cutoff);
    const auto kt = kt_low + z1 * (kt_upp - kt_low);
    const auto sqrt_s_per_kt = sqrt(mand_s) / kt;
    kt2 = kt * kt;

    const auto y1_upper = acosh(sqrt_s_per_kt / 2);
    // const auto y1_lower = -y1_upper;
    y1 = (-1 + (2 * z2)) * y1_upper; // y1_lower + z2 * (y1_upper - y1_lower)

    const auto y2_upper = log(sqrt_s_per_kt - exp(y1));
    const auto y2_lower = -log(sqrt_s_per_kt - exp(-y1));
    y2 = y2_lower + z3 * (y2_upper - y2_lower);

    // jacobian = ((mand_s / 4) - kt2_lower_cutoff) * (2 * y1_upper) * (y2_upper - y2_lower);
    jacobian = 2.0 * kt * (kt_upp - kt_low) * (2 * y1_upper) * (y2_upper - y2_lower);
}

auto pqcd::s_hat_from_ys(
    const double &y1,
    const double &y2,
    const double &kt2) noexcept -> double
{
    return 2 * kt2 * (1 + cosh(y1 - y2));
}

auto pqcd::t_hat_from_ys(
    const double &y1,
    const double &y2,
    const double &kt2) noexcept -> double
{
    return -kt2 * (1 + exp(-y1 + y2));
}

auto pqcd::u_hat_from_ys(
    const double &y1,
    const double &y2,
    const double &kt2) noexcept -> double
{
    return -kt2 * (1 + exp(y1 - y2));
}
