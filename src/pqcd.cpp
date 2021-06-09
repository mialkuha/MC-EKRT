//Copyright (c) 2021 Mikko Kuha

#include "pqcd.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846 /* pi */
#endif
#include "eps09.cxx"

void pqcd::generate_bin_NN_coll(nn_coll * coll) noexcept
{
    //TODO
}

xsectval pqcd::calculate_sigma_jet(std::shared_ptr<LHAPDF::GridPDF> p_pdf, const momentum *const p_mand_s, const momentum *const p_kt2_lower_cutoff, const pqcd::sigma_jet_params *const p_params) noexcept
{
    xsectval sigma_jet, error;
    const double upper_limits[3] = {1, 1, 1};
    const double lower_limits[3] = {0, 0, 0};
    const unsigned fdim = 1;
    std::tuple<std::shared_ptr<LHAPDF::GridPDF>, const momentum *const, const momentum *const, const pqcd::sigma_jet_params *const> fdata =
        {p_pdf, p_mand_s, p_kt2_lower_cutoff, p_params};

    int not_success;

    not_success = hcubature(fdim,                               //Integrand dimension
                            pqcd::sigma_jet_integrand,          //Integrand function
                            &fdata,                             //Pointer to additional arguments
                            3,                                  //Variable dimension
                            lower_limits,                       //Variables minimum
                            upper_limits,                       //Variables maximum
                            0,                                  //Max n:o of function evaluations
                            0,                                  //Required absolute error
                            pqcd::g_error_tolerance,    //Required relative error
                            ERROR_INDIVIDUAL,                   //Enumerate of which norm is used on errors
                            &sigma_jet,                            //Pointer to output
                            &error);                           //Pointer to error output

    if (not_success != 0)
    {
        std::cout << "Problem with integration" << std::endl;
        return -1;
    }

    return sigma_jet;
}

xsectval pqcd::diff_sigma::sigma_qiqj_qiqj(const momentum *const p_s_hat, const momentum *const p_t_hat, const momentum *const p_u_hat, const momentum alpha_s) noexcept
{
    const momentum s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s)) * 4 * ((s * s + u * u) / (t * t)) / 9;
}

xsectval pqcd::diff_sigma::sigma_qiqi_qiqi(const momentum *const p_s_hat, const momentum *const p_t_hat, const momentum *const p_u_hat, const momentum alpha_s) noexcept
{
    const momentum s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s)) * 4 * ((s * s + u * u) / (t * t) + (s * s + t * t) / (u * u) - (2 * s * s) / (3 * t * u)) / 9;
}

xsectval pqcd::diff_sigma::sigma_qiaqi_qjaqj(const momentum *const p_s_hat, const momentum *const p_t_hat, const momentum *const p_u_hat, const momentum alpha_s) noexcept
{
    const momentum s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s)) * 4 * ((t * t + u * u) / (s * s)) / 9;
}

xsectval pqcd::diff_sigma::sigma_qiaqi_qiaqi(const momentum *const p_s_hat, const momentum *const p_t_hat, const momentum *const p_u_hat, const momentum alpha_s) noexcept
{
    const momentum s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s)) * 4 * ((s * s + u * u) / (t * t) + (t * t + u * u) / (s * s) - (2 * u * u) / (3 * s * t)) / 9;
}

xsectval pqcd::diff_sigma::sigma_qiaqi_gg(const momentum *const p_s_hat, const momentum *const p_t_hat, const momentum *const p_u_hat, const momentum alpha_s) noexcept
{
    const momentum s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s)) * 8 * (t * t + u * u) * (4 / (9 * t * u) - 1 / (s * s)) / 3;
}

xsectval pqcd::diff_sigma::sigma_gg_qaq(const momentum *const p_s_hat, const momentum *const p_t_hat, const momentum *const p_u_hat, const momentum alpha_s) noexcept
{
    const momentum s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s)) * 3 * (t * t + u * u) * (4 / (9 * t * u) - 1 / (s * s)) / 8;
}

xsectval pqcd::diff_sigma::sigma_gq_gq(const momentum *const p_s_hat, const momentum *const p_t_hat, const momentum *const p_u_hat, const momentum alpha_s) noexcept
{
    const momentum s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s)) * (s * s + u * u) * (1 / (t * t) - 4 / (9 * s * u));
}

xsectval pqcd::diff_sigma::sigma_gg_gg(const momentum *const p_s_hat, const momentum *const p_t_hat, const momentum *const p_u_hat, const momentum alpha_s) noexcept
{
    const momentum s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s)) * 4.5 * (3.0 - (u * t) / (s * s) - (u * s) / (t * t) - (s * t) / (u * u));
}

xsectval pqcd::diff_sigma::sigma_jet(const rapidity *const p_x1, const rapidity *const p_x2, const momentum *const p_q2, std::shared_ptr<LHAPDF::GridPDF> p_pdf, const momentum *const p_s_hat, const momentum *const p_t_hat, const momentum *const p_u_hat, const pqcd::diff_sigma::params *const p_params) noexcept
{
    const rapidity x1 = *p_x1, x2 = *p_x2;
    const momentum q2 = *p_q2;
    xsectval sum = 0;
    const double alphas = p_pdf->alphasQ2(q2);
    const int numFlavours = std::stoi(p_pdf->info().get_entry("NumFlavors"));

    double ruv = 1.0, rdv = 1.0, ru = 1.0, rd = 1.0, rs = 1.0, rc = 1.0, rb = 1.0, rg = 1.0;

    if (p_params->projectile_with_npdfs)
    {
        eps09(1, p_params->npdf_setnumber, 208, x1, sqrt(q2), ruv, rdv, ru, rd, rs, rc, rb, rg);
    }
    double f_i_x1[] = {rg * p_pdf->xfxQ2(0, x1, q2),
                       ruv * p_pdf->xfxQ2(1, x1, q2) + (ru - ruv) * p_pdf->xfxQ2(-1, x1, q2),
                       rdv * p_pdf->xfxQ2(2, x1, q2) + (rd - rdv) * p_pdf->xfxQ2(-2, x1, q2),
                       rs * p_pdf->xfxQ2(3, x1, q2),
                       rc * p_pdf->xfxQ2(4, x1, q2),
                       rb * p_pdf->xfxQ2(5, x1, q2),
                       rb * p_pdf->xfxQ2(6, x1, q2)};
    double f_ai_x1[] = {rg * p_pdf->xfxQ2(0, x1, q2),
                        ru * p_pdf->xfxQ2(-1, x1, q2),
                        rd * p_pdf->xfxQ2(-2, x1, q2),
                        rs * p_pdf->xfxQ2(-3, x1, q2),
                        rc * p_pdf->xfxQ2(-4, x1, q2),
                        rb * p_pdf->xfxQ2(-5, x1, q2),
                        rb * p_pdf->xfxQ2(-6, x1, q2)};

    if (p_params->target_with_npdfs)
    {
        eps09(1, p_params->npdf_setnumber, 208, x2, sqrt(q2), ruv, rdv, ru, rd, rs, rc, rb, rg);
    }
    double f_i_x2[] = {rg * p_pdf->xfxQ2(0, x2, q2),
                       ruv * p_pdf->xfxQ2(1, x2, q2) + (ru - ruv) * p_pdf->xfxQ2(-1, x2, q2),
                       rdv * p_pdf->xfxQ2(2, x2, q2) + (rd - rdv) * p_pdf->xfxQ2(-2, x2, q2),
                       rs * p_pdf->xfxQ2(3, x2, q2),
                       rc * p_pdf->xfxQ2(4, x2, q2),
                       rb * p_pdf->xfxQ2(5, x2, q2),
                       rb * p_pdf->xfxQ2(6, x2, q2)};
    double f_ai_x2[] = {rg * p_pdf->xfxQ2(0, x2, q2),
                        ru * p_pdf->xfxQ2(-1, x2, q2),
                        rd * p_pdf->xfxQ2(-2, x2, q2),
                        rs * p_pdf->xfxQ2(-3, x2, q2),
                        rc * p_pdf->xfxQ2(-4, x2, q2),
                        rb * p_pdf->xfxQ2(-5, x2, q2),
                        rb * p_pdf->xfxQ2(-6, x2, q2)};

    double ZoA = 82.0 / 208.0, NoA = 126.0 / 208.0;
    double u, ub, d, db;
    ////ISOSCALAR NUCLEONS
    if (p_params->isoscalar_projectile)
    {
        u = f_i_x1[1];
        ub = f_ai_x1[1];
        d = f_i_x1[2];
        db = f_ai_x1[2];
        f_i_x1[1] = ZoA * u + NoA * d;
        f_i_x1[2] = ZoA * d + NoA * u;
        f_ai_x1[1] = ZoA * ub + NoA * db;
        f_ai_x1[2] = ZoA * db + NoA * ub;
    }
    if (p_params->isoscalar_target)
    {
        u = f_i_x2[1];
        ub = f_ai_x2[1];
        d = f_i_x2[2];
        db = f_ai_x2[2];
        f_i_x2[1] = ZoA * u + NoA * d;
        f_i_x2[2] = ZoA * d + NoA * u;
        f_ai_x2[1] = ZoA * ub + NoA * db;
        f_ai_x2[2] = ZoA * db + NoA * ub;
    }

    ///* GG->XX
    sum += 0.5 * f_i_x1[0] * f_i_x2[0] * pqcd::diff_sigma::sigma_gg_gg(p_s_hat, p_t_hat, p_u_hat, alphas);
    sum += numFlavours * f_i_x1[0] * f_i_x2[0] * pqcd::diff_sigma::sigma_gg_qaq(p_s_hat, p_t_hat, p_u_hat, alphas);
    //*/
    ///* GQ->XX
    for (int flavour = 1; flavour <= numFlavours; ++flavour)
    {
        sum += f_i_x1[0] * gsl::at(f_i_x2,flavour) * pqcd::diff_sigma::sigma_gq_gq(p_s_hat, p_t_hat, p_u_hat, alphas);
        sum += f_i_x1[0] * gsl::at(f_ai_x2,flavour) * pqcd::diff_sigma::sigma_gq_gq(p_s_hat, p_t_hat, p_u_hat, alphas);
        sum += gsl::at(f_i_x1,flavour) * f_i_x2[0] * pqcd::diff_sigma::sigma_gq_gq(p_s_hat, p_t_hat, p_u_hat, alphas);
        sum += gsl::at(f_ai_x1,flavour) * f_i_x2[0] * pqcd::diff_sigma::sigma_gq_gq(p_s_hat, p_t_hat, p_u_hat, alphas);
    }
    //*/
    ///* QQ->XX
    for (int flavour = 1; flavour <= numFlavours; ++flavour)
    {
        sum += 0.5 * gsl::at(f_i_x1,flavour) * gsl::at(f_i_x2,flavour) * pqcd::diff_sigma::sigma_qiqi_qiqi(p_s_hat, p_t_hat, p_u_hat, alphas);
        sum += 0.5 * gsl::at(f_ai_x1,flavour) * gsl::at(f_ai_x2,flavour) * pqcd::diff_sigma::sigma_qiqi_qiqi(p_s_hat, p_t_hat, p_u_hat, alphas);
    }

    for (int flavour1 = 1; flavour1 <= numFlavours; ++flavour1)
    {
        for (int flavour2 = 1; flavour2 <= numFlavours; ++flavour2)
        {
            if (flavour1 != flavour2)
            {
                sum += gsl::at(f_i_x1,flavour1) * gsl::at(f_i_x2,flavour2) * pqcd::diff_sigma::sigma_qiqj_qiqj(p_s_hat, p_t_hat, p_u_hat, alphas);
                sum += gsl::at(f_ai_x1,flavour1) * gsl::at(f_ai_x2,flavour2) * pqcd::diff_sigma::sigma_qiqj_qiqj(p_s_hat, p_t_hat, p_u_hat, alphas);
                sum += gsl::at(f_i_x1,flavour1) * gsl::at(f_ai_x2,flavour2) * pqcd::diff_sigma::sigma_qiqj_qiqj(p_s_hat, p_t_hat, p_u_hat, alphas);
                sum += gsl::at(f_ai_x1,flavour1) * gsl::at(f_i_x2,flavour2) * pqcd::diff_sigma::sigma_qiqj_qiqj(p_s_hat, p_t_hat, p_u_hat, alphas);
            }
        }
    }

    for (int flavour = 1; flavour <= numFlavours; ++flavour)
    {
        sum += gsl::at(f_i_x1,flavour) * gsl::at(f_ai_x2,flavour) * (pqcd::diff_sigma::sigma_qiaqi_qiaqi(p_s_hat, p_t_hat, p_u_hat, alphas) + 0.5 * pqcd::diff_sigma::sigma_qiaqi_gg(p_s_hat, p_t_hat, p_u_hat, alphas) + (numFlavours - 1) * pqcd::diff_sigma::sigma_qiaqi_qjaqj(p_s_hat, p_t_hat, p_u_hat, alphas));
        sum += gsl::at(f_ai_x1,flavour) * gsl::at(f_i_x2,flavour) * (pqcd::diff_sigma::sigma_qiaqi_qiaqi(p_s_hat, p_t_hat, p_u_hat, alphas) + 0.5 * pqcd::diff_sigma::sigma_qiaqi_gg(p_s_hat, p_t_hat, p_u_hat, alphas) + (numFlavours - 1) * pqcd::diff_sigma::sigma_qiaqi_qjaqj(p_s_hat, p_t_hat, p_u_hat, alphas));
    }
    //*/
    return sum;
}

int pqcd::sigma_jet_integrand(unsigned ndim, const double *p_x, void *p_fdata, unsigned fdim, double *p_fval) noexcept
{
    (void)ndim;
    (void)fdim; //To silence "unused" warnings
    std::tuple<std::shared_ptr<LHAPDF::GridPDF>, const momentum *const, const momentum *const, const pqcd::sigma_jet_params *const> fdata =
        *(static_cast<std::tuple<std::shared_ptr<LHAPDF::GridPDF>, const momentum *const, const momentum *const, const pqcd::sigma_jet_params *const> *>(p_fdata));
    std::shared_ptr<LHAPDF::GridPDF> p_pdf = std::get<0>(fdata);
    const momentum mand_s = *std::get<1>(fdata);
    const momentum p02 = *std::get<2>(fdata);
    auto p_params = std::get<3>(fdata);

    momentum kt2;
    rapidity y1, y2;
    xsectval jacobian;
    pqcd::scale_limits_from_0_1(p_x[0], p_x[1], p_x[2], p02, mand_s, &kt2, &y1, &y2, &jacobian);
    const auto sqrt_s_per_kt = sqrt(mand_s / kt2);

    momentum fac_scale;

    switch (p_params->scale_c)
    {
    case scaled_from_kt:
        fac_scale = pow(p_params->scalar, 2) * kt2;
        break;
    case constant:
        fac_scale = pow(p_params->scalar, 2);
        break;
    default:
        fac_scale = kt2;
        break;
    }

    const auto x1 = (exp(y1) + exp(y2)) / sqrt_s_per_kt;
    const auto x2 = (exp(-y1) + exp(-y2)) / sqrt_s_per_kt;

    const auto s_hat = pqcd::s_hat_from_ys(&y1, &y2, &kt2);
    const auto t_hat = pqcd::t_hat_from_ys(&y1, &y2, &kt2);
    const auto u_hat = pqcd::u_hat_from_ys(&y1, &y2, &kt2);

    //SES
    if (p_params->use_ses)
    {
        const auto subprocess_cs = pqcd::diff_sigma::sigma_gg_gg(&s_hat, &t_hat, &u_hat, p_pdf->alphasQ2(fac_scale));
        p_fval[0] = 0.5 * pqcd::f_ses(&x1, &fac_scale, p_pdf) * pqcd::f_ses(&x2, &fac_scale, p_pdf) * subprocess_cs * jacobian * 10 / pow(FMGEV, 2); //UNITS: mb
    }
    else //FULL SUMMATION
    {
        p_fval[0] = 0.5 * (pqcd::diff_sigma::sigma_jet(&x1, &x2, &fac_scale, p_pdf, &s_hat, &t_hat, &u_hat, p_params->p_d_params) + pqcd::diff_sigma::sigma_jet(&x2, &x1, &fac_scale, p_pdf, &s_hat, &u_hat, &t_hat, p_params->p_d_params)) * jacobian * 10 / pow(FMGEV, 2); //UNITS: mb
    }

    return 0; // success
}

xsectval pqcd::f_ses(const double *const p_x, const double *const p_q2, std::shared_ptr<LHAPDF::GridPDF> p_pdf) noexcept
{
    xsectval sum = 0;
    for (int flavour = -5; flavour <= 5; ++flavour)
    {
        sum += (flavour == 0) ? p_pdf->xfxQ2(flavour, *p_x, *p_q2) : 4 * p_pdf->xfxQ2(flavour, *p_x, *p_q2) / 9;
    }
    //cout<<"x="<<*p_x<<", q2="<<*p_q2<<" SES="<<sum<<endl;
    return sum;
}

void pqcd::scale_limits_from_0_1(const rapidity & z1, const rapidity & z2, const rapidity & z3,                     //variables between 0 and 1
                                           const momentum & kt2_lower_cutoff, const momentum & mand_s,                               //parameters
                                           momentum *const p_kt2, rapidity *const p_y1, rapidity *const p_y2, xsectval *const p_jacobian) noexcept //output
{
    *p_kt2 = kt2_lower_cutoff + z1 * ((mand_s / 4) - kt2_lower_cutoff);
    const auto sqrt_s_per_kt = sqrt(mand_s / *p_kt2);

    const auto y1_upper = acosh(sqrt_s_per_kt / 2);
    //const auto y1_lower = -y1_upper;
    *p_y1 = (-1 + (2 * z2)) * y1_upper; //y1_lower + z2 * (y1_upper - y1_lower)

    const auto y2_upper = log(sqrt_s_per_kt - exp(*p_y1));
    const auto y2_lower = -log(sqrt_s_per_kt - exp(-*p_y1));
    *p_y2 = y2_lower + z3 * (y2_upper - y2_lower);

    *p_jacobian = ((mand_s / 4) - kt2_lower_cutoff) * (2 * y1_upper) * (y2_upper - y2_lower);
}

momentum pqcd::s_hat_from_ys(const rapidity *const p_y1, const rapidity *const p_y2, const momentum *const p_kt2) noexcept
{
    const rapidity y1 = *p_y1, y2 = *p_y2;
    const momentum kt2 = *p_kt2;
    return 2 * kt2 * (1 + cosh(y1 - y2));
}

momentum pqcd::t_hat_from_ys(const rapidity *const p_y1, const rapidity *const p_y2, const momentum *const p_kt2) noexcept
{
    const rapidity y1 = *p_y1, y2 = *p_y2;
    const momentum kt2 = *p_kt2;
    return -kt2 * (1 + exp(-y1 + y2));
}

momentum pqcd::u_hat_from_ys(const rapidity *const p_y1, const rapidity *const p_y2, const momentum *const p_kt2) noexcept
{
    const rapidity y1 = *p_y1, y2 = *p_y2;
    const momentum kt2 = *p_kt2;
    return -kt2 * (1 + exp(y1 - y2));
}
