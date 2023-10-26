//Copyright (c) 2023 Mikko Kuha

#ifndef HIGH_LEVEL_CALCS_HPP
#define HIGH_LEVEL_CALCS_HPP

#include <algorithm>
#include <fstream>
#include <future>
#include <gsl/gsl_multimin.h>
#include <iostream>
#include <numeric>
#include <memory>
#include <omp.h>
#include <optional>
#include <sstream>
#include <variant>
#include <tuple>


#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshorten-64-to-32"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "LHAPDF/GridPDF.h"
#pragma GCC diagnostic pop

#include "generic_helpers.hpp"
#include "linear_interpolator.hpp"
#include "linterp.h"
#include "nucleus_generator.hpp"
#include "pqcd.hpp"
#include "Tpp_builder.hpp"
#include "typedefs.hpp"


struct AA_collision_params
{
  bool mc_glauber_mode;
  bool pp_scattering;
  bool pA_scattering;
  bool spatial_pdfs;
  bool calculate_end_state;
  bool use_nn_b2_max;
  double sigma_inel;
  Tpp_builder *const Tpp;
  B2_normalization_mode normalize_to;
  double sqrt_s;
  double energy_threshold;
  double nn_b2_max;
  double T_AA_0;
};

using variant_sigma_jet = std::variant<InterpMultilinear<3, double>, InterpMultilinear<2, double>, linear_interpolator, double>;
using variant_envelope_pars = std::variant<linear_interpolator, envelope_func>;

class calcs
{
public:
    /**
     * @brief Finds the maximum of the differential sigma_jet in terms of pT, y1 and y2. 
     * 
     * @param kt
     * @param sqrt_s 
     * @param p_pdf 
     * @param params 
     * @return auto 
     */
    static auto find_max_dsigma
    (
        const double &kt, 
        const double &sqrt_s,
        std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
        pqcd::sigma_jet_params params,
        const double &ave_T_AA_0
    ) noexcept -> std::tuple<double,double>
    {
        double max_dsigma;
        double error_est;

        if (params.d_params.npdfs_spatial)
        {
            if (params.snPDFs_linear)
            {
                //r=(1+cT)
                const uint_fast16_t A = params.d_params.A,
                                    B = params.d_params.B;
                //c=A*(R-1)/TAA(0)
                const double scaA = static_cast<double>(A) * 3.0, 
                            intA = 1.0 - scaA;
                const std::function<double(double const&)> 
                    rA_spatial_ = [&](double const &r)
                        {
                            auto dummy = r*scaA + intA;
                            return (dummy < 1.0 ) ? 1.0 : dummy; // looking for max: only antishadowing allowed
                        }; //r_s=1+c*sum(Tpp)

                const double scaB = static_cast<double>(B) * 3.0, 
                            intB = 1.0 - scaB;
                const std::function<double(double const&)> 
                    rB_spatial_ = [&](double const &r)
                        {
                            auto dummy = r*scaB + intB;
                            return (dummy < 1.0 ) ? 1.0 : dummy; // looking for max: only antishadowing allowed
                        };

                params.d_params.rA_spatial = rA_spatial_;
                params.d_params.rB_spatial = rB_spatial_;
            }
            else
            {
                //r=exp(cT)

                const std::function<double(double const&)> 
                    rA_spatial_ = [c_func=params.c_A_func](double const &r)
                        {
                            if (r>0.0)
                            {
                                double c = c_func.value_at(r);
                                auto dummy = std::exp(c * 0.4);
                                return (dummy < 1.0 ) ? 1.0 : dummy; // looking for max: only antishadowing allowed
                            }
                            else
                            {
                                return 0.0;
                            }
                        };

                const std::function<double(double const&)> 
                    rB_spatial_ = [c_func=params.c_A_func](double const &r)
                        {
                            if (r>0.0)
                            {
                                double c = c_func.value_at(r);
                                auto dummy = std::exp(c * 0.4);
                                return (dummy < 1.0 ) ? 1.0 : dummy; // looking for max: only antishadowing allowed
                            }
                            else
                            {
                                return 0.0;
                            }
                        };

                params.d_params.rA_spatial = rA_spatial_;
                params.d_params.rB_spatial = rB_spatial_;
            }
        }

        struct f_params fparams = {kt, sqrt_s, p_pdf, params};
        //Use Simplex algorithm by Nelder and Mead to find the maximum of dsigma
        const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
        gsl_multimin_fminimizer *minimizer = nullptr;
        gsl_vector *ys, *init_step_size;
        gsl_multimin_function minfunc;
        ys = gsl_vector_alloc(2);
        gsl_vector_set_all(ys, 0);
        init_step_size = gsl_vector_alloc(2);
        gsl_vector_set_all(init_step_size, 0.5);

        uint_fast8_t iter = 0;
        int status;
        double simplex_size = 0;

        minfunc.n = 2;
        minfunc.f = fdsigma;
        minfunc.params = &fparams;

        minimizer = gsl_multimin_fminimizer_alloc(T, 2);
        gsl_multimin_fminimizer_set(minimizer, &minfunc, ys, init_step_size);
        do
        {
            iter++;
            status = gsl_multimin_fminimizer_iterate(minimizer);

            // The error code GSL_ENOPROG signifies that the minimizer is unable 
            // to improve on its current estimate, either due to numerical 
            // difficulty or because a genuine local minimum has been reached.
            if (status == GSL_ENOPROG)
            {
                break;
            }

            simplex_size = gsl_multimin_fminimizer_size(minimizer);
            status = gsl_multimin_test_size(simplex_size, 1e-3);

        } while (status == GSL_CONTINUE && iter < 100);

        max_dsigma = static_cast<double>(-minimizer->fval);
        error_est = static_cast<double>(simplex_size);

        gsl_vector_free(ys);
        gsl_vector_free(init_step_size);
        gsl_multimin_fminimizer_free(minimizer);

        return std::make_tuple(std::move(max_dsigma), std::move(error_est));
    }
 
    /**
     * @brief Calculates the envelope function, its primitive and its primitives
     * inverse for dsigma_jet/dkt.
     * The envelope is of the form
     * f(kt) = norm1/kt             if  kt <= switch_kt,
     * f(kt) = norm2*kt^power       if  kt >  switch_kt,
     * if kt0 < 2.0. If kt0 >= 2.0, the ~1/kt part is not needed.
     * 
     * @param kt0 k_T cutoff value for the dsigma (GeV)
     * @param sqrt_s sqrt(s) for the dsigma (GeV)
     * @param p_pdf LHAPDF PDF object to calculate the dsigma with
     * @param jet_params collection of parameters for the dsigma 
     * @return auto an object holding the envelope function and its inverse
     */
    static auto calculate_envelope_params
    (
        const double &kt0, 
        const double &sqrt_s,
        std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
        pqcd::sigma_jet_params jet_params,
        const double &ave_T_AA_0
    ) noexcept
    {
        // How tight we want the envelope to be, lower values == faster but more prone to error
        double extra = 1.05;

        double env_min_kt=kt0;
        double env_norm1=0;
        double env_norm2=0;
        double env_power=0;
        double env_switch_kt=0;
        double env_prim_integ_constant{0};
        double env_prim_switch_y{0};
        std::function<double(const double&)> env_func;
        std::function<double(const double&)> env_prim;
        std::function<double(const double&)> env_prim_inv;


        if (kt0 < 2.0)
        {
            //Calculate the normalization for the ~1/kt part
            double dummy_kt0 = 1.0;
            auto [max_dsigma, err] = calcs::find_max_dsigma(dummy_kt0, sqrt_s, p_pdf, jet_params, ave_T_AA_0);
            env_norm1 = (max_dsigma + fabs(err)) * extra;

            //Calculate parameters for the a*kt^b part
            double kt1 = 2.0;
            double kt2 = 3.0;
            auto [max_dsigma1, err1] = calcs::find_max_dsigma(kt1, sqrt_s, p_pdf, jet_params, ave_T_AA_0);
            auto [max_dsigma2, err2] = calcs::find_max_dsigma(kt2, sqrt_s, p_pdf, jet_params, ave_T_AA_0);
            double logkt1 = std::log(kt1);
            double logkt2 = std::log(kt2);
            double logy1 = std::log((max_dsigma1 + fabs(err1)) * extra);
            double logy2 = std::log((max_dsigma2 + fabs(err2)) * extra);
            env_norm2 = std::exp(-(logkt2*logy1 - logkt1*logy2)/(logkt1 - logkt2)); 
            env_power = (logy1 - logy2)/(logkt1 - logkt2); 

            //Calculate the location of the knee
            env_switch_kt = std::pow(env_norm2/env_norm1, -1/(1+env_power));

            env_func = [ktx=env_switch_kt, A=env_norm1, a=env_norm2, b=env_power](const double &kt)
            {
                if (kt <= ktx)
                {
                    return A/kt;
                }
                else
                {
                    return a*std::pow(kt,b);
                }
            };

            //Calculate the primitive for the envelope
            //The integration constant to make the primitive continuous
            auto dummy_a = (env_norm2/(1+env_power))*std::pow(env_switch_kt,1+env_power);
            auto dummy_b = env_norm1*std::log(env_switch_kt/env_min_kt);
            env_prim_integ_constant = dummy_a - dummy_b;

            env_prim = [ktx=env_switch_kt, A=env_norm1, a=env_norm2, b=env_power, lb=env_min_kt, c=env_prim_integ_constant]
            (const double &kt)
            {
                if (kt <= ktx)
                {
                    return A*std::log(kt/lb);
                }
                else
                {
                    return (a/(1+b))*std::pow(kt,1+b) - c;
                }
            };

            //Calculate the inverse of the primitive
            env_prim_switch_y = env_prim(env_switch_kt);

            env_prim_inv = [sx=env_prim_switch_y, A=env_norm1, a=env_norm2, b=env_power, lb=env_min_kt, c=env_prim_integ_constant]
            (const double &s)
            {
                if (s <= sx)
                {
                    return lb*std::exp(s/A);
                }
                else
                {
                    return std::pow((c+s)*(1+b)/a, 1.0/(1+b));
                }
            };
        }
        else //If kt0>2.0, we don't need the ~1/kt part
        {
            //Calculate parameters for the a*kt^b part
            double kt1 = env_min_kt;
            double kt2 = env_min_kt + 1.0;
            auto [max_dsigma1, err1] = calcs::find_max_dsigma(kt1, sqrt_s, p_pdf, jet_params, ave_T_AA_0);
            auto [max_dsigma2, err2] = calcs::find_max_dsigma(kt2, sqrt_s, p_pdf, jet_params, ave_T_AA_0);
            double logkt1 = std::log(kt1);
            double logkt2 = std::log(kt2);
            double logy1 = std::log((max_dsigma1 + fabs(err1)) * extra);
            double logy2 = std::log((max_dsigma2 + fabs(err2)) * extra);
            env_norm2 = std::exp(-(logkt2*logy1 - logkt1*logy2)/(logkt1 - logkt2)); 
            env_power = (logy1 - logy2)/(logkt1 - logkt2); 

            env_func = [a=env_norm2, b=env_power](const double &kt)
            {
                return a*std::pow(kt,b);
            };

            //Calculate the primitive for the envelope
            //The integration constant so that G(kt0)=0
            env_prim_integ_constant = (env_norm2/(1+env_power))*std::pow(env_min_kt,1+env_power);

            env_prim = [a=env_norm2, b=env_power, c=env_prim_integ_constant](const double &kt)
            {
                return (a/(1+b))*std::pow(kt,1+b) - c;
            };

            //Calculate the inverse of the primitive
            env_prim_inv = [a=env_norm2, b=env_power, c=env_prim_integ_constant](const double &s)
            {
                return std::pow((c+s)*(1+b)/a, 1.0/(1+b));
            };
        }

        return envelope_func
            {
                env_min_kt,
                env_norm1,
                env_norm2,
                env_power,
                env_switch_kt,
                env_prim_integ_constant,
                env_prim_switch_y,
                env_func,
                env_prim,
                env_prim_inv
            };
    }
    
    /*
    * Struct needed by find_max_dsigma
    */
    struct f_params
    {
        const double &kt;
        const double &sqrt_s;
        std::shared_ptr<LHAPDF::GridPDF> p_pdf;
        pqcd::sigma_jet_params sigma_params;
    };

    /*
    *Function of diff_cross_section_2jet for maximization with respect to y1 and y2
    */
    static double fdsigma
    (
        const gsl_vector *v, 
        void *params
    )
    {
        auto *fparams = static_cast<struct f_params *>(params);
        auto kt = fparams->kt;
        auto sqrt_s = fparams->sqrt_s;
        auto p_pdf = fparams->p_pdf;
        auto sigma_params = fparams->sigma_params;
        double y1, y2;
        y1 = gsl_vector_get(v, 0);
        y2 = gsl_vector_get(v, 1);

        auto xsection = pqcd::diff_cross_section_2jet(sqrt_s, kt, y1, y2, p_pdf, sigma_params, false, false);
        double total_xsection = 0;
        for (auto xsect : xsection)
        {
            total_xsection += xsect.sigma;
        }

        return -total_xsection;
    }

    static auto prepare_sigma_jets
    (
        const bool read_sigmajets_from_file,
        std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
        const double &mand_s,
        const double &sqrt_s,
        const double &kt02, 
        const double &kt0, 
        const double &power_law,
        pqcd::sigma_jet_params jet_params,
        const double &ave_T_AA_0,
        const std::string &s_jet_file_name
    ) ->
    std::tuple
    <
        double,
        variant_sigma_jet,
        std::optional<std::vector<double> >,
        variant_envelope_pars,
        std::optional<std::vector<double> >
    >
    {
        const bool use_npdfs = (jet_params.d_params.projectile_with_npdfs 
                                || jet_params.d_params.target_with_npdfs);
        const bool spatial_pdfs = jet_params.d_params.npdfs_spatial;
        double dijet_norm = 0;

        if //proton PDFs
        (
            !use_npdfs && !spatial_pdfs
        )
        {
            std::cout<<"Calculating sigma_jet..."<<std::flush;
            double sigma_jet = pqcd::calculate_sigma_jet(p_pdf, &mand_s, &kt02, jet_params);
            std::cout<<"done!"<<std::endl;

            dijet_norm = sigma_jet;

            std::cout<<"Calculating envelope..."<<std::flush;
            auto env_params = calcs::calculate_envelope_params(kt0, sqrt_s, p_pdf, jet_params, ave_T_AA_0);
            std::cout<<"done!"<<std::endl;

            return std::make_tuple
                (
                    dijet_norm,
                    variant_sigma_jet(sigma_jet),
                    std::nullopt,
                    variant_envelope_pars(env_params),
                    std::nullopt
                );
        }
        else if //nPDFS in use
        (
            !spatial_pdfs
        )
        {
            std::cout<<"Calculating sigma_jet..."<<std::flush;
            double sigma_jet = pqcd::calculate_sigma_jet(p_pdf, &mand_s, &kt02, jet_params);
            std::cout<<"done!"<<std::endl;

            dijet_norm = sigma_jet;

            std::cout<<"Calculating envelope..."<<std::flush;
            auto env_params = calcs::calculate_envelope_params(kt0, sqrt_s, p_pdf, jet_params, ave_T_AA_0);
            std::cout<<"done!"<<std::endl;

            return std::make_tuple
                (
                    dijet_norm,
                    variant_sigma_jet(sigma_jet),
                    std::nullopt,
                    variant_envelope_pars(env_params),
                    std::nullopt
                );
        }
        else //spatial nPDFs
        {
            if (read_sigmajets_from_file)
            {
                std::cout<<"Reading spatial sigma_jets..."<<std::flush;
                variant_sigma_jet sigma_jet = calcs::read_sigma_jets_spatial(s_jet_file_name, jet_params.d_params.K_factor);
                std::cout<<"done!"<<std::endl;

                auto other_params = jet_params;
                other_params.d_params.npdfs_spatial = false;
                dijet_norm = pqcd::calculate_sigma_jet(p_pdf, &mand_s, &kt02, other_params); //sigma_jet with ave nPDF
                //dijet_norm = 93.7604;
                std::cout<<"dijet_norm = "<<dijet_norm<<std::endl;

                std::cout<<"Calculating envelope..."<<std::flush;
                auto env_params = calcs::calculate_envelope_params(kt0, sqrt_s, p_pdf, jet_params, ave_T_AA_0);
                std::cout<<"done!"<<std::endl;

                return std::make_tuple
                    (
                        dijet_norm,
                        std::move(sigma_jet),
                        std::nullopt,
                        variant_envelope_pars(env_params),
                        std::nullopt
                    );
            }
            else
            {
                double tolerance=0.05,
                    upper_sumTpp_limit=/*0.21033,*/0.44,//0.61, 
                    lower_sumTpp_limit=0.0;//0.01;

                std::cout<<"Calculating spatial sigma_jets..."<<std::endl;
                variant_sigma_jet sigma_jet 
                    = calcs::calculate_spatial_sigma_jets
                        (
                            s_jet_file_name,
                            tolerance, 
                            p_pdf, 
                            mand_s, 
                            kt02, 
                            jet_params, 
                            upper_sumTpp_limit, 
                            lower_sumTpp_limit,
                            ave_T_AA_0
                        );
                std::cout<<"done!"<<std::endl;


                std::cout<<"Calculating dijet_norm..."<<std::endl;
                auto other_params = jet_params;
                other_params.d_params.npdfs_spatial = false;
                dijet_norm = pqcd::calculate_sigma_jet(p_pdf, &mand_s, &kt02, other_params); //sigma_jet with ave nPDF
                //dijet_norm = 93.7604;
                std::cout<<"dijet_norm = "<<dijet_norm<<std::endl;

                std::cout<<"Calculating envelope..."<<std::flush;
                auto env_params = calcs::calculate_envelope_params(kt0, sqrt_s, p_pdf, jet_params, ave_T_AA_0);
                std::cout<<"done!"<<std::endl;

                return std::make_tuple
                    (
                        dijet_norm,
                        std::move(sigma_jet),
                        std::nullopt,
                        variant_envelope_pars(env_params),
                        std::nullopt
                    );
            }
        }
    }

    //collide_nuclei<false>
    //template<bool temp_bool, typename std::enable_if <!temp_bool> :: type* = nullptr>
    static auto collide_nuclei_no_snPDF
    (
        std::vector<nucleon> &pro, 
        std::vector<nucleon> &tar, 
        std::vector<nn_coll> &binary_collisions, 
        variant_sigma_jet &sigma_jets, 
        std::uniform_real_distribution<double> unirand, 
        std::shared_ptr<std::mt19937> eng, 
        const AA_collision_params &AA_params,
        const pqcd::sigma_jet_params &dsigma_params,
        const double &kt0,
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf,
        const double &power_law,
        variant_envelope_pars &env_func,
        const bool &verbose
    ) noexcept -> void
    {
        uint_fast32_t n_pairs = 0, mombroke = 0, skipped=0, nof_softs = 0;
        
        std::vector<std::tuple<nucleon* const, nucleon* const> > binary_pairs;
        
        for (auto &A : pro)
        {
            for (auto &B : tar)
            {
                binary_pairs.push_back(std::make_pair(&A, &B));
            }
        }
        
        std::vector<uint_fast64_t> pair_indexes(binary_pairs.size());
        std::iota(pair_indexes.begin(), pair_indexes.end(), 0);
        
        std::shuffle(pair_indexes.begin(), pair_indexes.end(), *eng);
        
        for (auto & ind : pair_indexes)
        {
            auto & [A, B] = binary_pairs[ind];

            if (AA_params.use_nn_b2_max && A->calculate_bsquared(*B) > AA_params.nn_b2_max)
            {
                skipped++;
                continue;
            }
            
            n_pairs++;
            if ((A->mom < AA_params.energy_threshold) || (B->mom < AA_params.energy_threshold))
            {
                mombroke++;
                continue;
            }
            
            nn_coll newpair(A, B, 2 * sqrt(A->mom * B->mom));
            
            if (AA_params.mc_glauber_mode)
            {
                // "ball" diameter = distance at which two nucleons interact
                const double d2 = AA_params.sigma_inel/(M_PI*10); // in fm^2
                const double dij2 = newpair.getcr_bsquared();
                
                if (dij2 > d2) //no collision
                {
                    continue;
                }
                //collision
                if (AA_params.calculate_end_state)
                {
                    double sigma_jet;
                    // if (AA_params.reduce_nucleon_energies)
                    // {
                    //     sigma_jet = std::get<linear_interpolator>(sigma_jets).value_at(pow(newpair.getcr_sqrt_s(), 2));
                    //     double envelope_maximum = std::get<linear_interpolator>(env_func).value_at(newpair.getcr_sqrt_s());
                    
                    //     pqcd::generate_bin_NN_coll
                    //     (
                    //         newpair, 
                    //         sigma_jet, 
                    //         AA_params.Tpp(newpair.getcr_bsquared()), 
                    //         kt0,
                    //         unirand, 
                    //         eng,
                    //         p_p_pdf,
                    //         dsigma_params,
                    //         power_law,
                    //         envelope_maximum
                    //     );
                    // }
                    // else //Single sigma_jet
                    {
                        sigma_jet = std::get<double>(sigma_jets);
                        auto env_func_ = std::get<envelope_func>(env_func);

                        std::poisson_distribution<uint_fast16_t> dist(sigma_jet*AA_params.Tpp->at(newpair));
                        uint_fast16_t nof_dijets = dist(*eng);
                        newpair.dijets.reserve(nof_dijets);                        

                        const double sqrt_s = newpair.getcr_sqrt_s();

                        for (uint_fast16_t i=0; i < nof_dijets; i++)
                        {
                            newpair.dijets.push_back(pqcd::generate_2_to_2_scatt
                            (
                                sqrt_s,
                                kt0,
                                sqrt_s / 2.0,
                                unirand,
                                eng,
                                p_p_pdf,
                                dsigma_params,
                                power_law,
                                env_func_,
                                B->is_neutron,
                                A->is_neutron
                            ));
                        }

                        // pqcd::generate_bin_NN_coll
                        // (
                        //     newpair, 
                        //     sigma_jet, 
                        //     AA_params.Tpp(newpair.getcr_bsquared()), 
                        //     kt0,
                        //     unirand, 
                        //     eng,
                        //     p_p_pdf,
                        //     dsigma_params,
                        //     power_law,
                        //     env_func_,
                        //     B->is_neutron,
                        //     A->is_neutron
                        // );
                    }
                    
                    // if (AA_params.reduce_nucleon_energies)
                    // {
                    //     newpair.reduce_energy_and_push_end_states_to_collider_frame();
                    // }
                    // else
                    {
                        newpair.push_end_states_to_collider_frame();
                    }
                }
                newpair.wound();
                binary_collisions.push_back(std::move(newpair));
            }
            else
            {
                double sigma_jet;
                // if (AA_params.reduce_nucleon_energies)
                // {
                //     sigma_jet = std::get<linear_interpolator>(sigma_jets).value_at(pow(newpair.getcr_sqrt_s(), 2));
                // }
                // else //Single sigma_jet
                {
                    sigma_jet = std::get<double>(sigma_jets);
                }
                
                // newpair.calculate_xsects(sigma_jet, AA_params.Tpp, newpair.getcr_bsquared());
                // auto ran = unirand(*eng);
                // switch (AA_params.normalize_to)
                // {
                // case B2_normalization_mode::total:
                //     ran *= newpair.getcr_max_tot_xsect();
                //     break;
                
                // case B2_normalization_mode::inelastic:
                //     ran *= newpair.getcr_max_inel_xsect();
                //     break;
                
                // case B2_normalization_mode::nothing:
                //     break;
                
                // default:
                //     ran *= newpair.getcr_max_inel_xsect();
                //     break;
                // }

                // if (ran > newpair.getcr_effective_inel_xsect())
                // {
                //     if (ran > newpair.getcr_effective_tot_xsect())
                //     {
                //         continue;
                //     }
                //     nof_softs++;
                //     continue;
                // }
                if (AA_params.calculate_end_state)
                {
                    // double envelope_maximum;
                    // if (AA_params.reduce_nucleon_energies)
                    // {
                    //     double envelope_maximum = std::get<linear_interpolator>(env_func).value_at(newpair.getcr_sqrt_s());
                    
                    //     pqcd::generate_bin_NN_coll
                    //     (
                    //         newpair, 
                    //         sigma_jet, 
                    //         AA_params.Tpp(newpair.getcr_bsquared()), 
                    //         kt0,
                    //         unirand, 
                    //         eng,
                    //         p_p_pdf,
                    //         dsigma_params,
                    //         power_law,
                    //         envelope_maximum
                    //     );
                    // }
                    // else //Single sigma_jet
                    {
                        auto env_func_ = std::get<envelope_func>(env_func);

                        std::poisson_distribution<uint_fast8_t> dist(sigma_jet*AA_params.Tpp->at(newpair));
                        uint_fast8_t nof_dijets = dist(*eng);
                        newpair.dijets.reserve(nof_dijets);

                        const double sqrt_s = newpair.getcr_sqrt_s();

                        for (uint_fast8_t i=0; i < nof_dijets; i++)
                        {
                            newpair.dijets.push_back(pqcd::generate_2_to_2_scatt
                            (
                                sqrt_s,
                                kt0,
                                sqrt_s / 2.0,
                                unirand,
                                eng,
                                p_p_pdf,
                                dsigma_params,
                                power_law,
                                env_func_,
                                B->is_neutron,
                                A->is_neutron
                            ));
                        }

                        // pqcd::generate_bin_NN_coll
                        // (
                        //     newpair, 
                        //     sigma_jet, 
                        //     AA_params.Tpp(newpair.getcr_bsquared()), 
                        //     kt0,
                        //     unirand, 
                        //     eng,
                        //     p_p_pdf,
                        //     dsigma_params,
                        //     power_law,
                        //     env_func_,
                        //     B->is_neutron,
                        //     A->is_neutron
                        // );
                    }
                    
                    // if (AA_params.reduce_nucleon_energies)
                    // {
                    //     newpair.reduce_energy_and_push_end_states_to_collider_frame();
                    // }
                    // else
                    {
                        newpair.push_end_states_to_collider_frame();
                    }
                }
                newpair.wound();
                binary_collisions.push_back(std::move(newpair));
            }
        }

        if (verbose)
        {
            std::cout << "Bruteforced " << n_pairs << " pairs, got " << binary_collisions.size()+nof_softs << " collisions, of which softs "<< nof_softs<< " and hards "<< binary_collisions.size()<<" , double threshold broke " << mombroke << " times, skipped "<< skipped << " pairs that were too far apart" << std::endl;
        }
    }

    //collide_nuclei<true>
    //template<bool temp_bool, typename std::enable_if <temp_bool> :: type* = nullptr>
    static auto collide_nuclei_snPDF
    (
        std::vector<nucleon> &pro, 
        std::vector<nucleon> &tar, 
        std::vector<nn_coll> &binary_collisions, 
        variant_sigma_jet &sigma_jets,
        std::uniform_real_distribution<double> unirand, 
        std::shared_ptr<std::mt19937> eng, 
        const AA_collision_params &AA_params,
        pqcd::sigma_jet_params dsigma_params,
        const double &kt0,
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf,
        const double &power_law,
        variant_envelope_pars &env_func,
        const bool &verbose,
        const double &M_factor,
        const double &proton_width,
        const uint_fast16_t &is_sat_y_dep
    ) noexcept -> void
    {

        uint_fast32_t n_pairs = 0, mombroke = 0, skipped=0, nof_softs = 0;

        if (AA_params.pp_scattering)
        {
            std::cout<<"SPATIAL NPDFS REQUESTED IN PROTON-PROTON COLLISION!!!"<<std::endl;
            exit(1);
        }
        else if (AA_params.pA_scattering)
        {
            std::cout<<"SPATIAL NPDFS REQUESTED IN PROTON-NUCLEUS COLLISION!!! (NOT IMPLEMENTED)"<<std::endl;
            exit(1);
        }

        double tAA_0 = AA_params.T_AA_0;
        double tBB_0 = AA_params.T_AA_0;
        
        if (verbose)
        {
            std::cout << "T_AA(0)= " << tAA_0 << ", T_BB(0)= " << tBB_0 << std::endl;
        }

        std::vector<std::tuple
        <
            std::tuple<nucleon* const, const double* const>, 
            std::tuple<nucleon* const, const double* const> 
        > > binary_pairs;

        std::vector
        <
            std::tuple<nucleon* const, const double>
        > pro_spatial, tar_spatial;

        for (auto &A : pro)
        {
            pro_spatial.push_back(std::make_pair(&A, AA_params.Tpp->calculate_sum_tpp(A, pro)));
        }

        for (auto &B : tar)
        {
            tar_spatial.push_back(std::make_pair(&B, AA_params.Tpp->calculate_sum_tpp(B, tar)));
        }

        for (auto & [A, sum_A] : pro_spatial)
        {
            for (auto & [B, sum_B] : tar_spatial)
            {
                binary_pairs.push_back(std::make_pair(std::make_pair(A, &sum_A), std::make_pair(B, &sum_B)));
            }
        }
        
        std::vector<uint_fast64_t> pair_indexes(binary_pairs.size());
        std::iota(pair_indexes.begin(), pair_indexes.end(), 0);

        std::shuffle(pair_indexes.begin(), pair_indexes.end(), *eng);

        for (auto & ind : pair_indexes)
        {
            auto & [A_pair, B_pair] = binary_pairs[ind];

            auto & [A, sum_tppa] = A_pair;
            auto & [B, sum_tppb] = B_pair;

            if (AA_params.use_nn_b2_max && A->calculate_bsquared(*B) > AA_params.nn_b2_max)
            {
                skipped++;
                continue;
            }

            n_pairs++;
            if ((A->mom < AA_params.energy_threshold) || (B->mom < AA_params.energy_threshold))
            {
                mombroke++;
                continue;
            }

            nn_coll new_pair(A, B, 2 * sqrt(A->mom * B->mom));
            
            if (AA_params.mc_glauber_mode)
            {
                // "ball" diameter = distance at which two nucleons interact
                const double d2 = AA_params.sigma_inel/(M_PI*10); // in fm^2
                const double dij2 = new_pair.getcr_bsquared();
                
                if (dij2 > d2) //no collision
                {
                    continue;
                }
                //collision
                if (AA_params.calculate_end_state)
                {
                    double sigma_jet;
                    // if (AA_params.reduce_nucleon_energies)
                    // {
                    //     array<double,3> args{*sum_tppa, *sum_tppb, pow(newpair.getcr_sqrt_s(), 2)};
                    //     sigma_jet = std::get<InterpMultilinear<3, double> >(sigma_jets).interp(args.begin());
                    // }
                    // else //Single sigma_jet
                    {
                        array<double,2> args{*sum_tppa, *sum_tppb};
                        sigma_jet = std::get<InterpMultilinear<2, double> >(sigma_jets).interp(args.begin());
                        //pqcd::calculate_spatial_sigma_jet(p_p_pdf, p_n_pdf, &mand_s, &kt02, &jet_params, &sum_tppa, &sum_tppb, &tAA_0, &tBB_0);
                    }
                    
                    if (dsigma_params.snPDFs_linear)
                    {
                        //r=(1+cT)

                        const uint_fast16_t NA = dsigma_params.d_params.A, 
                                            NB = dsigma_params.d_params.B;
                                            
                        const auto spatial_cutoff = dsigma_params.d_params.spatial_cutoff;
                        //c=A*(R-1)/TAA(0)
                        const double scaA = static_cast<double>(NA) * *sum_tppa / tAA_0, 
                                    intA = 1.0 - scaA;
                        const std::function<double(double const&)> 
                            rA_spatial_ = [&,co=spatial_cutoff](double const &r)
                                {
                                    auto dummy = r*scaA + intA;
                                    // return dummy; // no cutoff
                                    // return (dummy < 0.0 ) ? 0.0 : dummy; // cutoff at 1+cT < 0
                                    return (dummy < co ) ? co : dummy; // cutoff at 1+cT < spatial_cutoff
                                    // return (dummy < 1/static_cast<double>(NA) ) ? 1/static_cast<double>(NA) : dummy; // cutoff at 1+cT < A
                                    
                                }; //r_s=1+c*sum(Tpp)

                        const double scaB = static_cast<double>(NB) * *sum_tppb / tBB_0, 
                                    intB = 1.0 - scaB;
                        const std::function<double(double const&)> 
                            rB_spatial_ = [&,co=spatial_cutoff](double const &r)
                                {
                                    auto dummy = r*scaB + intB;
                                    // return dummy; // no cutoff
                                    // return (dummy < 0.0 ) ? 0.0 : dummy; // cutoff at 1+cT < 0
                                    return (dummy < co ) ? co : dummy; // cutoff at 1+cT < spatial_cutoff
                                    // return (dummy < 1/static_cast<double>(NB) ) ? 1/static_cast<double>(NB) : dummy; // cutoff at 1+cT < B
                                };
                                
                        dsigma_params.d_params.rA_spatial = rA_spatial_;
                        dsigma_params.d_params.rB_spatial = rB_spatial_;
                    }
                    else
                    {
                        //r=exp(cT)

                        const std::function<double(double const&)> 
                            rA_spatial_ = [TAi=*sum_tppa,c_func=dsigma_params.c_A_func](double const &r)
                                {
                                    if (r>0.0)
                                    {
                                        double c = c_func.value_at(r);
                                        return std::exp(c * TAi);
                                    }
                                    else
                                    {
                                        return 0.0;
                                    }
                                };

                        const std::function<double(double const&)> 
                            rB_spatial_ = [TBi=*sum_tppb,c_func=dsigma_params.c_A_func](double const &r)
                                {
                                    if (r>0.0)
                                    {
                                        double c = c_func.value_at(r);
                                        return std::exp(c * TBi);
                                    }
                                    else
                                    {
                                        return 0.0;
                                    }
                                };

                        dsigma_params.d_params.rA_spatial = rA_spatial_;
                        dsigma_params.d_params.rB_spatial = rB_spatial_;
                    }

                    // if (AA_params.reduce_nucleon_energies)
                    // {
                    //     double envelope_maximum = std::get<linear_interpolator>(env_func).value_at(newpair.getcr_sqrt_s());
                    //     pqcd::generate_bin_NN_coll
                    //     (
                    //         newpair, 
                    //         sigma_jet, 
                    //         AA_params.Tpp(newpair.getcr_bsquared()), 
                    //         kt0,
                    //         unirand, 
                    //         eng,
                    //         p_p_pdf,
                    //         dsigma_params,
                    //         power_law,
                    //         envelope_maximum
                    //     );
                    // }
                    // else
                    {
                        auto env_func_ = std::get<envelope_func>(env_func);

                        std::poisson_distribution<uint_fast8_t> dist(sigma_jet*AA_params.Tpp->at(new_pair));
                        uint_fast16_t nof_dijets = dist(*eng);
                        new_pair.dijets.reserve(nof_dijets);

                        const double sqrt_s = new_pair.getcr_sqrt_s();
                        for (uint_fast16_t i=0; i < nof_dijets; i++)
                        {
                            new_pair.dijets.push_back(pqcd::generate_2_to_2_scatt
                            (
                                sqrt_s,
                                kt0,
                                sqrt_s / 2.0,
                                unirand,
                                eng,
                                p_p_pdf,
                                dsigma_params,
                                power_law,
                                env_func_,
                                B->is_neutron,
                                A->is_neutron
                            ));
                        }

                        // pqcd::generate_bin_NN_coll
                        // (
                        //     newpair, 
                        //     sigma_jet, 
                        //     AA_params.Tpp(newpair.getcr_bsquared()), 
                        //     kt0,
                        //     unirand, 
                        //     eng,
                        //     p_p_pdf,
                        //     dsigma_params,
                        //     power_law,
                        //     env_func_,
                        //     B->is_neutron,
                        //     A->is_neutron
                        // );
                    }

                    // if (AA_params.reduce_nucleon_energies)
                    // {
                    //     newpair.reduce_energy_and_push_end_states_to_collider_frame();
                    // }
                    // else
                    {
                        new_pair.push_end_states_to_collider_frame();
                    }
                }
                new_pair.wound();
                binary_collisions.push_back(std::move(new_pair));
            }
            else
            {
                double sigma_jet;
                // if (AA_params.reduce_nucleon_energies) 
                // {
                //     array<double,3> args{*sum_tppa, *sum_tppb, pow(newpair.getcr_sqrt_s(), 2)};
                //     sigma_jet = std::get<InterpMultilinear<3, double> >(sigma_jets).interp(args.begin());
                // }
                // else
                {
                    array<double,2> args{*sum_tppa, *sum_tppb};
                    sigma_jet = std::get<InterpMultilinear<2, double> >(sigma_jets).interp(args.begin());
                    //pqcd::calculate_spatial_sigma_jet(p_p_pdf, p_n_pdf, &mand_s, &kt02, &jet_params, &sum_tppa, &sum_tppb, &tAA_0, &tBB_0);
                }
                
                // newpair.calculate_xsects(sigma_jet, AA_params.Tpp, newpair.getcr_bsquared());
                // auto ran = unirand(*eng);
                // switch (AA_params.normalize_to)
                // {
                // case B2_normalization_mode::total:
                //     ran *= newpair.getcr_max_tot_xsect();
                //     break;
                
                // case B2_normalization_mode::inelastic:
                //     ran *= newpair.getcr_max_inel_xsect();
                //     break;
                
                // case B2_normalization_mode::nothing:
                //     break;
                
                // default:
                //     ran *= newpair.getcr_max_inel_xsect();
                //     break;
                // }

                // if (ran > newpair.getcr_effective_inel_xsect())
                // {
                //     if (ran > newpair.getcr_effective_tot_xsect())
                //     {
                //         continue;
                //     }
                //     nof_softs++;
                //     continue;
                // }
                if (AA_params.calculate_end_state)
                {
                    
                    if (dsigma_params.snPDFs_linear)
                    {
                        //r=(1+cT)

                        const uint_fast16_t NA = dsigma_params.d_params.A, 
                                            NB = dsigma_params.d_params.B;
                                            
                        const auto spatial_cutoff = dsigma_params.d_params.spatial_cutoff;
                        //c=A*(R-1)/TAA(0)
                        const double scaA = static_cast<double>(NA) * *sum_tppa / tAA_0, 
                                    intA = 1.0 - scaA;
                        const std::function<double(double const&)> 
                            rA_spatial_ = [&,co=spatial_cutoff](double const &r)
                                {
                                    auto dummy = r*scaA + intA;
                                    // return dummy; // no cutoff
                                    // return (dummy < 0.0 ) ? 0.0 : dummy; // cutoff at 1+cT < 0
                                    return (dummy < co ) ? co : dummy; // cutoff at 1+cT < spatial_cutoff
                                    // return (dummy < 1/static_cast<double>(NA) ) ? 1/static_cast<double>(NA) : dummy; // cutoff at 1+cT < A
                                    
                                }; //r_s=1+c*sum(Tpp)

                        const double scaB = static_cast<double>(NB) * *sum_tppb / tBB_0, 
                                    intB = 1.0 - scaB;
                        const std::function<double(double const&)> 
                            rB_spatial_ = [&,co=spatial_cutoff](double const &r)
                                {
                                    auto dummy = r*scaB + intB;
                                    // return dummy; // no cutoff
                                    // return (dummy < 0.0 ) ? 0.0 : dummy; // cutoff at 1+cT < 0
                                    return (dummy < co ) ? co : dummy; // cutoff at 1+cT < spatial_cutoff
                                    // return (dummy < 1/static_cast<double>(NB) ) ? 1/static_cast<double>(NB) : dummy; // cutoff at 1+cT < B
                                };
                                
                        dsigma_params.d_params.rA_spatial = rA_spatial_;
                        dsigma_params.d_params.rB_spatial = rB_spatial_;
                    }
                    else
                    {
                        //r=exp(cT)

                        const std::function<double(double const&)> 
                            rA_spatial_ = [TAi=*sum_tppa,c_func=dsigma_params.c_A_func](double const &r)
                                {
                                    if (r>0.0)
                                    {
                                        double c = c_func.value_at(r);
                                        return std::exp(c * TAi);
                                    }
                                    else
                                    {
                                        return 0.0;
                                    }
                                };

                        const std::function<double(double const&)> 
                            rB_spatial_ = [TBi=*sum_tppb,c_func=dsigma_params.c_A_func](double const &r)
                                {
                                    if (r>0.0)
                                    {
                                        double c = c_func.value_at(r);
                                        return std::exp(c * TBi);
                                    }
                                    else
                                    {
                                        return 0.0;
                                    }
                                };

                        dsigma_params.d_params.rA_spatial = rA_spatial_;
                        dsigma_params.d_params.rB_spatial = rB_spatial_;
                    }

                    // if (AA_params.reduce_nucleon_energies)
                    // {
                    //     double envelope_maximum = std::get<linear_interpolator>(env_func).value_at(newpair.getcr_sqrt_s());
                    //     pqcd::generate_bin_NN_coll
                    //     (
                    //         newpair, 
                    //         sigma_jet, 
                    //         AA_params.Tpp(newpair.getcr_bsquared()), 
                    //         kt0,
                    //         unirand, 
                    //         eng,
                    //         p_p_pdf,
                    //         dsigma_params,
                    //         power_law,
                    //         envelope_maximum
                    //     );
                    // }
                    // else
                    {
                        auto env_func_ = std::get<envelope_func>(env_func);

                        std::poisson_distribution<uint_fast8_t> dist(sigma_jet*AA_params.Tpp->at(new_pair));
                        uint_fast8_t nof_dijets = dist(*eng);
                        new_pair.dijets.reserve(nof_dijets);

                        const double sqrt_s = new_pair.getcr_sqrt_s();

                        for (uint_fast8_t i=0; i < nof_dijets; i++)
                        {
                            auto new_dijet = pqcd::generate_2_to_2_scatt
                            (
                                sqrt_s,
                                kt0,
                                sqrt_s / 2.0,
                                unirand,
                                eng,
                                p_p_pdf,
                                dsigma_params,
                                power_law,
                                env_func_,
                                B->is_neutron,
                                A->is_neutron
                            );


                            if (is_sat_y_dep == 6)
                            {
                                auto kt2 = std::pow(new_dijet.kt,2);
                                auto dijet_area = p_p_pdf->alphasQ2(kt2)/kt2;

                                auto param = std::normal_distribution<double>::param_type{0., proton_width};
                                std::normal_distribution<double> normal_dist(0,0);
                                auto dx = normal_dist(*eng,param);
                                auto dy = normal_dist(*eng,param);

                                auto dijet_x = 0.5*(A->co.x + B->co.x + M_SQRT2*dx);
                                auto dijet_y = 0.5*(A->co.y + B->co.y + M_SQRT2*dy);
                                nucleon dummy_nucleon{coords{dijet_x, dijet_y, 0.0}, 0.0};
                                auto dijet_tppa = AA_params.Tpp->calculate_sum_tpp(dummy_nucleon, pro);
                                auto dijet_tppb = AA_params.Tpp->calculate_sum_tpp(dummy_nucleon, tar);

                                if (dijet_tppa * new_dijet.pro_pdf * dijet_area > M_factor)
                                {
                                    continue;
                                }
                                else if (dijet_tppb * new_dijet.pro_pdf * dijet_area > M_factor)
                                {
                                    continue;
                                }
                                else
                                {
                                    new_pair.dijets.push_back(new_dijet);
                                }
                            }
                            else
                            {
                                new_pair.dijets.push_back(new_dijet);
                            }
                        }

                        // pqcd::generate_bin_NN_coll
                        // (
                        //     newpair, 
                        //     sigma_jet, 
                        //     AA_params.Tpp(newpair.getcr_bsquared()), 
                        //     kt0,
                        //     unirand, 
                        //     eng,
                        //     p_p_pdf,
                        //     dsigma_params,
                        //     power_law,
                        //     env_func_,
                        //     B->is_neutron,
                        //     A->is_neutron
                        // );
                    }

                    // if (AA_params.reduce_nucleon_energies)
                    // {
                    //     newpair.reduce_energy_and_push_end_states_to_collider_frame();
                    // }
                    // else
                    {
                        new_pair.push_end_states_to_collider_frame();
                    }
                }
                new_pair.wound();
                binary_collisions.push_back(std::move(new_pair));
            }
        }

        if (verbose)
        {
            std::cout << "Bruteforced " << n_pairs << " pairs, got " << binary_collisions.size()+nof_softs << " collisions, of which softs "<< nof_softs<< " and hards "<< binary_collisions.size()<<" , double threshold broke " << mombroke << " times, skipped "<< skipped << " pairs that were too far apart" << std::endl;
        }
    }

    static auto collide_nuclei
    (
        std::vector<nucleon> &pro, 
        std::vector<nucleon> &tar, 
        std::vector<nn_coll> &binary_collisions, 
        variant_sigma_jet &sigma_jets,
        std::uniform_real_distribution<double> unirand, 
        std::shared_ptr<std::mt19937> eng, 
        const AA_collision_params &AA_params,
        pqcd::sigma_jet_params dsigma_params,
        const double &kt0,
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf,
        const double &power_law,
        variant_envelope_pars &env_func,
        const bool &verbose,
        const double &M_factor,
        const double &proton_width,
        const uint_fast16_t &is_sat_y_dep
    ) noexcept -> void
    {
        if (AA_params.spatial_pdfs)
        {
            collide_nuclei_snPDF
            (
                pro, 
                tar, 
                binary_collisions, 
                sigma_jets,
                unirand, 
                eng, 
                AA_params,
                dsigma_params,
                kt0,
                p_p_pdf,
                power_law,
                env_func,
                verbose,
                M_factor,
                proton_width,
                is_sat_y_dep
            );
        }
        else
        {
            collide_nuclei_no_snPDF
            (
                pro, 
                tar, 
                binary_collisions, 
                sigma_jets,
                unirand, 
                eng, 
                AA_params,
                dsigma_params,
                kt0,
                p_p_pdf,
                power_law,
                env_func,
                verbose
            );
        }
    }

    static auto read_nucleon_configs_from_file() noexcept
    {
        std::ifstream input("mc_initial.dat");

        auto pro = std::make_unique<std::vector<coords>>();
        auto tar = std::make_unique<std::vector<coords>>();

        if (input.is_open())
        {
            std::string line;

            for (uint_fast8_t i = 0; i < 8; i++)
            {
                std::getline(input, line); //Skip the 8 unimportant rows
            }

            while (true) //Loop for first nucleus
            {
                std::getline(input, line);

                if (line.empty()) //The nuclei are separated by empty line, followed by 4 comment lines
                {
                    for (uint_fast8_t i = 0; i < 4; i++)
                    {
                        std::getline(input, line);
                    }
                    break;
                }

                std::istringstream line_stream(line);
                std::vector<double> read_coords;
                double num = 0.0;

                while (line_stream >> num) 
                {
                    read_coords.push_back(num);
                }

                pro->push_back({read_coords[0], read_coords[1], read_coords[2]});
            }

            while (std::getline(input, line)) //Loop for other nucleus
            {
                std::istringstream line_stream(line);
                std::vector<double> read_coords;
                double num = 0.0;

                while (line_stream >> num) 
                {
                    read_coords.push_back(num);
                }
                
                tar->push_back({read_coords[0], read_coords[1], read_coords[2]});
            }
        }

        std::cout << "Read " << pro->size() << " projectile nuclei and " << tar->size() << " target nuclei!" << std::endl;
        return std::make_tuple(std::move(pro), std::move(tar));
    }

    static auto generate_nuclei
    (
        const nucleus_generator::nucleus_params &nuc_params,
        const double &sqrt_s,
        const double &impact_parameter,
        std::shared_ptr<std::mt19937> eng,
        const bool verbose,
        const bool read_nuclei_from_file=false
    )
    {
        std::vector<nucleon> pro, tar;
        if (read_nuclei_from_file)
        {
            if (verbose) std::cout<<"Reading nuclei..."<<std::flush;
            auto [pro_coords, tar_coords] = calcs::read_nucleon_configs_from_file();
            uint_fast16_t index = 0;
            for (auto &co : *pro_coords)
            {
                pro.emplace_back(co, sqrt_s / 2.0, index++);
            }
            index = 0;
            for (auto &co : *tar_coords)
            {
                tar.emplace_back(co, sqrt_s / 2.0, index++);
            }
            if (verbose) std::cout<<"Done!"<<std::endl;
        }
        else
        {
            if (verbose) std::cout<<"Generating nuclei..."<<std::flush;
            bool bugged;
            do //while (!bugged)
            {
                bugged = false;
                try
                {
                    pro = nucleus_generator::generate_nucleus
                        (
                            nuc_params,
                            false,
                            sqrt_s/2.0, 
                            -impact_parameter/2., 
                            eng
                        );
                    tar = nucleus_generator::generate_nucleus
                        (
                            nuc_params, 
                            true,
                            sqrt_s/2.0, 
                            impact_parameter/2., 
                            eng
                        );
                }
                catch(const std::exception& e)
                {
                    std::cout << e.what() << " in calcs, trying again"<<std::endl;
                    bugged = true;
                }
            } while (bugged);

            if (verbose) std::cout<<"Done!"<<std::endl;
        }

        return std::make_tuple(pro, tar);
    }

    // static auto ta_ws_folded
    // (
    //     const double &r
    // ) noexcept -> double
    // {
    //     if ((r > 0)&&(r < 2.99327))
    //     {
    //         return 2.1033+0.0387749*r-0.0641807*std::pow(r,1.5);
    //     }
    //     else if (r < 6.72925)
    //     {
    //         return 5.07336-2.16827*r+2.72582*std::pow(std::log(r),2);
    //     }
    //     else if (r < 9.12769)
    //     {
    //         return 13.4654+1.24837*r-11.2333*std::log(r);
    //     }
    //     else if (r < 15.6368)
    //     {
    //         return 0.00718685-0.000521874*r;
    //     }
    //     else 
    //     {
    //         return 0.0;
    //     }
    // }

private:

    static auto read_sigma_jets_spatial
    (
        const std::string &filename,
        const double &K_factor
    ) noexcept -> InterpMultilinear<2, double>
    {

        std::ifstream input(filename);

        std::array<uint_fast16_t,2> dim_Ns;
        std::vector<double> grid1, grid2;
        std::vector< std::vector<double>::iterator > grid_iter_list;
        std::vector<double> f_values;

        if (input.is_open())
        {
            std::string line;
            std::getline(input, line); //#1 Don't need anything from here

            std::getline(input, line); //#2
            std::istringstream line_stream(line);
            uint_fast64_t num_elements;
            line_stream.ignore(256,'=');
            line_stream >> num_elements;

            std::getline(input, line); //#3
            line_stream = std::istringstream(line);
            line_stream.ignore(256,'=');
            line_stream >> dim_Ns[0];
            std::getline(input, line); //#4
            line_stream = std::istringstream(line);
            line_stream.ignore(256,'=');
            line_stream >> dim_Ns[1];

            std::getline(input, line); //#5 empty
            std::getline(input, line); //#6 Don't need anything from here
            std::getline(input, line); //#7
            line_stream = std::istringstream(line);
            double num;
            while (line_stream >> num)
            {
                grid1.push_back(num);
            }
            std::getline(input, line); //#8 empty
            std::getline(input, line); //#9 Don't need anything from here
            std::getline(input, line); //#10
            line_stream = std::istringstream(line);
            while (line_stream >> num)
            {
                grid2.push_back(num);
            }
            grid_iter_list.push_back(grid1.begin());
            grid_iter_list.push_back(grid2.begin());
            
            std::getline(input, line); //#11 empty
            std::getline(input, line); //#12 Don't need anything from here
            f_values.reserve(num_elements);
            uint_fast64_t j=0;
            double sigma_jet;

            for (uint_fast64_t i=0; i<dim_Ns[0]; i++)
            {
                std::getline(input, line);
                line_stream = std::istringstream(line);
                while (line_stream >> sigma_jet)
                {
                    f_values[i*dim_Ns[1] + j] = K_factor*sigma_jet;
                    j++;
                }
                j=0;
            }

            return InterpMultilinear<2, double>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data() + num_elements);
        }

        std::cout<<"ERROR READING SIGMA_JETS"<<std::endl;
        
        return InterpMultilinear<2, double>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data());
    }

    // static auto calculate_spatial_sigma_jets_MC
    // (
    //     const double &tolerance, 
    //     std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
    //     /*std::shared_ptr<LHAPDF::GridPDF> p_n_pdf,*/ 
    //     const double &mand_s, 
    //     const double &kt02, 
    //     const pqcd::sigma_jet_params &jet_params, 
    //     const double &upper_sumTpp_limit, 
    //     const double &lower_sumTpp_limit
    // ) noexcept -> InterpMultilinear<3, double>
    // {
    //     const double marginal = 1.2; //20% more divisions than the tolerance gives us on the edges
    //     std::array<uint_fast16_t,3> dim_Ns{0}; //How many points to calculate in each dimension
    //     std::array<double,8> corners{0};
    //     const double tAA_0 = 29.5494;//30.5;
    //     const double tBB_0 = 29.5494;//30.5;

    //     auto sigma_jet_function = [=](const double sum_tppa, const double sum_tppb, const double mand_s_)
    //     {
    //         return pqcd::calculate_spatial_sigma_jet(p_p_pdf,/* p_n_pdf,*/ &mand_s_, &kt02, jet_params, sum_tppa, sum_tppb, tAA_0, tBB_0);
    //     };

    //     std::array<std::future<double>, 8> corner_futures{};
        
    //     //First calculation in a single thread, so the PDF gets fully initialized thread-safe
    //     corners[0]         =                                sigma_jet_function( upper_sumTpp_limit, upper_sumTpp_limit, mand_s);
    //     corner_futures[1]  = std::async(std::launch::async, sigma_jet_function, upper_sumTpp_limit, upper_sumTpp_limit, kt02  );
    //     corner_futures[2]  = std::async(std::launch::async, sigma_jet_function, upper_sumTpp_limit, lower_sumTpp_limit, mand_s);
    //     corner_futures[3]  = std::async(std::launch::async, sigma_jet_function, upper_sumTpp_limit, lower_sumTpp_limit, kt02  );
    //     corner_futures[4]  = std::async(std::launch::async, sigma_jet_function, lower_sumTpp_limit, upper_sumTpp_limit, mand_s);
    //     corner_futures[5]  = std::async(std::launch::async, sigma_jet_function, lower_sumTpp_limit, upper_sumTpp_limit, kt02  );
    //     corner_futures[6]  = std::async(std::launch::async, sigma_jet_function, lower_sumTpp_limit, lower_sumTpp_limit, mand_s);
    //     corner_futures[7]  = std::async(std::launch::async, sigma_jet_function, lower_sumTpp_limit, lower_sumTpp_limit, kt02  );

    //     for (uint_fast8_t i=1; i<8; i++)
    //     {
    //         corners[i] = corner_futures[i].get();
    //     }

    //     //Determine the grid spacings in all directions

    //     double max_corner = *std::max_element(corners.begin(), corners.end());

    //     //sum_Tpp_A
    //     std::array<double,4> differences
    //     {
    //         abs(corners[0]-corners[4]), 
    //         abs(corners[1]-corners[5]), 
    //         abs(corners[2]-corners[6]), 
    //         abs(corners[3]-corners[7])
    //     };
    //     double max_diff = *std::max_element(differences.begin(), differences.end());
    //     dim_Ns[0] = static_cast<uint_fast16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

    //     //sum_Tpp_B
    //     differences = std::array<double,4>
    //     ({
    //         abs(corners[0]-corners[2]), 
    //         abs(corners[1]-corners[3]), 
    //         abs(corners[4]-corners[6]), 
    //         abs(corners[5]-corners[7])
    //     });
    //     max_diff = *std::max_element(differences.begin(), differences.end());
    //     dim_Ns[1] = static_cast<uint_fast16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

    //     //mand_s
    //     differences = std::array<double,4>
    //     ({
    //         abs(corners[0]-corners[1]), 
    //         abs(corners[2]-corners[3]), 
    //         abs(corners[4]-corners[5]), 
    //         abs(corners[6]-corners[7])
    //     });
    //     max_diff = *std::max_element(differences.begin(), differences.end());
    //     dim_Ns[2] = static_cast<uint_fast16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

    //     for (auto & n : dim_Ns)
    //     {
    //         if (!std::isnormal(n) || n<2)
    //         {
    //             n=2;
    //         }
    //     }

    //     // construct the grid in each dimension. 
    //     // note that we will pass in a sequence of iterators pointing to the beginning of each grid
    //     std::vector<double> grid1 = helpers::linspace(lower_sumTpp_limit, upper_sumTpp_limit, dim_Ns[0]);
    //     std::vector<double> grid2 = helpers::linspace(lower_sumTpp_limit, upper_sumTpp_limit, dim_Ns[1]);
    //     std::vector<double> grid3 = helpers::linspace(kt02, mand_s, dim_Ns[2]);
    //     std::vector< std::vector<double>::iterator > grid_iter_list;
    //     grid_iter_list.push_back(grid1.begin());
    //     grid_iter_list.push_back(grid2.begin());
    //     grid_iter_list.push_back(grid3.begin());
    
    //     // total number of elements
    //     uint_fast64_t num_elements = dim_Ns[0] * dim_Ns[1] * dim_Ns[2]; 

    //     std::ofstream sigma_jet_grid_file;
    //     sigma_jet_grid_file.open("sigma_jet_grid_spagtial_MC.dat", std::ios::out);
    //     sigma_jet_grid_file << "%p_pdf=" << p_p_pdf->info().get_entry("SetIndex") << std::endl;
    //     sigma_jet_grid_file << "%num_elements=" << num_elements << std::endl;
    //     sigma_jet_grid_file << "%num_sum_T_pp_A=" << dim_Ns[0] << std::endl;
    //     sigma_jet_grid_file << "%num_sum_T_pp_B=" << dim_Ns[1] << std::endl;
    //     sigma_jet_grid_file << "%num_mand_s=" << dim_Ns[2] << std::endl;
    //     sigma_jet_grid_file << std::endl;
    //     sigma_jet_grid_file << "%sum_T_pp_A" << std::endl;
    //     for (auto g : grid1) sigma_jet_grid_file << g << ' ';
    //     sigma_jet_grid_file << std::endl << std::endl;
    //     sigma_jet_grid_file << "%sum_T_pp_B" << std::endl;
    //     for (auto g : grid2) sigma_jet_grid_file << g << ' ';
    //     sigma_jet_grid_file << std::endl << std::endl;
    //     sigma_jet_grid_file << "%mand_s" << std::endl;
    //     for (auto g : grid3) sigma_jet_grid_file << g << ' ';
    //     sigma_jet_grid_file << std::endl << std::endl;
    //     sigma_jet_grid_file << "%sigma_jet_values" << std::endl;

    //     // fill in the values of f(x) at the gridpoints. 
    //     // we will pass in a contiguous sequence, values are assumed to be laid out C-style
    //     std::vector<uint_fast64_t> c_style_indexes(num_elements);
    //     std::iota(c_style_indexes.begin(), c_style_indexes.end(), 0); //generates the list as {0,1,2,3,...}
    //     const uint_fast64_t rad1 = dim_Ns[1]*dim_Ns[2], 
    //                 rad2 = dim_Ns[2]; //These will help untangle the C-style index into coordinates
    //     // c_index = ii*rad1 + jj*rad2 + kk
    //     std::vector<double> f_values(num_elements);
    //     uint_fast64_t running_count{num_elements};
    //     std::mutex count_mutex;
        
    //     std::for_each
    //     (
    //         std::execution::par, 
    //         c_style_indexes.begin(), 
    //         c_style_indexes.end(), 
    //         [=, &f_values, &running_count, &count_mutex]
    //         (const uint_fast64_t index) 
    //         {
    //             uint_fast64_t kk = index % rad1 % rad2;
    //             uint_fast64_t i_dummy = index - kk; 
    //             uint_fast64_t jj = (i_dummy % rad1) / rad2;
    //             uint_fast64_t ii = (i_dummy - jj*rad2) / rad1;
    //             double dummy = pqcd::calculate_spatial_sigma_jet
    //                             (
    //                                 p_p_pdf, 
    //                                 /*p_n_pdf,*/ 
    //                                 &grid3[kk], 
    //                                 &kt02, 
    //                                 jet_params, 
    //                                 grid1[ii], 
    //                                 grid2[jj], 
    //                                 tAA_0, 
    //                                 tBB_0
    //                             );
    //             f_values[index] = dummy;
    //             {
    //                 const std::lock_guard<std::mutex> lock(count_mutex);
    //                 std::cout <<'\r'<<--running_count<<" left of "
    //                     <<num_elements<<" grid points to be calculated"<<std::flush;
    //             }
    //         }
    //     );
        
    //     for (uint_fast64_t i=0; i<dim_Ns[0]; i++)
    //     {
    //         for (uint_fast64_t j=0; j<dim_Ns[1]; j++)
    //         {
    //             for (uint_fast64_t k=0; k<dim_Ns[2]; k++)
    //             {
    //                 sigma_jet_grid_file << f_values[i*rad1 + j*rad2 + k] << ' ';
    //             }
    //             sigma_jet_grid_file << std::endl;
    //         }
    //         sigma_jet_grid_file << std::endl;
    //     }
    //     sigma_jet_grid_file.close();
    //     std::cout<<std::endl;

    //     return InterpMultilinear<3, double>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data() + num_elements);
    // }

    static auto calculate_spatial_sigma_jets
    (
        const std::string filename,
        const double &tolerance, 
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
        /*std::shared_ptr<LHAPDF::GridPDF> p_n_pdf,*/ 
        const double &mand_s, 
        const double &kt02, 
        const pqcd::sigma_jet_params &jet_params, 
        const double &upper_sumTpp_limit, 
        const double &lower_sumTpp_limit, 
        const double &ave_T_AA_0
    ) noexcept -> InterpMultilinear<2, double>
    {
        const double marginal = 1.2; //20% more divisions than the tolerance gives us on the edges
        std::array<uint_fast16_t,2> dim_Ns{0}; //How many points to calculate in each dimension
        std::array<double,4> corners{0};
        const double tAA_0 = ave_T_AA_0;//29.5494;//30.5
        const double tBB_0 = ave_T_AA_0;//29.5494;//30.5

        auto sigma_jet_function = [=](const double sum_tppa, const double sum_tppb)
        {
            return pqcd::calculate_spatial_sigma_jet(p_p_pdf,/* p_n_pdf,*/ &mand_s, &kt02, jet_params, sum_tppa, sum_tppb, tAA_0, tBB_0);
        };

        std::array<std::future<double>, 8> corner_futures{};
        
        //First calculation in a single thread, so the PDF gets fully initialized thread-safe
        corners[0]         =                                sigma_jet_function( upper_sumTpp_limit, upper_sumTpp_limit);
        corner_futures[1]  = std::async(std::launch::async, sigma_jet_function, upper_sumTpp_limit, lower_sumTpp_limit);
        corner_futures[2]  = std::async(std::launch::async, sigma_jet_function, lower_sumTpp_limit, upper_sumTpp_limit);
        corner_futures[3]  = std::async(std::launch::async, sigma_jet_function, lower_sumTpp_limit, lower_sumTpp_limit);

        for (uint_fast8_t i=1; i<4; i++)
        {
            corners[i] = corner_futures[i].get();
        }

        //Determine the grid spacings in all directions

        double max_corner = *std::max_element(corners.begin(), corners.end());

        //sum_Tpp_A
        std::array<double,2> differences
        {
            abs(corners[0]-corners[2]), 
            abs(corners[1]-corners[3])
        };
        double max_diff = *std::max_element(differences.begin(), differences.end());
        dim_Ns[0] = static_cast<uint_fast16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

        //sum_Tpp_B
        differences = std::array<double,2>
        ({
            abs(corners[0]-corners[1]), 
            abs(corners[2]-corners[3])
        });
        max_diff = *std::max_element(differences.begin(), differences.end());
        dim_Ns[1] = static_cast<uint_fast16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

        for (auto & n : dim_Ns)
        {
            if (!std::isnormal(n) || n<2)
            {
                n=2;
            }
        }

        // construct the grid in each dimension. 
        // note that we will pass in a sequence of iterators pointing to the beginning of each grid
        std::vector<double> grid1 = helpers::linspace(lower_sumTpp_limit, upper_sumTpp_limit, dim_Ns[0]);
        std::vector<double> grid2 = helpers::linspace(lower_sumTpp_limit, upper_sumTpp_limit, dim_Ns[1]);
        std::vector< std::vector<double>::iterator > grid_iter_list;
        grid_iter_list.push_back(grid1.begin());
        grid_iter_list.push_back(grid2.begin());
    
        // total number of elements
        uint_fast64_t num_elements = dim_Ns[0] * dim_Ns[1]; 

        std::ofstream sigma_jet_grid_file;
        sigma_jet_grid_file.open(filename, std::ios::out);
        sigma_jet_grid_file << "%mand_s=" << mand_s << " kt02=" << kt02 << " p_pdf=" << p_p_pdf->info().get_entry("SetIndex") << std::endl;
        sigma_jet_grid_file << "%num_elements=" << num_elements << std::endl;
        sigma_jet_grid_file << "%num_sum_T_pp_A=" << dim_Ns[0] << std::endl;
        sigma_jet_grid_file << "%num_sum_T_pp_B=" << dim_Ns[1] << std::endl;
        sigma_jet_grid_file << std::endl;
        sigma_jet_grid_file << "%sum_T_pp_A" << std::endl;
        for (auto g : grid1) sigma_jet_grid_file << g << ' ';
        sigma_jet_grid_file << std::endl << std::endl;
        sigma_jet_grid_file << "%sum_T_pp_B" << std::endl;
        for (auto g : grid2) sigma_jet_grid_file << g << ' ';
        sigma_jet_grid_file << std::endl << std::endl;
        sigma_jet_grid_file << "%sigma_jet_values" << std::endl;

        // fill in the values of f(x) at the gridpoints. 
        // we will pass in a contiguous sequence, values are assumed to be laid out C-style
        std::vector<uint_fast64_t> c_style_indexes(num_elements);
        std::iota(c_style_indexes.begin(), c_style_indexes.end(), 0); //generates the list as {0,1,2,3,...}
        const uint_fast64_t rad = dim_Ns[1]; //This will help untangle the C-style index into coordinates
        // c_index = ii*rad1 + jj
        std::vector<double> f_values(num_elements);
        uint_fast64_t running_count{num_elements};
        std::mutex count_mutex;
        
        #pragma omp parallel for
        for (auto it = c_style_indexes.begin(); it < c_style_indexes.end(); it++)
        {
            uint_fast64_t jj = *it % rad;
            uint_fast64_t ii = (*it - jj) / rad;
            double dummy = pqcd::calculate_spatial_sigma_jet
                            (
                                p_p_pdf, 
                                /*p_n_pdf,*/ 
                                &mand_s, 
                                &kt02, 
                                jet_params, 
                                grid1[ii], 
                                grid2[jj], 
                                tAA_0, 
                                tBB_0
                            );
            f_values[*it] = dummy;
            {
                const std::lock_guard<std::mutex> lock(count_mutex);
                std::cout <<'\r'<<--running_count<<" left of "
                    <<num_elements<<" grid points to be calculated"<<std::flush;
            }
        }
        
        double K_fac = jet_params.d_params.K_factor;

        for (uint_fast64_t i=0; i<dim_Ns[0]; i++)
        {
            for (uint_fast64_t j=0; j<dim_Ns[1]; j++)
            {
                sigma_jet_grid_file << f_values[i*rad + j]/K_fac<< ' ';
            }
            sigma_jet_grid_file << std::endl;
        }
        sigma_jet_grid_file.close();
        std::cout<<std::endl;

        return InterpMultilinear<2, double>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data() + num_elements);
    }
};

#endif // HIGH_LEVEL_CALCS_HPP
