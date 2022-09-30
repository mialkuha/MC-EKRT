//Copyright (c) 2022 Mikko Kuha

#ifndef HIGH_LEVEL_CALCS_HPP
#define HIGH_LEVEL_CALCS_HPP

#include <algorithm>
#include <execution>
#include <fstream>
#include <future>
#include <gsl/gsl_multimin.h>
#include <iostream>
#include <memory>
#include <sstream>
#include <tuple>

#include "generic_helpers.hpp"
#include "LHAPDF/GridPDF.h"
#include "linear_interpolator.hpp"
#include "linterp.h"
#include "nucleus_generator.hpp"
#include "pqcd.hpp"
#include "typedefs.hpp"

using variant_sigma_jet = std::variant<InterpMultilinear<3, xsectval>, InterpMultilinear<2, xsectval>, linear_interpolator, xsectval>;
using variant_envelope_max = std::variant<linear_interpolator, xsectval>;

class calcs
{
public:
    //class PrintThread: public std::ostringstream
    //{
    //public:
    //    PrintThread() = default;
//
    //    ~PrintThread()
    //    {
    //        std::lock_guard<std::mutex> guard(_mutexPrint);
    //        std::cout << this->str();
    //    }
//
    //private:
    //    static std::mutex _mutexPrint;
    //};
//
    //std::mutex PrintThread::_mutexPrint{};
    
    static auto find_max_dsigma
    (
        const momentum &kt, 
        const momentum &sqrt_s,
        std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
        pqcd::sigma_jet_params params
    ) noexcept
    {
        xsectval max_dsigma;
        xsectval error_est;

        if (params.d_params.npdfs_spatial)
        {
            //c=A*(R-1)/TAA(0)
            const double scaA = 208 * 0.01 / 29.5494,//30.5, 
                        intA = 1.0 - scaA;
            const std::function<double(double const&)> 
                rA_spatial_ = [&](double const &r)
                    {return r*scaA + intA;}; //r_s=1+c*sum(Tpp)
            params.d_params.rA_spatial = rA_spatial_;
            params.d_params.rB_spatial = rA_spatial_;
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

        uint8_t iter = 0;
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

            //if (status == GSL_SUCCESS)
            //{
            //    printf ("converged to maximum at\n");
            //    printf ("%5i %10.3e %10.3e f() = %10.3e size = %.3f\n",
            //            int(iter),
            //            gsl_vector_get (minimizer->x, 0),
            //            gsl_vector_get (minimizer->x, 1),
            //            -minimizer->fval, simplex_size);
            //}

        } while (status == GSL_CONTINUE && iter < 100);

        max_dsigma = static_cast<xsectval>(-minimizer->fval);
        error_est = static_cast<xsectval>(simplex_size);

        //  double y1 = gsl_vector_get(minimizer->x, 0);
        //  double y2 = gsl_vector_get(minimizer->x, 1);

        gsl_vector_free(ys);
        gsl_vector_free(init_step_size);
        gsl_multimin_fminimizer_free(minimizer);

        return std::make_tuple(std::move(max_dsigma), std::move(error_est));
    }

    static auto calculate_sigma_jets_for_MC
    (
        std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
        const momentum &mand_s, 
        const momentum &kt02,  
        pqcd::sigma_jet_params params,
        const uint16_t n_divisions = 100
    )
    {
        const std::vector<double> mand_ss{helpers::linspace(kt02, mand_s, n_divisions)};
        std::vector<std::tuple<double, double> > ret;
        std::mutex ret_mutex; 
        ret.push_back(std::make_tuple(kt02, 0.0));
        ret.push_back(std::make_tuple(mand_ss[1], pqcd::calculate_sigma_jet(p_pdf, &mand_ss[1], &kt02, params)));

        std::for_each
        (
            std::execution::par, 
            mand_ss.begin()+2, 
            mand_ss.end(), 
            [&](const double ss) 
            {
                auto sigma_jet = pqcd::calculate_sigma_jet(p_pdf, &ss, &kt02, params);
                const std::lock_guard<std::mutex> lock(ret_mutex);
                ret.push_back(std::make_tuple(ss, sigma_jet));
            }
        );
        
        //Sort the pairs by mand_s
        std::sort
        (
            ret.begin(), 
            ret.end(),
            [](std::tuple<double, double> a, std::tuple<double, double> b)
                {
                    auto [ a0, ph1 ] = a;
                    auto [ b0, ph2 ] = b;
                    return a0 < b0;
                }
        );

        std::vector<double> xs;
        std::vector<double> ys;

        for (auto p : ret)
        {
            auto [x, y] = p;
            xs.push_back(x);
            ys.push_back(y);
        }

        return std::make_tuple(xs, ys);
    }

    /*
    * Struct needed by find_max_dsigma
    */
    struct f_params
    {
        const momentum &kt;
        const momentum &sqrt_s;
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
        rapidity y1, y2;
        y1 = gsl_vector_get(v, 0);
        y2 = gsl_vector_get(v, 1);

        auto xsection = pqcd::diff_cross_section_2jet(sqrt_s, kt, y1, y2, p_pdf, sigma_params);
        xsectval total_xsection = 0;
        for (auto xsect : xsection)
        {
            total_xsection += xsect.sigma;
        }

        return -total_xsection;
    }

    static auto calculate_max_dsigmas_for_MC
    (
        const momentum &kt0,
        const momentum &sqrt_s,
        std::shared_ptr<LHAPDF::GridPDF> p_pdf,    
        pqcd::sigma_jet_params params,
        const uint16_t n_divisions = 100
    )
    {
        const std::vector<double> sqrt_ss{helpers::linspace(kt0, sqrt_s, n_divisions)};
        std::vector<std::tuple<double, std::tuple<double, double> > > ret;
        std::mutex ret_mutex; 
        ret.push_back(std::make_tuple(kt0, std::make_tuple(0.0, 0.0)));

        std::for_each
        (
            std::execution::par, 
            sqrt_ss.begin()+1, 
            sqrt_ss.end(), 
            [&](const double ss) 
            {
                auto max_dsigma = calcs::find_max_dsigma(kt0, ss, p_pdf, params);
                const std::lock_guard<std::mutex> lock(ret_mutex);
                ret.push_back(std::make_tuple(ss, max_dsigma));
            }
        );
        
        //Sort the pairs by sqrt_s
        std::sort
        (
            ret.begin(), 
            ret.end(),
            [](std::tuple<double, std::tuple<double, double> > a, std::tuple<double, std::tuple<double, double> > b)
                {
                    auto [ a0, ph1 ] = a;
                    auto [ b0, ph2 ] = b;
                    return a0 < b0;
                }
        );

        std::vector<double> xs;
        std::vector<std::tuple<double, double> > ys;

        for (auto p : ret)
        {
            auto [x, y] = p;
            xs.push_back(x);
            ys.push_back(y);
        }

        return std::make_tuple(xs, ys);
    }

    static auto prepare_sigma_jets
    (
        const bool use_npdfs,
        const bool spatial_pdfs,
        const bool reduce_nucleon_energies,
        const bool read_sigmajets_from_file,
        std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
        const momentum &mand_s,
        const momentum &sqrt_s,
        const momentum &kt02, 
        const momentum &kt0,
        pqcd::sigma_jet_params jet_params
    ) ->
    std::tuple
    <
        xsectval,
        double,
        variant_sigma_jet,
        std::optional<std::vector<double> >,
        variant_envelope_max,
        std::optional<std::vector<double> >
    >
    {
        double power_law = 0;
        xsectval dijet_norm = 0;

        if //proton PDFs
        (
            !use_npdfs && !spatial_pdfs
        )
        {
            if (!reduce_nucleon_energies) //Only one sigma_jet
            {
                power_law = 3.1;

                std::cout<<"Calculating sigma_jet..."<<std::flush;
                xsectval sigma_jet = pqcd::calculate_sigma_jet(p_pdf, &mand_s, &kt02, jet_params);
                std::cout<<"done!"<<std::endl;

                dijet_norm = sigma_jet;

                std::cout<<"Calculating envelope..."<<std::flush;
                auto [max_dsigma, err] = calcs::find_max_dsigma(kt0, sqrt_s, p_pdf, jet_params);
                xsectval envelope_maximum = (max_dsigma + fabs(err))*pow(kt0,power_law);
                std::cout<<"done!"<<std::endl;

                return std::make_tuple
                    (
                        dijet_norm,
                        power_law,
                        variant_sigma_jet(sigma_jet),
                        std::nullopt,
                        variant_envelope_max(envelope_maximum),
                        std::nullopt
                    );
            }
            else //sigma_jet depends on energy
            {
                power_law = 2.8;

                std::cout<<"Calculating sigma_jets..."<<std::flush;
                auto [mand_ss, sigmas] = calcs::calculate_sigma_jets_for_MC(p_pdf, mand_s, kt02, jet_params);
                auto sigma_jet = linear_interpolator(mand_ss, sigmas);
                std::cout<<"done!"<<std::endl;

                dijet_norm = sigmas.back();

                std::cout<<"Calculating envelope..."<<std::flush;
                auto [sqrt_ss, max_dsigmas] = calcs::calculate_max_dsigmas_for_MC(kt0, sqrt_s, p_pdf, jet_params);
                std::vector<double> envelopes;
                for (auto [ds, err] : max_dsigmas)
                {
                    envelopes.push_back((ds + fabs(err))*1.05*pow(kt0,power_law));
                }
                auto envelope_maximum = linear_interpolator(sqrt_ss, envelopes);
                std::cout<<"done!"<<std::endl;

                return std::make_tuple
                    (
                        dijet_norm,
                        power_law,
                        variant_sigma_jet(sigma_jet),
                        mand_ss,
                        variant_envelope_max(envelope_maximum),
                        sqrt_ss
                    );
            }
        }
        else if //nPDFS in use
        (
            !spatial_pdfs
        )
        {
            if (!reduce_nucleon_energies) //Only one sigma_jet
            {
                power_law = 2.5;

                std::cout<<"Calculating sigma_jet..."<<std::flush;
                xsectval sigma_jet = pqcd::calculate_sigma_jet(p_pdf, &mand_s, &kt02, jet_params);
                std::cout<<"done!"<<std::endl;

                dijet_norm = sigma_jet;

                std::cout<<"Calculating envelope..."<<std::flush;
                auto [max_dsigma, err] = calcs::find_max_dsigma(kt0, sqrt_s, p_pdf, jet_params);
                xsectval envelope_maximum = (max_dsigma + fabs(err))*pow(kt0,power_law);
                std::cout<<"done!"<<std::endl;

                return std::make_tuple
                    (
                        dijet_norm,
                        power_law,
                        variant_sigma_jet(sigma_jet),
                        std::nullopt,
                        variant_envelope_max(envelope_maximum),
                        std::nullopt
                    );
            }
            else //sigma_jet depends on energy
            {
                power_law = 2.2;

                std::cout<<"Calculating sigma_jets..."<<std::flush;
                auto [mand_ss, sigmas] = calcs::calculate_sigma_jets_for_MC(p_pdf, mand_s, kt02, jet_params);
                auto sigma_jet = linear_interpolator(mand_ss, sigmas);
                std::cout<<"done!"<<std::endl;

                dijet_norm = sigmas.back();

                std::cout<<"Calculating envelope..."<<std::flush;
                auto [sqrt_ss, max_dsigmas] = calcs::calculate_max_dsigmas_for_MC(kt0, sqrt_s, p_pdf, jet_params);
                std::vector<double> envelopes;
                for (auto [ds, err] : max_dsigmas)
                {
                    envelopes.push_back((ds + fabs(err))*1.05*pow(kt0,power_law));
                }
                auto envelope_maximum = linear_interpolator(sqrt_ss, envelopes);
                std::cout<<"done!"<<std::endl;

                return std::make_tuple
                    (
                        dijet_norm,
                        power_law,
                        variant_sigma_jet(sigma_jet),
                        mand_ss,
                        variant_envelope_max(envelope_maximum),
                        sqrt_ss
                    );
            }
        }
        else //spatial nPDFs
        {
            if (!reduce_nucleon_energies && read_sigmajets_from_file) //sigma_jet does not depend on energy
            {
                power_law = 2.0;

                std::cout<<"Reading spatial sigma_jets..."<<std::flush;
                variant_sigma_jet sigma_jet = calcs::read_sigma_jets_mf("sigma_jet_grid_mf.dat");
                std::cout<<"done!"<<std::endl;

                dijet_norm = 93.7604; //sigma_jet with ave nPDF

                std::cout<<"Calculating envelope..."<<std::flush;
                auto [max_dsigma, err] = calcs::find_max_dsigma(kt0, sqrt_s, p_pdf, jet_params);
                xsectval envelope_maximum = (max_dsigma + fabs(err))*1.05*pow(kt0,power_law);
                std::cout<<"done!"<<std::endl;

                return std::make_tuple
                    (
                        dijet_norm,
                        power_law,
                        std::move(sigma_jet),
                        std::nullopt,
                        variant_envelope_max(envelope_maximum),
                        std::nullopt
                    );
            }
            else if (!reduce_nucleon_energies)//sigma_jet does not depend on energy
            {
                power_law = 2.0;

                double tolerance=0.01,
                    upper_sumTpp_limit=0.61, 
                    lower_sumTpp_limit=0.01;

                std::cout<<"Calculating spatial sigma_jets..."<<std::endl;
                variant_sigma_jet sigma_jet 
                    = calcs::calculate_spatial_sigma_jets_mf
                        (
                            tolerance, 
                            p_pdf, 
                            mand_s, 
                            kt02, 
                            jet_params, 
                            upper_sumTpp_limit, 
                            lower_sumTpp_limit
                        );
                std::cout<<"done!"<<std::endl;

                dijet_norm = 93.7604; //sigma_jet with ave nPDF

                std::cout<<"Calculating envelope..."<<std::flush;
                auto [max_dsigma, err] = calcs::find_max_dsigma(kt0, sqrt_s, p_pdf, jet_params);
                xsectval envelope_maximum = (max_dsigma + fabs(err))*1.05*pow(kt0,power_law);
                std::cout<<"done!"<<std::endl;

                return std::make_tuple
                    (
                        dijet_norm,
                        power_law,
                        std::move(sigma_jet),
                        std::nullopt,
                        variant_envelope_max(envelope_maximum),
                        std::nullopt
                    );
            }
            else if (read_sigmajets_from_file)//sigma_jet depends on energy
            {
                power_law = 2.0;

                std::cout<<"Reading spatial sigma_jets..."<<std::flush;
                variant_sigma_jet sigma_jet = calcs::read_sigma_jets_mf_MC("sigma_jet_grid_mf_MC.dat");
                std::cout<<"done!"<<std::endl;

                dijet_norm = 93.7604; //sigma_jet with ave nPDF

                std::cout<<"Calculating envelope..."<<std::flush;
                auto [sqrt_ss, max_dsigmas] = calcs::calculate_max_dsigmas_for_MC(kt0, sqrt_s, p_pdf, jet_params);
                std::vector<double> envelopes;
                for (auto [ds, err] : max_dsigmas)
                {
                    envelopes.push_back((ds + fabs(err))*1.05*pow(kt0,power_law));
                }
                auto envelope_maximum = linear_interpolator(sqrt_ss, envelopes);
                std::cout<<"done!"<<std::endl;

                return std::make_tuple
                    (
                        dijet_norm,
                        power_law,
                        std::move(sigma_jet),
                        std::nullopt,
                        variant_envelope_max(envelope_maximum),
                        sqrt_ss
                    );
            }
            else //sigma_jet depends on energy
            {
                power_law = 2.0;

                double tolerance=0.02,
                    upper_sumTpp_limit=0.61, 
                    lower_sumTpp_limit=0.01;

                std::cout<<"Calculating spatial sigma_jets..."<<std::endl;
                variant_sigma_jet sigma_jet 
                    = calcs::calculate_spatial_sigma_jets_mf_MC
                        (
                            tolerance, 
                            p_pdf, 
                            mand_s, 
                            kt02, 
                            jet_params, 
                            upper_sumTpp_limit, 
                            lower_sumTpp_limit
                        );
                std::cout<<"done!"<<std::endl;

                dijet_norm = 93.7604; //sigma_jet with ave nPDF

                std::cout<<"Calculating envelope..."<<std::flush;
                auto [sqrt_ss, max_dsigmas] = calcs::calculate_max_dsigmas_for_MC(kt0, sqrt_s, p_pdf, jet_params);
                std::vector<double> envelopes;
                for (auto [ds, err] : max_dsigmas)
                {
                    envelopes.push_back((ds + fabs(err))*1.05*pow(kt0,power_law));
                }
                auto envelope_maximum = linear_interpolator(sqrt_ss, envelopes);
                std::cout<<"done!"<<std::endl;

                return std::make_tuple
                    (
                        dijet_norm,
                        power_law,
                        std::move(sigma_jet),
                        std::nullopt,
                        variant_envelope_max(envelope_maximum),
                        sqrt_ss
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
        const momentum &kt0,
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf,
        const double &power_law,
        variant_envelope_max &envelope_maximums,
        const bool &verbose
    ) noexcept -> void
    {
        uint n_pairs = 0, mombroke = 0, nof_softs = 0;
        
        std::vector<std::tuple<nucleon* const, nucleon* const> > binary_pairs;
        
        for (auto &A : pro)
        {
            for (auto &B : tar)
            {
                binary_pairs.push_back(std::make_pair(&A, &B));
            }
        }
        
        std::vector<uint64_t> pair_indexes(binary_pairs.size());
        std::iota(pair_indexes.begin(), pair_indexes.end(), 0);
        
        std::shuffle(pair_indexes.begin(), pair_indexes.end(), *eng);
        
        for (auto & ind : pair_indexes)
        {
            auto & [A, B] = binary_pairs[ind];
            
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
                const spatial d2 = AA_params.sigma_inel_for_glauber/(M_PI*10); // in fm^2
                const spatial dij2 = newpair.getcr_bsquared();
                
                if (dij2 > d2) //no collision
                {
                    continue;
                }
                //collision
                if (AA_params.calculate_end_state)
                {
                    xsectval sigma_jet;
                    xsectval envelope_maximum;
                    if (AA_params.reduce_nucleon_energies)
                    {
                        sigma_jet = std::get<linear_interpolator>(sigma_jets).value_at(pow(newpair.getcr_sqrt_s(), 2));
                        envelope_maximum = std::get<linear_interpolator>(envelope_maximums).value_at(newpair.getcr_sqrt_s());
                    }
                    else //Single sigma_jet
                    {
                        sigma_jet = std::get<xsectval>(sigma_jets);
                        envelope_maximum = std::get<xsectval>(envelope_maximums);
                    }
                    
                    pqcd::generate_bin_NN_coll
                    (
                        newpair, 
                        sigma_jet, 
                        AA_params.Tpp(newpair.getcr_bsquared()), 
                        kt0,
                        unirand, 
                        eng,
                        p_p_pdf,
                        dsigma_params,
                        power_law,
                        envelope_maximum
                    );
                    
                    if (AA_params.reduce_nucleon_energies)
                    {
                        newpair.reduce_energy_and_push_end_states_to_collider_frame();
                    }
                    else
                    {
                        newpair.push_end_states_to_collider_frame();
                    }
                }
                newpair.wound();
                binary_collisions.push_back(std::move(newpair));
            }
            else
            {
                xsectval sigma_jet;
                if (AA_params.reduce_nucleon_energies)
                {
                    sigma_jet = std::get<linear_interpolator>(sigma_jets).value_at(pow(newpair.getcr_sqrt_s(), 2));
                }
                else //Single sigma_jet
                {
                    sigma_jet = std::get<xsectval>(sigma_jets);
                }
                
                newpair.calculate_xsects(sigma_jet, AA_params.Tpp, newpair.getcr_bsquared(), AA_params.normalize_to);
                auto ran = unirand(*eng)*M_PI;
                
                if (ran > newpair.getcr_effective_inel_xsect())
                {
                    if (ran > newpair.getcr_effective_tot_xsect())
                    {
                        continue;
                    }
                    nof_softs++;
                    continue;
                }
                if (AA_params.calculate_end_state)
                {
                    xsectval envelope_maximum;
                    
                    if (AA_params.reduce_nucleon_energies)
                    {
                        envelope_maximum = std::get<linear_interpolator>(envelope_maximums).value_at(newpair.getcr_sqrt_s());
                    }
                    else
                    {
                        envelope_maximum = std::get<xsectval>(envelope_maximums);
                    }
                    
                    pqcd::generate_bin_NN_coll
                    (
                        newpair, 
                        sigma_jet, 
                        AA_params.Tpp(newpair.getcr_bsquared()), 
                        kt0,
                        unirand, 
                        eng,
                        p_p_pdf,
                        dsigma_params,
                        power_law,
                        envelope_maximum
                    );
                    
                    if (AA_params.reduce_nucleon_energies)
                    {
                        newpair.reduce_energy_and_push_end_states_to_collider_frame();
                    }
                    else
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
            std::cout << "Bruteforced " << n_pairs << " pairs, got " << binary_collisions.size()+nof_softs << " collisions, of which softs "<< nof_softs<< " and hards "<< binary_collisions.size()<<" , momentum threshold broke " << mombroke << " times" << std::endl;
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
        const momentum &kt0,
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf,
        const double &power_law,
        variant_envelope_max &envelope_maximums,
        const bool &verbose
    ) noexcept -> void
    {

        uint n_pairs = 0, mombroke = 0, skipped=0, nof_softs = 0;

        spatial tAA_0 = (AA_params.pA_scattering||AA_params.pp_scattering)? AA_params.Tpp(0.) : 29.5494;//30.5//calculate_tAB({0,0,0}, pro, pro, AA_params.Tpp);
        spatial tBB_0 = (AA_params.pp_scattering)? AA_params.Tpp(0.) : 29.5494;//30.5//calculate_tAB({0,0,0}, tar, tar, AA_params.Tpp);
        
        if (verbose)
        {
            std::cout << "T_AA(0)= " << tAA_0 << ", T_BB(0)= " << tBB_0 << std::endl;
        }

        std::vector<std::tuple
        <
            std::tuple<nucleon* const, const spatial* const>, 
            std::tuple<nucleon* const, const spatial* const> 
        > > binary_pairs;

        std::vector
        <
            std::tuple<nucleon* const, const spatial>
        > pro_spatial, tar_spatial;

        for (auto &A : pro)
        {
            pro_spatial.push_back(std::make_pair(&A, calculate_sum_tpp(A, pro, AA_params.Tpp)));
        }

        for (auto &B : tar)
        {
            tar_spatial.push_back(std::make_pair(&B, calculate_sum_tpp(B, tar, AA_params.Tpp)));
        }

        for (auto & [A, sum_A] : pro_spatial)
        {
            for (auto & [B, sum_B] : tar_spatial)
            {
                binary_pairs.push_back(std::make_pair(std::make_pair(A, &sum_A), std::make_pair(B, &sum_B)));
            }
        }
        
        std::vector<uint64_t> pair_indexes(binary_pairs.size());
        std::iota(pair_indexes.begin(), pair_indexes.end(), 0);

        std::shuffle(pair_indexes.begin(), pair_indexes.end(), *eng);

        for (auto & ind : pair_indexes)
        {
            auto & [A_pair, B_pair] = binary_pairs[ind];

            auto & [A, sum_tppa] = A_pair;
            auto & [B, sum_tppb] = B_pair;

            if (A->calculate_bsquared(*B) > 35.)
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
            xsectval sigma_jet;
            if (AA_params.reduce_nucleon_energies) 
            {
                array<double,3> args{*sum_tppa, *sum_tppb, pow(newpair.getcr_sqrt_s(), 2)};
                sigma_jet = std::get<InterpMultilinear<3, xsectval> >(sigma_jets).interp(args.begin());
            }
            else
            {
                array<spatial,2> args{*sum_tppa, *sum_tppb};
                sigma_jet = std::get<InterpMultilinear<2, xsectval> >(sigma_jets).interp(args.begin());
                //pqcd::calculate_spatial_sigma_jet_mf(p_p_pdf, p_n_pdf, &mand_s, &kt02, &jet_params, &sum_tppa, &sum_tppb, &tAA_0, &tBB_0);
            }
                
            newpair.calculate_xsects(sigma_jet, AA_params.Tpp, newpair.getcr_bsquared(), AA_params.normalize_to);

            if (verbose)
            {
                std::cout << "<T_pp>_i= " << *sum_tppa << ", <T_pp>_j= " << *sum_tppb << ", sigma_jet= " 
                        << sigma_jet <<  ", sigma_inel_eff= " << newpair.getcr_effective_inel_xsect() 
                        << ", sigma_tot_eff= " << newpair.getcr_effective_tot_xsect() << std::endl;
            }

            auto ran = unirand(*eng)*M_PI;
                
            if (ran > newpair.getcr_effective_inel_xsect())
            {
                if (ran > newpair.getcr_effective_tot_xsect())
                {
                    continue;
                }
                nof_softs++;
                continue;
            }
            if (AA_params.calculate_end_state)
            {
                xsectval envelope_maximum;

                if (AA_params.reduce_nucleon_energies)
                {
                    envelope_maximum = std::get<linear_interpolator>(envelope_maximums).value_at(newpair.getcr_sqrt_s());
                }
                else
                {
                    envelope_maximum = std::get<xsectval>(envelope_maximums);
                }

                const int NA = dsigma_params.d_params.A, 
                        NB = dsigma_params.d_params.B;
                //c=A*(R-1)/TAA(0)
                const double scaA = NA * *sum_tppa / tAA_0, 
                            intA = 1.0 - scaA;
                const std::function<double(double const&)> 
                    rA_spatial_ = [&](double const &r)
                        {return r*scaA + intA;}; //r_s=1+c*sum(Tpp)

                const double scaB = NB * *sum_tppb / tBB_0, 
                            intB = 1.0 - scaB;
                const std::function<double(double const&)> 
                    rB_spatial_ = [&](double const &r)
                        {return r*scaB + intB;};

                dsigma_params.d_params.rA_spatial = rA_spatial_;
                dsigma_params.d_params.rB_spatial = rB_spatial_;

                pqcd::generate_bin_NN_coll
                (
                    newpair, 
                    sigma_jet, 
                    AA_params.Tpp(newpair.getcr_bsquared()), 
                    kt0,
                    unirand, 
                    eng,
                    p_p_pdf,
                    dsigma_params,
                    power_law,
                    envelope_maximum
                );

                if (AA_params.reduce_nucleon_energies)
                {
                    newpair.reduce_energy_and_push_end_states_to_collider_frame();
                }
                else
                {
                    newpair.push_end_states_to_collider_frame();
                }
            }
            newpair.wound();
            binary_collisions.push_back(std::move(newpair));
        }

        if (verbose)
        {
            std::cout << "Bruteforced " << n_pairs << " pairs, got " << binary_collisions.size()+nof_softs << " collisions, of which softs "<< nof_softs<< " and hards "<< binary_collisions.size()<<" , momentum threshold broke " << mombroke << " times, skipped "<< skipped << " pairs that were too far apart" << std::endl;
        }
    }


    static auto collide_nuclei
    (
        const bool is_snPDF,
        std::vector<nucleon> &pro, 
        std::vector<nucleon> &tar, 
        std::vector<nn_coll> &binary_collisions, 
        variant_sigma_jet &sigma_jets,
        std::uniform_real_distribution<double> unirand, 
        std::shared_ptr<std::mt19937> eng, 
        const AA_collision_params &AA_params,
        pqcd::sigma_jet_params dsigma_params,
        const momentum &kt0,
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf,
        const double &power_law,
        variant_envelope_max &envelope_maximums,
        const bool &verbose
    ) noexcept -> void
    {
        if (is_snPDF)
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
                envelope_maximums,
                verbose
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
                envelope_maximums,
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

            for (uint8_t i = 0; i < 8; i++)
            {
                std::getline(input, line); //Skip the 8 unimportant rows
            }

            while (true) //Loop for first nucleus
            {
                std::getline(input, line);

                if (line.empty()) //The nuclei are separated by empty line, followed by 4 comment lines
                {
                    for (uint8_t i = 0; i < 4; i++)
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
        const momentum &sqrt_s,
        const spatial &impact_parameter,
        std::shared_ptr<std::mt19937> eng,
        std::shared_ptr<ars> radial_sampler,
        const bool &read_nuclei_from_file,
        const bool &verbose
    )
    {
        std::vector<nucleon> pro, tar;
        if (read_nuclei_from_file)
        {
            if (verbose) std::cout<<"Reading nuclei..."<<std::flush;
            auto [pro_coords, tar_coords] = calcs::read_nucleon_configs_from_file();
            for (auto &co : *pro_coords)
            {
                pro.emplace_back(co, sqrt_s / 2.0);
            }
            for (auto &co : *tar_coords)
            {
                tar.emplace_back(co, sqrt_s / 2.0);
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
                            eng, 
                            radial_sampler
                        );
                    tar = nucleus_generator::generate_nucleus
                        (
                            nuc_params, 
                            true,
                            sqrt_s/2.0, 
                            impact_parameter/2., 
                            eng, 
                            radial_sampler
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

        //if (verbose) std::cout<<"Shuffling nuclei..."<<std::flush;
        //std::shuffle(pro.begin(), pro.end(), *eng);
        //std::shuffle(tar.begin(), tar.end(), *eng);
        //if (verbose) std::cout<<"Done!"<<std::endl;

        return std::make_tuple(pro, tar);
    }

    
private:

    static auto calculate_sigma_1jet_analytical
    (
        const momentum &mand_s, 
        const std::vector<double> &kt_bin_walls,
        const std::vector<double> &y_bin_walls,
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
        pqcd::sigma_jet_params params, 
        const std::string name_postfix
    ) noexcept -> void
    {
        std::string dataname = "sigma1jet_analytical"+name_postfix+".dat";
        std::ofstream dat_file;
        dat_file.open(dataname);

        const auto n_kt_bins = kt_bin_walls.size() - 1;
        const auto n_y_bins = y_bin_walls.size() - 1;

        dat_file << "///y bin walls:  ";
        for (auto y : y_bin_walls)
        {
            dat_file<<y<<' ';
        }
        dat_file << std::endl;

        for (uint8_t i = 0; i < n_kt_bins; i++)
        {
            dat_file << kt_bin_walls[i] << ' ' << kt_bin_walls[i+1] << ' ';
            const auto kt_bin_size = kt_bin_walls[i+1] - kt_bin_walls[i];

            for (uint8_t j = 0; j < n_y_bins; j++)
            {
                const auto y_bin_size = y_bin_walls[j+1] - y_bin_walls[j];
                const auto bin = std::make_tuple(kt_bin_walls[i], kt_bin_walls[i+1], y_bin_walls[j], y_bin_walls[j+1]);
                
                const auto diff_xsect = pqcd::calculate_sigma_1jet_binned
                    (
                        p_p_pdf, 
                        &mand_s,
                        &bin,
                        params
                    )/(kt_bin_size*y_bin_size);
                dat_file << diff_xsect << ' ';
            }
            dat_file << std::endl;
        }
        std::cout<<"printed to "<<dataname<<std::endl;
    }

    static auto calculate_sigma_jet_analytical
    (
        const momentum &mand_s, 
        const std::vector<double> &kt_bin_walls,
        const std::vector<double> &y_bin_walls,
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
        pqcd::sigma_jet_params params
    ) noexcept -> void
    {
        std::string dataname = "sigmajet_analytical.dat";
        std::ofstream dat_file;
        dat_file.open(dataname);

        const auto n_kt_bins = kt_bin_walls.size() - 1;
        const auto n_y_bins = y_bin_walls.size() - 1;

        dat_file << "///y bin walls:  ";
        for (auto y : y_bin_walls)
        {
            dat_file<<y<<' ';
        }
        dat_file << std::endl;

        for (uint8_t i = 0; i < n_kt_bins; i++)
        {
            dat_file << kt_bin_walls[i] << ' ' << kt_bin_walls[i+1] << ' ';
            const auto kt_bin_size = kt_bin_walls[i+1] - kt_bin_walls[i];

            for (uint8_t j = 0; j < n_y_bins; j++)
            {
                const auto y_bin_size = y_bin_walls[j+1] - y_bin_walls[j];
                const auto bin = std::make_tuple(kt_bin_walls[i], kt_bin_walls[i+1], y_bin_walls[j], y_bin_walls[j+1]);
                
                dat_file << 
                    pqcd::calculate_sigma_jet_binned
                    (
                        p_p_pdf, 
                        &mand_s,
                        &bin,
                        params
                    )/(kt_bin_size*y_bin_size) << ' ';
            }
            dat_file << std::endl;
        }
        std::cout<<"printed to "<<dataname<<std::endl;
    }

    static auto calculate_sigma_dijet_analytical
    (
        const momentum &mand_s, 
        const std::vector<double> &pt_bin_walls,
        const std::vector<double> &eta_bin_walls,
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
        pqcd::sigma_jet_params params, 
        const std::string name_postfix
    ) noexcept -> void
    {
        std::string dataname = "sigmadijet_analytical"+name_postfix+".dat";
        std::ofstream dat_file;
        dat_file.open(dataname);

        const auto n_pt_bins = pt_bin_walls.size() - 1;
        const auto n_eta_bins = eta_bin_walls.size() - 1;

        dat_file << "///eta bin walls:  ";
        for (auto eta : eta_bin_walls)
        {
            dat_file<<eta<<' ';
        }
        dat_file << std::endl;

        for (uint8_t i = 0; i < n_pt_bins; i++)
        {
            dat_file << pt_bin_walls[i] << ' ' << pt_bin_walls[i+1] << ' ';
            const auto pt_bin_size = pt_bin_walls[i+1] - pt_bin_walls[i];

            for (uint8_t j = 0; j < n_eta_bins; j++)
            {
                const auto eta_bin_size = eta_bin_walls[j+1] - eta_bin_walls[j];
                const auto bin = std::make_tuple(pt_bin_walls[i], pt_bin_walls[i+1], eta_bin_walls[j], eta_bin_walls[j+1]);
                
                dat_file << 
                    pqcd::calculate_sigma_dijet_binned
                    (
                        p_p_pdf, 
                        &mand_s,
                        &bin,
                        params
                    )/(pt_bin_size*eta_bin_size) << ' ';
            }
            dat_file << std::endl;
        }
        std::cout<<"printed to "<<dataname<<std::endl;
    }

    static auto read_sigma_jets_mf_MC
    (
        const std::string &filename
    ) noexcept -> InterpMultilinear<3, xsectval>
    {

        std::ifstream input(filename);

        std::array<uint16_t,3> dim_Ns;
        std::vector<spatial> grid1, grid2, grid3;
        std::vector< std::vector<spatial>::iterator > grid_iter_list;
        std::vector<xsectval> f_values;

        if (input.is_open())
        {
            std::string line;
            std::getline(input, line); //#1 Don't need anything from here

            std::getline(input, line); //#2
            std::istringstream line_stream(line);
            uint64_t num_elements;
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
            std::getline(input, line); //#5
            line_stream = std::istringstream(line);
            line_stream.ignore(256,'=');
            line_stream >> dim_Ns[2];

            std::getline(input, line); //#6 empty
            std::getline(input, line); //#7 Don't need anything from here
            std::getline(input, line); //#8
            line_stream = std::istringstream(line);
            spatial num;
            while (line_stream >> num)
            {
                grid1.push_back(num);
            }
            std::getline(input, line); //#9 empty
            std::getline(input, line); //#10 Don't need anything from here
            std::getline(input, line); //#11
            line_stream = std::istringstream(line);
            while (line_stream >> num)
            {
                grid2.push_back(num);
            }
            std::getline(input, line); //#12 empty
            std::getline(input, line); //#13 Don't need anything from here
            std::getline(input, line); //#14
            line_stream = std::istringstream(line);
            while (line_stream >> num)
            {
                grid3.push_back(num);
            }
            grid_iter_list.push_back(grid1.begin());
            grid_iter_list.push_back(grid2.begin());
            grid_iter_list.push_back(grid3.begin());
            
            std::getline(input, line); //#15 empty
            std::getline(input, line); //#16 Don't need anything from here
            f_values.reserve(num_elements);
            uint64_t k=0;
            xsectval sigma_jet;

            for (uint64_t i=0; i<dim_Ns[0]; i++)
            {
                for (uint64_t j=0; j<dim_Ns[1]; j++)
                {
                    std::getline(input, line);
                    line_stream = std::istringstream(line);
                    while (line_stream >> sigma_jet)
                    {
                        f_values[i*dim_Ns[1]*dim_Ns[2] + j*dim_Ns[2] + k] = sigma_jet;
                        k++;
                    }
                    k=0;
                }
                std::getline(input, line); //empty
            }

            return InterpMultilinear<3, xsectval>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data() + num_elements);
        }

        std::cout<<"ERROR READING SIGMA_JETS"<<std::endl;
        
        return InterpMultilinear<3, xsectval>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data());
    }

    static auto read_sigma_jets_mf
    (
        const std::string &filename
    ) noexcept -> InterpMultilinear<2, xsectval>
    {

        std::ifstream input(filename);

        std::array<uint16_t,2> dim_Ns;
        std::vector<spatial> grid1, grid2;
        std::vector< std::vector<spatial>::iterator > grid_iter_list;
        std::vector<xsectval> f_values;

        if (input.is_open())
        {
            std::string line;
            std::getline(input, line); //#1 Don't need anything from here

            std::getline(input, line); //#2
            std::istringstream line_stream(line);
            uint64_t num_elements;
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
            spatial num;
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
            uint64_t j=0;
            xsectval sigma_jet;

            for (uint64_t i=0; i<dim_Ns[0]; i++)
            {
                std::getline(input, line);
                line_stream = std::istringstream(line);
                while (line_stream >> sigma_jet)
                {
                    f_values[i*dim_Ns[1] + j] = sigma_jet;
                    j++;
                }
                j=0;
            }

            return InterpMultilinear<2, xsectval>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data() + num_elements);
        }

        std::cout<<"ERROR READING SIGMA_JETS"<<std::endl;
        
        return InterpMultilinear<2, xsectval>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data());
    }

    static auto read_sigma_jets_full
    (
        const std::string &filename
    ) noexcept -> InterpMultilinear<5, xsectval>
    {

        std::ifstream input(filename);

        std::array<uint16_t,5> dim_Ns;
        std::vector<spatial> grid1, grid2, grid3, grid4, grid5;
        std::vector< std::vector<spatial>::iterator > grid_iter_list;
        std::vector<xsectval> f_values;

        if (input.is_open())
        {
            std::string line;
            std::getline(input, line); //#1 Don't need anything from here

            std::getline(input, line); //#2
            std::istringstream line_stream(line);
            uint64_t num_elements;
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
            std::getline(input, line); //#5
            line_stream = std::istringstream(line);
            line_stream.ignore(256,'=');
            line_stream >> dim_Ns[2];
            std::getline(input, line); //#6
            line_stream = std::istringstream(line);
            line_stream.ignore(256,'=');
            line_stream >> dim_Ns[3];
            std::getline(input, line); //#7
            line_stream = std::istringstream(line);
            line_stream.ignore(256,'=');
            line_stream >> dim_Ns[4];

            std::getline(input, line); //#8 empty
            std::getline(input, line); //#9 Don't need anything from here
            std::getline(input, line); //#10
            line_stream = std::istringstream(line);
            spatial num;
            while (line_stream >> num)
            {
                grid1.push_back(num);
            }
            std::getline(input, line); //#11 empty
            std::getline(input, line); //#12 Don't need anything from here
            std::getline(input, line); //#13
            line_stream = std::istringstream(line);
            while (line_stream >> num)
            {
                grid2.push_back(num);
            }
            std::getline(input, line); //#14 empty
            std::getline(input, line); //#15 Don't need anything from here
            std::getline(input, line); //#16
            line_stream = std::istringstream(line);
            while (line_stream >> num)
            {
                grid3.push_back(num);
            }
            std::getline(input, line); //#17 empty
            std::getline(input, line); //#18 Don't need anything from here
            std::getline(input, line); //#19
            line_stream = std::istringstream(line);
            while (line_stream >> num)
            {
                grid4.push_back(num);
            }
            std::getline(input, line); //#20 empty
            std::getline(input, line); //#21 Don't need anything from here
            std::getline(input, line); //#22
            line_stream = std::istringstream(line);
            while (line_stream >> num)
            {
                grid5.push_back(num);
            }
            grid_iter_list.push_back(grid1.begin());
            grid_iter_list.push_back(grid2.begin());
            grid_iter_list.push_back(grid3.begin());
            grid_iter_list.push_back(grid4.begin());
            grid_iter_list.push_back(grid5.begin());
            
            std::getline(input, line); //#23 empty
            std::getline(input, line); //#24 Don't need anything from here
            f_values.reserve(num_elements);
            uint16_t m=0;
            xsectval sigma_jet;

            for (uint16_t i=0; i<dim_Ns[0]; i++)
            {
                for (uint16_t j=0; j<dim_Ns[1]; j++)
                {
                    for (uint16_t k=0; k<dim_Ns[2]; k++)
                    {
                        for (uint16_t l=0; l<dim_Ns[3]; l++)
                        {
                            std::getline(input, line);
                            line_stream = std::istringstream(line);
                            while (line_stream >> sigma_jet)
                            {
                                f_values[ i*dim_Ns[1]*dim_Ns[2]*dim_Ns[3]*dim_Ns[4] 
                                        + j*dim_Ns[2]*dim_Ns[3]*dim_Ns[4] 
                                        + k*dim_Ns[3]*dim_Ns[4] 
                                        + l*dim_Ns[4] + m                            ] = sigma_jet;
                                m++;
                            }
                            m=0;
                        }
                        std::getline(input, line); //empty
                    }
                    std::getline(input, line); //empty
                }
                std::getline(input, line); //empty
            }

            return InterpMultilinear<5, xsectval>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data() + num_elements);
        }

        std::cout<<"ERROR READING SIGMA_JETS"<<std::endl;
        
        return InterpMultilinear<5, xsectval>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data());
    }

    static auto calculate_spatial_sigma_jets_mf_MC
    (
        const double &tolerance, 
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
        /*std::shared_ptr<LHAPDF::GridPDF> p_n_pdf,*/ 
        const momentum &mand_s, 
        const momentum &kt02, 
        const pqcd::sigma_jet_params &jet_params, 
        const double &upper_sumTpp_limit, 
        const double &lower_sumTpp_limit
    ) noexcept -> InterpMultilinear<3, xsectval>
    {
        const double marginal = 1.2; //20% more divisions than the tolerance gives us on the edges
        std::array<uint16_t,3> dim_Ns{0}; //How many points to calculate in each dimension
        std::array<xsectval,8> corners{0};
        const spatial tAA_0 = 29.5494;//30.5;
        const spatial tBB_0 = 29.5494;//30.5;

        auto sigma_jet_function = [=](const spatial sum_tppa, const spatial sum_tppb, const momentum mand_s_)
        {
            return pqcd::calculate_spatial_sigma_jet_mf(p_p_pdf,/* p_n_pdf,*/ &mand_s_, &kt02, jet_params, sum_tppa, sum_tppb, tAA_0, tBB_0);
        };

        std::array<std::future<xsectval>, 8> corner_futures{};
        
        //First calculation in a single thread, so the PDF gets fully initialized thread-safe
        corners[0]         =                                sigma_jet_function( upper_sumTpp_limit, upper_sumTpp_limit, mand_s);
        corner_futures[1]  = std::async(std::launch::async, sigma_jet_function, upper_sumTpp_limit, upper_sumTpp_limit, kt02  );
        corner_futures[2]  = std::async(std::launch::async, sigma_jet_function, upper_sumTpp_limit, lower_sumTpp_limit, mand_s);
        corner_futures[3]  = std::async(std::launch::async, sigma_jet_function, upper_sumTpp_limit, lower_sumTpp_limit, kt02  );
        corner_futures[4]  = std::async(std::launch::async, sigma_jet_function, lower_sumTpp_limit, upper_sumTpp_limit, mand_s);
        corner_futures[5]  = std::async(std::launch::async, sigma_jet_function, lower_sumTpp_limit, upper_sumTpp_limit, kt02  );
        corner_futures[6]  = std::async(std::launch::async, sigma_jet_function, lower_sumTpp_limit, lower_sumTpp_limit, mand_s);
        corner_futures[7]  = std::async(std::launch::async, sigma_jet_function, lower_sumTpp_limit, lower_sumTpp_limit, kt02  );

        for (uint8_t i=1; i<8; i++)
        {
            corners[i] = corner_futures[i].get();
        }

        //Determine the grid spacings in all directions

        xsectval max_corner = *std::max_element(corners.begin(), corners.end());

        //sum_Tpp_A
        std::array<xsectval,4> differences
        {
            abs(corners[0]-corners[4]), 
            abs(corners[1]-corners[5]), 
            abs(corners[2]-corners[6]), 
            abs(corners[3]-corners[7])
        };
        xsectval max_diff = *std::max_element(differences.begin(), differences.end());
        dim_Ns[0] = static_cast<uint16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

        //sum_Tpp_B
        differences = std::array<xsectval,4>
        ({
            abs(corners[0]-corners[2]), 
            abs(corners[1]-corners[3]), 
            abs(corners[4]-corners[6]), 
            abs(corners[5]-corners[7])
        });
        max_diff = *std::max_element(differences.begin(), differences.end());
        dim_Ns[1] = static_cast<uint16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

        //mand_s
        differences = std::array<xsectval,4>
        ({
            abs(corners[0]-corners[1]), 
            abs(corners[2]-corners[3]), 
            abs(corners[4]-corners[5]), 
            abs(corners[6]-corners[7])
        });
        max_diff = *std::max_element(differences.begin(), differences.end());
        dim_Ns[2] = static_cast<uint16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

        for (auto & n : dim_Ns)
        {
            if (!std::isnormal(n) || n<2)
            {
                n=2;
            }
        }

        // construct the grid in each dimension. 
        // note that we will pass in a sequence of iterators pointing to the beginning of each grid
        std::vector<spatial> grid1 = helpers::linspace(lower_sumTpp_limit, upper_sumTpp_limit, dim_Ns[0]);
        std::vector<spatial> grid2 = helpers::linspace(lower_sumTpp_limit, upper_sumTpp_limit, dim_Ns[1]);
        std::vector<spatial> grid3 = helpers::linspace(kt02, mand_s, dim_Ns[2]);
        std::vector< std::vector<spatial>::iterator > grid_iter_list;
        grid_iter_list.push_back(grid1.begin());
        grid_iter_list.push_back(grid2.begin());
        grid_iter_list.push_back(grid3.begin());
    
        // total number of elements
        uint64_t num_elements = dim_Ns[0] * dim_Ns[1] * dim_Ns[2]; 

        std::ofstream sigma_jet_grid_file;
        sigma_jet_grid_file.open("sigma_jet_grid_mf_MC.dat", std::ios::out);
        sigma_jet_grid_file << "%p_pdf=" << p_p_pdf->info().get_entry("SetIndex") << std::endl;
        sigma_jet_grid_file << "%num_elements=" << num_elements << std::endl;
        sigma_jet_grid_file << "%num_sum_T_pp_A=" << dim_Ns[0] << std::endl;
        sigma_jet_grid_file << "%num_sum_T_pp_B=" << dim_Ns[1] << std::endl;
        sigma_jet_grid_file << "%num_mand_s=" << dim_Ns[2] << std::endl;
        sigma_jet_grid_file << std::endl;
        sigma_jet_grid_file << "%sum_T_pp_A" << std::endl;
        for (auto g : grid1) sigma_jet_grid_file << g << ' ';
        sigma_jet_grid_file << std::endl << std::endl;
        sigma_jet_grid_file << "%sum_T_pp_B" << std::endl;
        for (auto g : grid2) sigma_jet_grid_file << g << ' ';
        sigma_jet_grid_file << std::endl << std::endl;
        sigma_jet_grid_file << "%mand_s" << std::endl;
        for (auto g : grid3) sigma_jet_grid_file << g << ' ';
        sigma_jet_grid_file << std::endl << std::endl;
        sigma_jet_grid_file << "%sigma_jet_values" << std::endl;

        // fill in the values of f(x) at the gridpoints. 
        // we will pass in a contiguous sequence, values are assumed to be laid out C-style
        std::vector<uint64_t> c_style_indexes(num_elements);
        std::iota(c_style_indexes.begin(), c_style_indexes.end(), 0); //generates the list as {0,1,2,3,...}
        const uint64_t rad1 = dim_Ns[1]*dim_Ns[2], 
                    rad2 = dim_Ns[2]; //These will help untangle the C-style index into coordinates
        // c_index = ii*rad1 + jj*rad2 + kk
        std::vector<xsectval> f_values(num_elements);
        std::atomic<uint64_t> running_count{num_elements};
        
        std::for_each
        (
            std::execution::par, 
            c_style_indexes.begin(), 
            c_style_indexes.end(), 
            [=, &f_values, &running_count](const uint64_t index) 
            {
                uint64_t kk = index % rad1 % rad2;
                uint64_t i_dummy = index - kk; 
                uint64_t jj = (i_dummy % rad1) / rad2;
                uint64_t ii = (i_dummy - jj*rad2) / rad1;
                xsectval dummy = pqcd::calculate_spatial_sigma_jet_mf
                                (
                                    p_p_pdf, 
                                    /*p_n_pdf,*/ 
                                    &grid3[kk], 
                                    &kt02, 
                                    jet_params, 
                                    grid1[ii], 
                                    grid2[jj], 
                                    tAA_0, 
                                    tBB_0
                                );
                f_values[index] = dummy;
                //PrintThread{} <<'\r'<<--running_count<<" left of "
                //            <<num_elements<<" grid points to be calculated";
            }
        );
        
        for (uint64_t i=0; i<dim_Ns[0]; i++)
        {
            for (uint64_t j=0; j<dim_Ns[1]; j++)
            {
                for (uint64_t k=0; k<dim_Ns[2]; k++)
                {
                    sigma_jet_grid_file << f_values[i*rad1 + j*rad2 + k] << ' ';
                }
                sigma_jet_grid_file << std::endl;
            }
            sigma_jet_grid_file << std::endl;
        }
        sigma_jet_grid_file.close();
        std::cout<<std::endl;

        return InterpMultilinear<3, xsectval>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data() + num_elements);
    }

    static auto calculate_spatial_sigma_jets_mf
    (
        const double &tolerance, 
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
        /*std::shared_ptr<LHAPDF::GridPDF> p_n_pdf,*/ 
        const momentum &mand_s, 
        const momentum &kt02, 
        const pqcd::sigma_jet_params &jet_params, 
        const double &upper_sumTpp_limit, 
        const double &lower_sumTpp_limit
    ) noexcept -> InterpMultilinear<2, xsectval>
    {
        const double marginal = 1.2; //20% more divisions than the tolerance gives us on the edges
        std::array<uint16_t,2> dim_Ns{0}; //How many points to calculate in each dimension
        std::array<xsectval,4> corners{0};
        const spatial tAA_0 = 29.5494;//30.5
        const spatial tBB_0 = 29.5494;//30.5

        auto sigma_jet_function = [=](const spatial sum_tppa, const spatial sum_tppb)
        {
            return pqcd::calculate_spatial_sigma_jet_mf(p_p_pdf,/* p_n_pdf,*/ &mand_s, &kt02, jet_params, sum_tppa, sum_tppb, tAA_0, tBB_0);
        };

        std::array<std::future<xsectval>, 8> corner_futures{};
        
        //First calculation in a single thread, so the PDF gets fully initialized thread-safe
        corners[0]         =                                sigma_jet_function( upper_sumTpp_limit, upper_sumTpp_limit);
        corner_futures[1]  = std::async(std::launch::async, sigma_jet_function, upper_sumTpp_limit, lower_sumTpp_limit);
        corner_futures[2]  = std::async(std::launch::async, sigma_jet_function, lower_sumTpp_limit, upper_sumTpp_limit);
        corner_futures[3]  = std::async(std::launch::async, sigma_jet_function, lower_sumTpp_limit, lower_sumTpp_limit);

        for (uint8_t i=1; i<4; i++)
        {
            corners[i] = corner_futures[i].get();
        }

        //Determine the grid spacings in all directions

        xsectval max_corner = *std::max_element(corners.begin(), corners.end());

        //sum_Tpp_A
        std::array<xsectval,2> differences
        {
            abs(corners[0]-corners[2]), 
            abs(corners[1]-corners[3])
        };
        xsectval max_diff = *std::max_element(differences.begin(), differences.end());
        dim_Ns[0] = static_cast<uint16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

        //sum_Tpp_B
        differences = std::array<xsectval,2>
        ({
            abs(corners[0]-corners[1]), 
            abs(corners[2]-corners[3])
        });
        max_diff = *std::max_element(differences.begin(), differences.end());
        dim_Ns[1] = static_cast<uint16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

        for (auto & n : dim_Ns)
        {
            if (!std::isnormal(n) || n<2)
            {
                n=2;
            }
        }

        // construct the grid in each dimension. 
        // note that we will pass in a sequence of iterators pointing to the beginning of each grid
        std::vector<spatial> grid1 = helpers::linspace(lower_sumTpp_limit, upper_sumTpp_limit, dim_Ns[0]);
        std::vector<spatial> grid2 = helpers::linspace(lower_sumTpp_limit, upper_sumTpp_limit, dim_Ns[1]);
        std::vector< std::vector<spatial>::iterator > grid_iter_list;
        grid_iter_list.push_back(grid1.begin());
        grid_iter_list.push_back(grid2.begin());
    
        // total number of elements
        uint64_t num_elements = dim_Ns[0] * dim_Ns[1]; 

        std::ofstream sigma_jet_grid_file;
        sigma_jet_grid_file.open("sigma_jet_grid_mf.dat", std::ios::out);
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
        std::vector<uint64_t> c_style_indexes(num_elements);
        std::iota(c_style_indexes.begin(), c_style_indexes.end(), 0); //generates the list as {0,1,2,3,...}
        const uint64_t rad = dim_Ns[1]; //This will help untangle the C-style index into coordinates
        // c_index = ii*rad1 + jj
        std::vector<xsectval> f_values(num_elements);
        std::atomic<uint64_t> running_count{num_elements};
        
        std::for_each
        (
            std::execution::par, 
            c_style_indexes.begin(), 
            c_style_indexes.end(), 
            [=, &f_values, &running_count](const uint64_t index) 
            {
                uint64_t jj = index % rad;
                uint64_t ii = (index - jj) / rad;
                xsectval dummy = pqcd::calculate_spatial_sigma_jet_mf
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
                f_values[index] = dummy;
                //PrintThread{} <<'\r'<<--running_count<<" left of "
                //            <<num_elements<<" grid points to be calculated";
            }
        );
        
        for (uint64_t i=0; i<dim_Ns[0]; i++)
        {
            for (uint64_t j=0; j<dim_Ns[1]; j++)
            {
                sigma_jet_grid_file << f_values[i*rad + j] << ' ';
            }
            sigma_jet_grid_file << std::endl;
        }
        sigma_jet_grid_file.close();
        std::cout<<std::endl;

        return InterpMultilinear<2, xsectval>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data() + num_elements);
    }

    static auto calculate_spatial_sigma_jets_full
    (
        const double &tolerance, 
        std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
        std::shared_ptr<LHAPDF::GridPDF> p_n_pdf, 
        const momentum &mand_s, 
        const momentum &kt02, 
        const pqcd::sigma_jet_params &jet_params, 
        const std::array<const double, 4> &lower_limits, 
        const std::array<const double, 4> &upper_limits
    ) noexcept -> InterpMultilinear<5, xsectval>
    {
        const double marginal = 1.2; //20% more divisions than the tolerance gives us on the edges
        std::array<uint16_t,5> dim_Ns{0}; //How many points to calculate in each dimension
        std::array<xsectval,32> corners{0};

        auto sigma_jet_function = [=](const double T_sums_1, const double T_sums_2, const double T_sums_3, const spatial tAA_0, const spatial tBB_0)
        {
            const std::array<const double, 3> T_sums{T_sums_1, T_sums_2, T_sums_3};
            return pqcd::calculate_spatial_sigma_jet_full(p_p_pdf, p_n_pdf, &mand_s, &kt02, jet_params, T_sums, tAA_0, tBB_0);
        };

        std::array<std::future<xsectval>, 32> corner_futures{};
        
        //First calculation in a single thread, so the PDF gets fully initialized thread-safe
        corners[0]         = sigma_jet_function(                                upper_limits[0], upper_limits[1], upper_limits[2], upper_limits[3], upper_limits[3]);
        std::cout<<corners[0]<<std::endl;
        //corner_futures[1]  = std::async(std::launch::async, sigma_jet_function, upper_limits[0], upper_limits[1], upper_limits[2], upper_limits[3], lower_limits[3]);
        //corner_futures[2]  = std::async(std::launch::async, sigma_jet_function, upper_limits[0], upper_limits[1], upper_limits[2], lower_limits[3], upper_limits[3]);
        //corner_futures[3]  = std::async(std::launch::async, sigma_jet_function, upper_limits[0], upper_limits[1], upper_limits[2], lower_limits[3], lower_limits[3]);
        corner_futures[4]  = std::async(std::launch::async, sigma_jet_function, upper_limits[0], upper_limits[1], lower_limits[2], upper_limits[3], upper_limits[3]);
        //corner_futures[5]  = std::async(std::launch::async, sigma_jet_function, upper_limits[0], upper_limits[1], lower_limits[2], upper_limits[3], lower_limits[3]);
        //corner_futures[6]  = std::async(std::launch::async, sigma_jet_function, upper_limits[0], upper_limits[1], lower_limits[2], lower_limits[3], upper_limits[3]);
        //corner_futures[7]  = std::async(std::launch::async, sigma_jet_function, upper_limits[0], upper_limits[1], lower_limits[2], lower_limits[3], lower_limits[3]);
        corner_futures[8]  = std::async(std::launch::async, sigma_jet_function, upper_limits[0], lower_limits[1], upper_limits[2], upper_limits[3], upper_limits[3]);
        //corner_futures[9]  = std::async(std::launch::async, sigma_jet_function, upper_limits[0], lower_limits[1], upper_limits[2], upper_limits[3], lower_limits[3]);
        //corner_futures[10] = std::async(std::launch::async, sigma_jet_function, upper_limits[0], lower_limits[1], upper_limits[2], lower_limits[3], upper_limits[3]);
        //corner_futures[11] = std::async(std::launch::async, sigma_jet_function, upper_limits[0], lower_limits[1], upper_limits[2], lower_limits[3], lower_limits[3]);
        corner_futures[12] = std::async(std::launch::async, sigma_jet_function, upper_limits[0], lower_limits[1], lower_limits[2], upper_limits[3], upper_limits[3]);
        //corner_futures[13] = std::async(std::launch::async, sigma_jet_function, upper_limits[0], lower_limits[1], lower_limits[2], upper_limits[3], lower_limits[3]);
        //corner_futures[14] = std::async(std::launch::async, sigma_jet_function, upper_limits[0], lower_limits[1], lower_limits[2], lower_limits[3], upper_limits[3]);
        //corner_futures[15] = std::async(std::launch::async, sigma_jet_function, upper_limits[0], lower_limits[1], lower_limits[2], lower_limits[3], lower_limits[3]);
        corner_futures[16] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], upper_limits[1], upper_limits[2], upper_limits[3], upper_limits[3]);
        //corner_futures[17] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], upper_limits[1], upper_limits[2], upper_limits[3], lower_limits[3]);
        //corner_futures[18] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], upper_limits[1], upper_limits[2], lower_limits[3], upper_limits[3]);
        //corner_futures[19] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], upper_limits[1], upper_limits[2], lower_limits[3], lower_limits[3]);
        corner_futures[20] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], upper_limits[1], lower_limits[2], upper_limits[3], upper_limits[3]);
        //corner_futures[21] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], upper_limits[1], lower_limits[2], upper_limits[3], lower_limits[3]);
        //corner_futures[22] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], upper_limits[1], lower_limits[2], lower_limits[3], upper_limits[3]);
        //corner_futures[23] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], upper_limits[1], lower_limits[2], lower_limits[3], lower_limits[3]);
        corner_futures[24] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], lower_limits[1], upper_limits[2], upper_limits[3], upper_limits[3]);
        //corner_futures[25] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], lower_limits[1], upper_limits[2], upper_limits[3], lower_limits[3]);
        //corner_futures[26] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], lower_limits[1], upper_limits[2], lower_limits[3], upper_limits[3]);
        //corner_futures[27] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], lower_limits[1], upper_limits[2], lower_limits[3], lower_limits[3]);
        corner_futures[28] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], lower_limits[1], lower_limits[2], upper_limits[3], upper_limits[3]);
        //corner_futures[29] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], lower_limits[1], lower_limits[2], upper_limits[3], lower_limits[3]);
        //corner_futures[30] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], lower_limits[1], lower_limits[2], lower_limits[3], upper_limits[3]);
        //corner_futures[31] = std::async(std::launch::async, sigma_jet_function, lower_limits[0], lower_limits[1], lower_limits[2], lower_limits[3], lower_limits[3]);
    uint8_t iiii=0;
        for (uint8_t i=4; i<32; i+=4)
        {
            corners[i] = corner_futures[i].get();
            iiii=1;
            corners[i+iiii++] = corners[i];
            corners[i+iiii++] = corners[i];
            corners[i+iiii++] = corners[i];
        }

        //Determine the grid spacings in all directions

        xsectval max_corner = *std::max_element(corners.begin(), corners.end());

        //sum_1
        std::array<xsectval,16> differences;
        for (uint8_t i=0; i<16; i++)
        {
            differences[i] = abs(corners[i]-corners[i+16]);
        }
        xsectval max_diff = *std::max_element(differences.begin(), differences.end());
        dim_Ns[0] = static_cast<uint16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

        //sum_2
        {
        uint8_t j=0;
        for (uint8_t i=0; i<8; i++)
        {
            differences[i] = abs(corners[j]-corners[j+8]);
            j++;
        }
        j=16;
        for (uint8_t i=8; i<16; i++)
        {
            differences[i] = abs(corners[j]-corners[j+8]);
            j++;
        }
        max_diff = *std::max_element(differences.begin(), differences.end());
        dim_Ns[1] = static_cast<uint16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

        //sum_3
        j=0;
        for (uint8_t i=0; i<16; i++)
        {
            differences[i++] = abs(corners[j]-corners[j+4]);
            differences[i++] = abs(corners[j+1]-corners[j+5]);
            differences[i++] = abs(corners[j+2]-corners[j+6]);
            differences[i] = abs(corners[j+3]-corners[j+7]);
            j+=8;
        }
        max_diff = *std::max_element(differences.begin(), differences.end());
        dim_Ns[2] = static_cast<uint16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));
        }
        //TAA(0)
        /*
        j=0;
        for (uint8_t i=0; i<16; i++)
        {
            differences[i++] = abs(corners[j]-corners[j+2]);
            differences[i] = abs(corners[j+1]-corners[j+3]);
            j+=4;
        }
        max_diff = *std::max_element(differences.begin(), differences.end());
        dim_Ns[3] = static_cast<uint16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));*/
        dim_Ns[3] = 1;
        
        //TBB(0)
        /*
        j=0;
        for (uint8_t i=0; i<16; i++)
        {
            differences[i] = abs(corners[j]-corners[j+1]);
            j+=2;
        }
        max_diff = *std::max_element(differences.begin(), differences.end());
        dim_Ns[4] = static_cast<uint16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));*/
        dim_Ns[4] = 1;

        for (auto & n : dim_Ns)
        {
            if (!std::isnormal(n) || n<1)
            {
                n=1;
            }
            std::cout<< n<< ' ';
        }
            std::cout<<std::endl;

        // construct the grid in each dimension. 
        // note that we will pass in a sequence of iterators pointing to the beginning of each grid
        std::vector<spatial> grid1 = helpers::linspace(lower_limits[0], upper_limits[0], dim_Ns[0]);
        std::vector<spatial> grid2 = helpers::linspace(lower_limits[1], upper_limits[1], dim_Ns[1]);
        std::vector<spatial> grid3 = helpers::linspace(lower_limits[2], upper_limits[2], dim_Ns[2]);
        std::vector<spatial> grid4 = /*linspace(*/{lower_limits[3]}/*, upper_limits[3], dim_Ns[3])*/;
        std::vector<spatial> grid5 = /*linspace(*/{lower_limits[3]}/*, upper_limits[3], dim_Ns[4])*/;
        std::vector< std::vector<spatial>::iterator > grid_iter_list;
        grid_iter_list.push_back(grid1.begin());
        grid_iter_list.push_back(grid2.begin());
        grid_iter_list.push_back(grid3.begin());
        grid_iter_list.push_back(grid4.begin());
        grid_iter_list.push_back(grid5.begin());
    
        // total number of elements
        uint64_t num_elements = dim_Ns[0] * dim_Ns[1] * dim_Ns[2] * dim_Ns[3] * dim_Ns[4]; 

        std::ofstream sigma_jet_grid_file;
        sigma_jet_grid_file.open("sigma_jet_full_grid.dat", std::ios::out);
        sigma_jet_grid_file << "%mand_s=" << mand_s << " kt02=" << kt02 << " p_pdf=" << p_p_pdf->info().get_entry("SetIndex") << std::endl;
        sigma_jet_grid_file << "%num_elements=" << num_elements << std::endl;
        sigma_jet_grid_file << "%num_sum_1=" << dim_Ns[0] << std::endl;
        sigma_jet_grid_file << "%num_sum_2=" << dim_Ns[1] << std::endl;
        sigma_jet_grid_file << "%num_sum_3=" << dim_Ns[2] << std::endl;
        sigma_jet_grid_file << "%num_T_AA(0)=" << dim_Ns[3] << std::endl;
        sigma_jet_grid_file << "%num_T_BB(0)=" << dim_Ns[4] << std::endl;
        sigma_jet_grid_file << std::endl;
        sigma_jet_grid_file << "%sum_1" << std::endl;
        for (auto g : grid1) sigma_jet_grid_file << g << ' ';
        sigma_jet_grid_file << std::endl << std::endl;
        sigma_jet_grid_file << "%sum_2" << std::endl;
        for (auto g : grid2) sigma_jet_grid_file << g << ' ';
        sigma_jet_grid_file << std::endl << std::endl;
        sigma_jet_grid_file << "%sum_3" << std::endl;
        for (auto g : grid3) sigma_jet_grid_file << g << ' ';
        sigma_jet_grid_file << std::endl << std::endl;
        sigma_jet_grid_file << "%T_AA(0)" << std::endl;
        for (auto g : grid4) sigma_jet_grid_file << g << ' ';
        sigma_jet_grid_file << std::endl << std::endl;
        sigma_jet_grid_file << "%T_BB(0)" << std::endl;
        for (auto g : grid5) sigma_jet_grid_file << g << ' ';
        sigma_jet_grid_file << std::endl << std::endl;
        sigma_jet_grid_file << "%sigma_jet_values" << std::endl;

        // fill in the values of f(x) at the gridpoints. 
        // we will pass in a contiguous sequence, values are assumed to be laid out C-style
        std::vector<uint64_t> c_style_indexes(num_elements);
        std::iota(c_style_indexes.begin(), c_style_indexes.end(), 0); //generates the list as {0,1,2,3,...}
        const uint64_t rad1 = dim_Ns[1]*dim_Ns[2]*dim_Ns[3]*dim_Ns[4], 
                    rad2 = dim_Ns[2]*dim_Ns[3]*dim_Ns[4], 
                    rad3 = dim_Ns[3]*dim_Ns[4], 
                    rad4 = dim_Ns[4]; //These will help untangle the C-style index into coordinates
        // c_index = ii*rad1 + jj*rad2 + kk*rad3 + ll*rad4 + mm
        std::vector<xsectval> f_values(num_elements);
        std::atomic<uint64_t> running_count{num_elements};
        
        std::for_each
        (
            std::execution::par, 
            c_style_indexes.begin(), 
            c_style_indexes.end(), 
            [=, &f_values, &running_count](const uint64_t index) 
            {
                uint64_t mm = index % rad1 % rad2 % rad3 % rad4;
                uint64_t i_dummy = (index - mm); 
                uint64_t ll = (i_dummy % rad1 % rad2 % rad3) / rad4;
                i_dummy = (i_dummy - ll*rad4); 
                uint64_t kk = (i_dummy % rad1 % rad2) / rad3;
                i_dummy = (i_dummy - kk*rad3); 
                uint64_t jj = (i_dummy % rad1) / rad2;
                uint64_t ii = (i_dummy - jj*rad2) / rad1;
                xsectval dummy = sigma_jet_function(grid1[ii], grid2[jj], grid3[kk], grid4[ll], grid5[mm]);
                //xsectval dummy = pqcd::calculate_spatial_sigma_jet_full(p_p_pdf, p_n_pdf, &mand_s, &kt02, &jet_params, &grid1[ii], &grid2[jj], &grid3[kk], &grid4[ll], &grid5[mm]);
                f_values[index] = dummy;
                //PrintThread{} <<'\r'<<--running_count<<" left of "<<num_elements<<" grid points to be calculated";
            }
        );


        for (uint64_t i=0; i<dim_Ns[0]; i++)
        {
            for (uint64_t j=0; j<dim_Ns[1]; j++)
            {
                for (uint64_t k=0; k<dim_Ns[2]; k++)
                {
                    for (uint64_t l=0; l<dim_Ns[3]; l++)
                    {
                        for (uint64_t m=0; m<dim_Ns[4]; m++)
                        {
                            sigma_jet_grid_file << f_values[i*rad1 + j*rad2 + k*rad3 + l*rad4 + m] << ' ';
                        }      
                        sigma_jet_grid_file << std::endl;
                    }
                    sigma_jet_grid_file << std::endl;
                }
                sigma_jet_grid_file << std::endl;
            }
            sigma_jet_grid_file << std::endl;
        }
        sigma_jet_grid_file.close();
        std::cout<<std::endl;

        return InterpMultilinear<5, xsectval>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data() + num_elements);
    }

    static auto calculate_tA
    (
        const coords &b, 
        const std::vector<nucleon> &nucleus,
        const std::function<spatial(const spatial&)> Tp
    ) noexcept -> spatial
    {
        spatial tA=0.0; //sum(T_p(b_ij + b))
        uint16_t A=static_cast<uint16_t>(nucleus.size());

        auto [bx, by, bz] = b;

        for (uint16_t i=0; i<A; i++)
        {
            //tAB1 += Tpp((pro.at(i).co - tar.at(j).co + b).magt2());
            auto [x1, y1, z1] = nucleus.at(i).co;
            tA += Tp(pow(x1-bx,2) + pow(y1-by,2));
        }
        return tA;
    }

    static auto calculate_tAB
    (
        const coords &b, 
        const std::vector<nucleon> &pro, 
        const std::vector<nucleon> &tar, 
        const std::function<spatial(const spatial&)> Tpp
    ) noexcept -> spatial
    {
        spatial tAB=0.0; //sum(T_pp(b_ij + b))
        uint16_t A=static_cast<uint16_t>(pro.size()), 
                B=static_cast<uint16_t>(tar.size());

        auto [bx, by, bz] = b;

        for (uint16_t i=0; i<A; i++)
        {
            for (uint16_t j=0; j<B; j++)
            {
                //tAB1 += Tpp((pro.at(i).co - tar.at(j).co + b).magt2());
                auto [x1, y1, z1] = pro.at(i).co;
                auto [x2, y2, z2] = tar.at(j).co;
                tAB += Tpp(pow(x1-x2+bx,2) + pow(y1-y2+by,2));
            }   
        }
        return tAB;
    }

    static auto calculate_sum_tpp
    (
        const nucleon &nuc, 
        const std::vector<nucleon> &nucleus, 
        const std::function<spatial(const spatial&)> Tpp
    ) noexcept -> spatial
    {
        spatial sum_tpp=0.0; //sum(T_pp(b_ii'))
        uint16_t A=static_cast<uint16_t>(nucleus.size());

        auto [x1, y1, z1] = nuc.co;

        for (uint16_t i=0; i<A; i++)
        {
            //sum_tpp += Tpp((nuc.co - nucleus.at(i).co).magt2());
            auto [x2, y2, z2] = nucleus.at(i).co;
            sum_tpp += Tpp(pow(x1-x2,2) + pow(y1-y2,2));
        }
        //std::cout<<sum_tpp1<<' '<<sum_tpp2<<std::endl;
        return sum_tpp;
    }

    static void calculate_and_save_nuclei_TAs_TAAs
    (
        const nucleus_generator::nucleus_params &nuc_params,
        std::shared_ptr<std::mt19937> eng,
        std::shared_ptr<ars> radial_sampler
    )
    {
        uint64_t num_nuclei = 10000;
        std::array<double,201> grid_xs;
        std::array<double,201> grid_ys;

        const spatial proton_width_2 = pow(0.573, 2);
        const std::function<spatial(const spatial&)> 
            Tp{[&proton_width_2](const spatial &bsquared)
            {
                return exp(-bsquared / (2 * proton_width_2)) / (20 * M_PI * proton_width_2); // 1/fm² = mb/fm² * 1/mb = 0.1 * 1/mb
            }}; 
        //const std::function<spatial(const spatial&)> 
        //    Tpp{[&proton_width_2](const spatial &bsquared)
        //    {
        //        return exp(-bsquared / (4 * proton_width_2)) / (40 * M_PI * proton_width_2); // 1/fm² = mb/fm² * 1/mb = 0.1 * 1/mb
        //    }};

        std::ofstream nuclei_file;
        std::ofstream TAs_file;
        //std::ofstream TAAs_file;
        std::mutex nuclei_file_mutex;
        std::mutex TAs_file_mutex;
        //std::mutex TAAs_file_mutex;
        std::mutex radial_sampler_mutex; 

        std::iota(grid_xs.begin(), grid_xs.end(), 0);
        std::iota(grid_ys.begin(), grid_ys.end(), 0);

        for (auto & x : grid_xs)
        {
            x = -10.0 + 0.1*x;
        }
        for (auto & y : grid_ys)
        {
            y = -10.0 + 0.1*y;
        }

        std::array<std::array<coords,201>,201> grid;
        for (uint16_t i=0; i<201; i++)
        {
            grid[i] = std::array<coords,201>();
            for (uint16_t j=0; j<201; j++)
            {
                grid[i][j] = coords({grid_xs[i], grid_ys[j], 0});
            }
        }
        nuclei_file.open("nuclei_shift.wl");
        TAs_file.open("TAs_shift.wl");
        //TAAs_file.open("TAAs_noshift.wl");

        nuclei_file<<"nucleiShift"<<num_nuclei<<" = {";
        TAs_file<<"TAsShift"<<num_nuclei<<" = {";
        //TAAs_file<<"TAAsNoShift"<<num_nuclei<<" = {";
        std::vector<uint64_t> indexes((num_nuclei/2)-1);
        std::iota(indexes.begin(), indexes.end(), 0); //generates the list as {0,1,2,3,...}

        calcs::generate_nuclei
        (
            nuc_params, 
            1000, 
            0, 
            eng, 
            radial_sampler, 
            false, 
            false
        );

        std::for_each
        (
            std::execution::par, 
            indexes.begin(), 
            indexes.end(), 
            [&](const uint64_t index) 
            {

                std::unique_lock<std::mutex> lock_rad(radial_sampler_mutex);
                auto [pro, tar] = calcs::generate_nuclei
                (
                    nuc_params, 
                    1000, 
                    0, 
                    eng, 
                    radial_sampler, 
                    false, 
                    false
                );
                lock_rad.~unique_lock();

                {
                    const std::lock_guard<std::mutex> lock(nuclei_file_mutex);
                    nuclei_file<<"{";
                    std::for_each(pro.begin(), pro.end()-1, [&nuclei_file](const nucleon &nuc) 
                    { 
                        nuclei_file<<"{"<<nuc.co.x<<","<<nuc.co.y<<","<<nuc.co.z<<"},\n    ";
                    });
                    nuclei_file<<"{"<<pro.back().co.x<<","<<pro.back().co.y<<","<<pro.back().co.z<<"}},\n    {";
                    std::for_each(tar.begin(), tar.end()-1, [&nuclei_file](const nucleon &nuc) 
                    { 
                        nuclei_file<<"{"<<nuc.co.x<<","<<nuc.co.y<<","<<nuc.co.z<<"},\n    ";
                    });
                    nuclei_file<<"{"<<tar.back().co.x<<","<<tar.back().co.y<<","<<tar.back().co.z<<"}},\n    ";
                }

                std::array<std::array<double,201>,201> grid_TAs;
                //std::array<std::array<double,201>,201> grid_TAAs;
                std::array<std::array<double,201>,201> grid_TBs;
                //std::array<std::array<double,201>,201> grid_TBBs;

                for (uint16_t i=0; i<201; i++)
                {
                    grid_TAs[i] = std::array<double,201>();
                    //grid_TAAs[i] = std::array<double,201>();
                    grid_TBs[i] = std::array<double,201>();
                    //grid_TBBs[i] = std::array<double,201>();
                    for (uint16_t j=0; j<201; j++)
                    {
                        grid_TAs[i][j] = calculate_tA(grid[i][j], pro, Tp);
                        //grid_TAAs[i][j] = calculate_tAB(grid[i][j], pro, pro, Tpp);
                        grid_TBs[i][j] = calculate_tA(grid[i][j], tar, Tp);
                        //grid_TBBs[i][j] = calculate_tAB(grid[i][j], tar, tar, Tpp);
                    }
                }

                {
                    const std::lock_guard<std::mutex> lock(TAs_file_mutex);
                    TAs_file<<"{";
                    std::for_each(grid_TAs.begin(), grid_TAs.end()-1, [&TAs_file](const std::array<double,201> &ys) 
                    { 
                        TAs_file<<"{";
                        std::for_each(ys.begin(), ys.end()-1, [&TAs_file](const double &val) 
                        { 
                            TAs_file<<val<<",";
                        });
                        TAs_file<<ys.back()<<"},\n    ";
                    });
                    TAs_file<<"{";
                    std::for_each(grid_TAs.back().begin(), grid_TAs.back().end()-1, [&TAs_file](const double &val) 
                    { 
                        TAs_file<<val<<",";
                    });
                    TAs_file<<grid_TAs.back().back()<<"}},\n    {";
                    std::for_each(grid_TBs.begin(), grid_TBs.end()-1, [&TAs_file](const std::array<double,201> &ys) 
                    { 
                        TAs_file<<"{";
                        std::for_each(ys.begin(), ys.end()-1, [&TAs_file](const double &val) 
                        { 
                            TAs_file<<val<<",";
                        });
                        TAs_file<<ys.back()<<"},\n    ";
                    });
                    TAs_file<<"{";
                    std::for_each(grid_TBs.back().begin(), grid_TBs.back().end()-1, [&TAs_file](const double &val) 
                    { 
                        TAs_file<<val<<",";
                    });
                    TAs_file<<grid_TBs.back().back()<<"}},\n    ";
                }

                //{
                //    const std::lock_guard<std::mutex> lock(TAAs_file_mutex);
                //    TAAs_file<<"{";
                //    std::for_each(grid_TAAs.begin(), grid_TAAs.end()-1, [&TAAs_file](const std::array<double,201> &ys) 
                //    { 
                //        TAAs_file<<"{";
                //        std::for_each(ys.begin(), ys.end()-1, [&TAAs_file](const double &val) 
                //        { 
                //            TAAs_file<<val<<",";
                //        });
                //        TAAs_file<<ys.back()<<"},\n    ";
                //    });
                //    TAAs_file<<"{";
                //    std::for_each(grid_TAAs.back().begin(), grid_TAAs.back().end()-1, [&TAAs_file](const double &val) 
                //    { 
                //        TAAs_file<<val<<",";
                //    });
                //    TAAs_file<<grid_TAAs.back().back()<<"}},\n    {";
                //    std::for_each(grid_TBBs.begin(), grid_TBBs.end()-1, [&TAAs_file](const std::array<double,201> &ys) 
                //    { 
                //        TAAs_file<<"{";
                //        std::for_each(ys.begin(), ys.end()-1, [&TAAs_file](const double &val) 
                //        { 
                //            TAAs_file<<val<<",";
                //        });
                //        TAAs_file<<ys.back()<<"},\n    ";
                //    });
                //    TAAs_file<<"{";
                //    std::for_each(grid_TBBs.back().begin(), grid_TBBs.back().end()-1, [&TAAs_file](const double &val) 
                //    { 
                //        TAAs_file<<val<<",";
                //    });
                //    TAAs_file<<grid_TBBs.back().back()<<"}},\n    ";
                //    std::cout<<index<<std::flush;
                //}
            }
        );

        auto [pro, tar] = generate_nuclei
        (
            nuc_params, 
            1000, 
            0, 
            eng, 
            radial_sampler, 
            false, 
            false
        );

        nuclei_file<<"{";
        std::for_each(pro.begin(), pro.end()-1, [&nuclei_file](const nucleon &nuc) 
        { 
            nuclei_file<<"{"<<nuc.co.x<<","<<nuc.co.y<<","<<nuc.co.z<<"},\n    ";
        });
        nuclei_file<<"{"<<pro.back().co.x<<","<<pro.back().co.y<<","<<pro.back().co.z<<"}},\n    {";
        std::for_each(tar.begin(), tar.end()-1, [&nuclei_file](const nucleon &nuc) 
        { 
            nuclei_file<<"{"<<nuc.co.x<<","<<nuc.co.y<<","<<nuc.co.z<<"},\n    ";
        });
        nuclei_file<<"{"<<tar.back().co.x<<","<<tar.back().co.y<<","<<tar.back().co.z<<"}}};\n\n";
        nuclei_file.close();

        std::array<std::array<double,201>,201> grid_TAs;
        //std::array<std::array<double,201>,201> grid_TAAs;
        std::array<std::array<double,201>,201> grid_TBs;
        //std::array<std::array<double,201>,201> grid_TBBs;

        for (uint16_t i=0; i<201; i++)
        {
            grid_TAs[i] = std::array<double,201>();
            //grid_TAAs[i] = std::array<double,201>();
            grid_TBs[i] = std::array<double,201>();
            //grid_TBBs[i] = std::array<double,201>();
            for (uint16_t j=0; j<201; j++)
            {
                grid_TAs[i][j] = calculate_tA(grid[i][j], pro, Tp);
                //grid_TAAs[i][j] = calculate_tAB(grid[i][j], pro, pro, Tpp);
                grid_TBs[i][j] = calculate_tA(grid[i][j], tar, Tp);
                //grid_TBBs[i][j] = calculate_tAB(grid[i][j], tar, tar, Tpp);
            }
        }

        TAs_file<<"{";
        std::for_each(grid_TAs.begin(), grid_TAs.end()-1, [&TAs_file](const std::array<double,201> &ys) 
        { 
            TAs_file<<"{";
            std::for_each(ys.begin(), ys.end()-1, [&TAs_file](const double &val) 
            { 
                TAs_file<<val<<",";
            });
            TAs_file<<ys.back()<<"},\n    ";
        });
        TAs_file<<"{";
        std::for_each(grid_TAs.back().begin(), grid_TAs.back().end()-1, [&TAs_file](const double &val) 
        { 
            TAs_file<<val<<",";
        });
        TAs_file<<grid_TAs.back().back()<<"}},\n    {";
        std::for_each(grid_TBs.begin(), grid_TBs.end()-1, [&TAs_file](const std::array<double,201> &ys) 
        { 
            TAs_file<<"{";
            std::for_each(ys.begin(), ys.end()-1, [&TAs_file](const double &val) 
            { 
                TAs_file<<val<<",";
            });
            TAs_file<<ys.back()<<"},\n    ";
        });
        TAs_file<<"{";
        std::for_each(grid_TBs.back().begin(), grid_TBs.back().end()-1, [&TAs_file](const double &val) 
        { 
            TAs_file<<val<<",";
        });
        TAs_file<<grid_TBs.back().back()<<"}}};\n\n";
        TAs_file.close();

        //TAAs_file<<"{";
        //std::for_each(grid_TAAs.begin(), grid_TAAs.end()-1, [&TAAs_file](const std::array<double,201> &ys) 
        //{ 
        //    TAAs_file<<"{";
        //    std::for_each(ys.begin(), ys.end()-1, [&TAAs_file](const double &val) 
        //    { 
        //        TAAs_file<<val<<",";
        //    });
        //    TAAs_file<<ys.back()<<"},\n    ";
        //});
        //TAAs_file<<"{";
        //std::for_each(grid_TAAs.back().begin(), grid_TAAs.back().end()-1, [&TAAs_file](const double &val) 
        //{ 
        //    TAAs_file<<val<<",";
        //});
        //TAAs_file<<grid_TAAs.back().back()<<"}},\n    {";
        //std::for_each(grid_TBBs.begin(), grid_TBBs.end()-1, [&TAAs_file](const std::array<double,201> &ys) 
        //{ 
        //    TAAs_file<<"{";
        //    std::for_each(ys.begin(), ys.end()-1, [&TAAs_file](const double &val) 
        //    { 
        //        TAAs_file<<val<<",";
        //    });
        //    TAAs_file<<ys.back()<<"},\n    ";
        //});
        //TAAs_file<<"{";
        //std::for_each(grid_TBBs.back().begin(), grid_TBBs.back().end()-1, [&TAAs_file](const double &val) 
        //{ 
        //    TAAs_file<<val<<",";
        //});
        //TAAs_file<<grid_TBBs.back().back()<<"}}};\n\n";
        //TAAs_file.close();
    }

    static void calculate_and_save_average_nuclei_TAs_TAAs
    (
        const nucleus_generator::nucleus_params &nuc_params,
        std::shared_ptr<std::mt19937> eng,
        std::shared_ptr<ars> radial_sampler
    )
    {
        uint64_t num_nuclei = 10000;
        std::array<double,201> grid_xs;
        std::array<double,201> grid_ys;

        const spatial proton_width_2 = pow(0.573, 2);
        const std::function<spatial(const spatial&)> 
            Tp{[&proton_width_2](const spatial &bsquared)
            {
                return exp(-bsquared / (2 * proton_width_2)) / (20 * M_PI * proton_width_2); // 1/fm² = mb/fm² * 1/mb = 0.1 * 1/mb
            }}; 

        std::ofstream TAs_file;
        std::mutex radial_sampler_mutex; 

        std::iota(grid_xs.begin(), grid_xs.end(), 0);
        std::iota(grid_ys.begin(), grid_ys.end(), 0);

        for (auto & x : grid_xs)
        {
            x = -10.0 + 0.1*x;
        }
        for (auto & y : grid_ys)
        {
            y = -10.0 + 0.1*y;
        }

        std::array<std::array<coords,201>,201> grid;
        for (uint16_t i=0; i<201; i++)
        {
            grid[i] = std::array<coords,201>();
            for (uint16_t j=0; j<201; j++)
            {
                grid[i][j] = coords({grid_xs[i], grid_ys[j], 0});
            }
        }

        std::array<std::array<double,201>,201> TA_grid;
        std::mutex TA_grid_mutex;
        for (uint16_t i=0; i<201; i++)
        {
            TA_grid[i] = std::array<double,201>();
            for (uint16_t j=0; j<201; j++)
            {
                TA_grid[i][j] = 0.0;
            }
        }

        TAs_file.open("TA_ave_no_shift.wl");

        TAs_file<<"TAaveNoShift"<<num_nuclei<<" = ";
        std::vector<uint64_t> indexes((num_nuclei/2)-1);
        std::iota(indexes.begin(), indexes.end(), 0); //generates the list as {0,1,2,3,...}

        calcs::generate_nuclei
        (
            nuc_params, 
            1000, 
            0, 
            eng, 
            radial_sampler, 
            false, 
            false
        );

        std::for_each
        (
            std::execution::par, 
            indexes.begin(), 
            indexes.end(), 
            [&](const uint64_t index) 
            {

                std::unique_lock<std::mutex> lock_rad(radial_sampler_mutex);
                auto [pro, tar] = calcs::generate_nuclei
                (
                    nuc_params, 
                    1000, 
                    0, 
                    eng, 
                    radial_sampler, 
                    false, 
                    false
                );
                lock_rad.~unique_lock();

                std::array<std::array<double,201>,201> TA_grid_dummy;
                for (uint16_t i=0; i<201; i++)
                {
                    TA_grid_dummy[i] = std::array<double,201>();
                    for (uint16_t j=0; j<201; j++)
                    {
                        TA_grid_dummy[i][j] = calcs::calculate_tA(grid[i][j], pro, Tp);
                        TA_grid_dummy[i][j] += calcs::calculate_tA(grid[i][j], tar, Tp);
                        TA_grid_dummy[i][j] /= static_cast<double>(num_nuclei);
                    }
                }

                {
                    const std::lock_guard<std::mutex> lock(TA_grid_mutex);
                    for (uint16_t i=0; i<201; i++)
                    {
                        for (uint16_t j=0; j<201; j++)
                        {
                            TA_grid[i][j] += TA_grid_dummy[i][j];
                        }
                    }
                }
            }
        );

        TAs_file<<"{";
        std::for_each(TA_grid.begin(), TA_grid.end()-1, [&TAs_file](const std::array<double,201> &ys) 
        { 
            TAs_file<<"{";
            std::for_each(ys.begin(), ys.end()-1, [&TAs_file](const double &val) 
            { 
                TAs_file<<val<<",";
            });
            TAs_file<<ys.back()<<"},\n    ";
        });
        TAs_file<<"{";
        std::for_each(TA_grid.back().begin(), TA_grid.back().end()-1, [&TAs_file](const double &val) 
        { 
            TAs_file<<val<<",";
        });
        TAs_file<<TA_grid.back().back()<<"}};\n\n";
        TAs_file.close();
    }
};

#endif // HIGH_LEVEL_CALCS_HPP
