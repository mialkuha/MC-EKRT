//Copyright (c) 2022 Mikko Kuha 
//TODO implement M-tree (https://github.com/erdavila/M-Tree ? Boost::Graph?)
//     Probably more like NN-tree.... 

#include <algorithm>
#include <atomic>
#include <chrono>
#include <csignal>
#include <execution>
#include <fstream>
#include <gsl/gsl_sf_expint.h>
#include <iostream>
#include <mutex>
#include <random>
#include <tuple>
#include <unordered_set>
#include <vector>

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshorten-64-to-32"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "LHAPDF/GridPDF.h"
#pragma GCC diagnostic pop

#include "ars.hpp"
#include "generic_helpers.hpp"
#include "high_level_calcs.hpp"
#include "histo.hpp"
#include "io_helpers.hpp"
#include "nucleus_generator.hpp"
#include "pqcd.hpp"
#include "typedefs.hpp"

/**
 * @brief Find k_T^2 cutoff that fits sigma_jet to a given 
 * target number using secant method.
 * 
 * @param kt02 k_T^2 cutoff will be saved to this (GeV^2)
 * @param mand_s mandelstam s for the process (GeV^2)
 * @param target the target value for sigma_jet (mb)
 * @param p_p_pdf pointer to LHAPDF PDF-object
 * @param jet_params the struct of jet parameters
 * @param verbose if true, prints all the intermediate results
 * @return double the found k_T^2 cutoff (GeV^2)
 */
auto fit_sigma_jet_p0_cutoff
(
    double &kt02, 
    const double &mand_s, 
    const double &target, 
    std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
    const pqcd::sigma_jet_params &jet_params, 
    const bool &verbose=true
) noexcept -> double
{
    double sigma_jet=0.0;
    kt02 = 2.0;

    auto difference_to_target = [&](const double &_kt02)
    {
        return pqcd::calculate_sigma_jet(p_p_pdf, &mand_s, &_kt02, jet_params) - target;
    };

    helpers::secant_method(&kt02, difference_to_target, 1e-3, &sigma_jet);

    if (verbose) std::cout<<sqrt(kt02)<<' '<<sigma_jet+target<<' '<<target<<std::endl;
    
    return kt02;
}

/**
 * @brief Find scale choice Q/k_T that fits sigma_jet to a given 
 * target number, for a given target k_T^2 cutoff, using secant method.
 * 
 * @param kt02 The target k_T^2 cutoff value (GeV^2)
 * @param mand_s mandelstam s for the process (GeV^2)
 * @param target the target value for sigma_jet (mb)
 * @param p_p_pdf pointer to LHAPDF PDF-object
 * @param jet_params the struct of jet parameters
 * @param verbose if true, prints all the intermediate results
 * @return double the found scale choice Q/k_T
 */
auto find_sigma_jet_cutoff_Q
(
    const double &kt02, 
    const double &mand_s, 
    const double &target, 
    std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
    pqcd::sigma_jet_params jet_params, 
    const bool &verbose=true
) noexcept -> double
{
    double sigma_jet=0.0;
    double scalar = 1.0;

    auto difference_to_target = [&](const double &scalar_)
    {
        double kt02dummy = 4.0;
        auto jet_params_ = pqcd::sigma_jet_params(
        /*d_params=                 */jet_params.d_params,
        /*scale_choice=             */jet_params.scale_c,
        /*scalar=                   */scalar_,
        /*use_ses=                  */jet_params.use_ses);
        auto kt02_ = fit_sigma_jet_p0_cutoff(kt02dummy, mand_s, target, p_p_pdf, jet_params_);
        return kt02_ - kt02;
    };

    helpers::secant_method(&scalar, difference_to_target, 1e-3, &sigma_jet);

    if (verbose) std::cout<<scalar<<' '<<sigma_jet+target<<' '<<target<<std::endl;
    
    return scalar;
}

/**
 * @brief struct for the end state filtering functions. In the constructor 
 * calculates the formation time t0 = y/max(u,t) = (1/Q)*(E/Q) (in the collider frame) 
 * for the jets.
 * 
 */
struct dijet_with_ns
{
    dijet_specs dijet;
    nucleon * pro_nucleon;
    nucleon * tar_nucleon;
    double t01;
    double t02;
    double t0; //Larger of the formation times := the formation time of the process
    dijet_with_ns(dijet_specs dijet_, nucleon * pro_nucleon_, nucleon * tar_nucleon_)
        : dijet(std::move(dijet_)), pro_nucleon(pro_nucleon_), tar_nucleon(tar_nucleon_)
    {
        auto y1 = dijet.y1;
        auto y2 = dijet.y2;
        auto max_t_u = dijet.kt*(1+std::exp(-std::abs(y1-y2)));
        t01 = std::cosh(std::abs(y1))/max_t_u;
        t02 = std::cosh(std::abs(y2))/max_t_u;
        t0 = std::max(t01, t02);
    }
};

/**
 * @brief Checks if the candidate circle overlaps with the previous circles. If
 * it does not, add it to the collection of accepted circles.
 * 
 * @param cand_circle (x,y,r,n), where x and y are coordinates (fm), r is 
 * radius (fm) and n is for supporting, say, two circles overlapping but no 
 * more (not implemented).
 * @param final_circles A collection of the previously accepted circles.
 * @return true The circle was accepted and added to the collection.
 * @return false The circle was not accepted.
 */
auto check_and_place_circle_among_others
(
    std::tuple<double, double, double, uint_fast16_t> cand_circle, 
    std::vector<std::tuple<double, double, double, uint_fast16_t> > &final_circles
) noexcept -> bool
{
    auto & [cand_x, cand_y, cand_r, cand_overlap] = cand_circle;
    
    for (auto & [circ_x, circ_y, circ_r, circ_overlap] : final_circles)
    {
        //If the any of the circles radii overlap with the candidate, just return false
        if ( pow((circ_x-cand_x),2) + pow((circ_y-cand_y),2) < pow((circ_r+cand_r),2) )
        {
            return false;
        }
    }

    //None overlapped, emplace the circle among others and return true
    final_circles.emplace_back(std::move(cand_circle));
        
    return true;
}

/**
 * @brief Sample the coordinates for a given dijet from the Gaussian distribution
 * that is the product of the two Gaussian distributions that are associated with
 * the mother nucleons.
 * 
 * @param dijet The dijet in question
 * @param proton_width Proton width (fm)
 * @param normal_dist Normal distribution object
 * @param random_generator Randdom generator object
 * @return std::tuple<double, double, double> (x,y,z=0) coordinates of the dijet (fm)
 */
auto throw_location_for_dijet
(
    const dijet_with_ns &dijet,
    const double &proton_width,
    std::normal_distribution<double> &normal_dist, 
    std::shared_ptr<std::mt19937> random_generator
) noexcept -> std::tuple<double, double, double>
{
    auto param = std::normal_distribution<double>::param_type{0., proton_width};
    auto dx = normal_dist(*random_generator,param);
    auto dy = normal_dist(*random_generator,param);
    auto dz = 0.0;

    auto retx = 0.5*(dijet.pro_nucleon->co.x + dijet.tar_nucleon->co.x + M_SQRT2*dx);
    auto rety = 0.5*(dijet.pro_nucleon->co.y + dijet.tar_nucleon->co.y + M_SQRT2*dy);

    return std::make_tuple(retx, rety, dz);
}

auto filter_end_state
(
    std::vector<nn_coll> &binary_collisions, 
    std::vector<dijet_with_coords> &filtered_scatterings,
    std::shared_ptr<std::mt19937> random_generator,
    const bool mom_cons = false,
    /*const bool mom_cons_local = false, NOT YET IMPLEMENTED*/
    const bool saturation = false,
    const bool deplete_nucleons = false,
    const double sqrt_s = 0,
    const double &maximum_overlap = 2.0,
    const double &proton_width = 1.0,
    /*const double &sigma_inel = 1.0, NOT YET IMPLEMENTED*/
    const bool calculate_tata = true,
    const std::vector<nucleon> &pro = {}, 
    const std::vector<nucleon> &tar = {}, 
    const std::function<double(const double&)> Tpp = nullptr
) noexcept -> void
{
    std::vector<dijet_with_ns> candidates;
    candidates.reserve(binary_collisions.size()*10); //10 events on average is just an overhead guess

    std::unordered_map<nucleon*, double> x1s;
    std::unordered_map<nucleon*, double> x2s;
    std::unordered_set<nucleon*> depleted_pro;
    std::unordered_set<nucleon*> depleted_tar;
    

    std::normal_distribution<double> normal_dist(0,0);

    for (auto &col : binary_collisions)
    {
        for (auto &dij : col.dijets)
        {
            candidates.emplace_back(std::move(dij), col.projectile, col.target);
        }
    }
    candidates.shrink_to_fit();
    
    std::vector<std::tuple<double, double, double, uint_fast16_t> > final_circles;
    final_circles.reserve(candidates.size());
    filtered_scatterings.reserve(filtered_scatterings.size()+candidates.size());

    //std::sort(collision_candidates.begin(), collision_candidates.end(), //Sort the candidates so that the one with the biggest kt is first
    //          [](colls_with_ns &s1, colls_with_ns &s2) { return (s1.kt > s2.kt); });
    std::sort(candidates.begin(), candidates.end(), //Sort the candidate events so that the one with the smallest t0 is first
              [](dijet_with_ns &s1, dijet_with_ns &s2) { return (s1.t0 < s2.t0); });
    
    double kt;
    double y1, y2;
    double i_x1_sum_to_be, i_x2_sum_to_be;
    double tata = 0.0;

    for (auto & cand : candidates)
    {
        kt = cand.dijet.kt;
        y1 = cand.dijet.y1;
        y2 = cand.dijet.y2;

        if (mom_cons)
        {
            //double CONSERVATION

            if (deplete_nucleons && (depleted_pro.contains(cand.pro_nucleon) || depleted_tar.contains(cand.tar_nucleon)))
            {
                continue;
            }

            auto x1 = (kt / sqrt_s) * (exp(y1) + exp(y2));
            auto x2 = (kt / sqrt_s) * (exp(-y1) + exp(-y2));

            i_x1_sum_to_be = x1;
            i_x2_sum_to_be = x2;

            auto i_x1_sum = x1s.find(cand.pro_nucleon);
            if (i_x1_sum != x1s.end())
            {
                i_x1_sum_to_be = i_x1_sum->second + x1;

                if (i_x1_sum_to_be > 1.0) //Energy budget broken --> discard
                {
                    if (deplete_nucleons)
                    {
                        depleted_pro.insert(cand.pro_nucleon);
                    }
                    continue;
                }
            }

            auto i_x2_sum = x2s.find(cand.tar_nucleon);
            if (i_x2_sum != x2s.end())
            {
                i_x2_sum_to_be = i_x2_sum->second + x2;

                if (i_x2_sum_to_be > 1.0) //Energy budget broken --> discard
                {
                    if (deplete_nucleons)
                    {
                        depleted_tar.insert(cand.tar_nucleon);
                    }
                    continue;
                }
            }
        }
        
        auto [cand_x, cand_y, cand_z] = throw_location_for_dijet(cand, proton_width, normal_dist, random_generator);

        if (saturation)
        {
            //SATURATION

            auto cand_circle = std::make_tuple(cand_x, cand_y, 1/(cand.dijet.kt*maximum_overlap*FMGEV), 1);

            if (!check_and_place_circle_among_others(std::move(cand_circle), final_circles))
            {
                continue; //Did not fit into saturated PS --> discard
            }
        }

        if (calculate_tata)
        {
            nucleon dummy{coords{cand_x, cand_y, cand_z}, 0};
            tata = calcs::calculate_sum_tpp(dummy, pro, Tpp) * calcs::calculate_sum_tpp(dummy, tar, Tpp);
        }

        if (mom_cons)
        {
            x1s.insert_or_assign(cand.pro_nucleon, i_x1_sum_to_be);
            x2s.insert_or_assign(cand.tar_nucleon, i_x2_sum_to_be);
        }
        filtered_scatterings.push_back({cand.dijet, coords{cand_x, cand_y, cand_z}, cand.t01, cand.t02, tata});
    }

    binary_collisions.clear();
}

//--------------------------------------------//
// THESE ARE A HACK TO SAVE THE WORK THIS FAR //
// IN THE CASE OF SUDDEN TERMINATION          //
//--------------------------------------------//
volatile std::atomic_bool user_aborted = false;
volatile std::atomic_bool g_bug_bool = false;
void abort_handler(int num)
{
    user_aborted = true;
    std::cout<<std::endl<<"Abort called "<<num<<std::endl;
}

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cout<<"Provide parameter file name as a parameter, see params_template"<<std::endl;
        return 1;
    }

    auto 
    [
        name_postfix,
        g_is_pp,
        g_is_pa,
        g_is_aa,
        g_use_npdfs,
        g_use_snpdfs,
        g_is_mom_cons,
        g_are_ns_depleted,
        g_is_saturation,
        g_is_mc_glauber,
        desired_N_events,
        b_min,
        b_max,
        pt0,
        K_factor,
        K_sat
    ] = io::read_conf(std::string(argv[1]));

    //Prints a lot of debugging info
    bool verbose = false;

    if (verbose) std::cout<<"Initializing..."<<std::flush;
    
    //General parameters for the simulation
    const bool    read_sigmajets_from_file  = false,
                  end_state_filtering       = true, 
                  sigma_inel_from_sigma_jet = true, 
                  save_endstate_jets        = true,
                  save_events_plaintext     = false,
                  calculate_end_state       = true,
                  calculate_tata            = true,
                  reduce_nucleon_energies   = false;

    uint_fast32_t AA_events_done            = 0;
    std::mutex AA_events_mutex;

    std::cout<<std::endl<<std::endl<<"Doing the run "<<name_postfix<<std::endl;

    //auto eng = std::make_shared<std::mt19937>(static_cast<ulong>(1));
    auto eng = std::make_shared<std::mt19937>(static_cast<ulong>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<double> unirand{0.0, 1.0};

    //Parameters for the nuclei
    const double rad_min=0,
                  rad_max=30;
    const std::function<double(const double&)> rad_pdf{[](const double & x)
    {
        return x*x/(1+exp((x-6.624)/0.549));
    }};
    auto radial_sampler{std::make_shared<ars>(rad_pdf, rad_min, rad_max)};
    //The adaptive algorithm in the sampler is not thread-safe,
    //so to run the program multithreaded let's first saturate the sampler
    do //while (radial_sampler->is_adaptive()) 
    {
        radial_sampler->throw_one(*eng);
    } while (radial_sampler->is_adaptive());

    nucleus_generator::nucleus_params nuc_params = 
    {
        /* .NA=                   */(g_is_aa)? 208u : 1u,
        /* .ZA=                   */(g_is_aa)? 82u : 1u,  
        /* .NB=                   */(g_is_pp)? 1u : 208u, 
        /* .ZB=                   */(g_is_pp)? 1u : 82u,  
        /* .min_distance=         */0.4, 
        /* .shift_cms=            */true, 
        /* .correct_overlap_bias= */true
    };
    
    //Parameters for the hard collisions
    const double proton_width = 0.573;
    const double proton_width_2 = pow(proton_width, 2);
    const std::function<double(const double&)> Tpp{[&proton_width_2](const double &bsquared)
    {
        return exp(-bsquared / (4 * proton_width_2)) / (40 * M_PI * proton_width_2); // 1/fm² = mb/fm² * 1/mb = 0.1 * 1/mb
    }}; 
    double sigma_inel_for_glauber       = 70;//mb
    const double sqrt_s                 = 5020;//GeV
    const double mand_s                 = pow(sqrt_s, 2);//GeV^2
    double kt0                          = pt0;//2.728321;//GeV
    double kt02                         = pow(kt0, 2);//GeV^2
    double ylim                         = static_cast<double>(log(sqrt_s / kt0));
    auto p_pdf = std::make_shared<LHAPDF::GridPDF>("CT14lo", 0);

    std::vector<io::Coll> collisions_for_reporting;
    std::vector<io::Coll> collisions_for_reporting_midrap;
    std::mutex colls_mutex; 
    std::vector<double> kt_bins{helpers::loglinspace(kt0, 100, 16u)};
    auto kt_v_dummy = helpers::loglinspace(200, sqrt_s/2.0, 5u);
    for (auto ktdummy : kt_v_dummy)
    {
        kt_bins.push_back(ktdummy);
    }
    const std::vector<double> y_bins{helpers::linspace(-ylim, ylim, 40u)};
    const std::vector<double> b_bins{helpers::linspace(b_min, b_max, 21u)};
    std::vector<double> et_bins{helpers::loglinspace(2*kt0, 30000, 31u)};

    //sigma_jet parameters
    pqcd::diff_sigma::params diff_params = pqcd::diff_sigma::params(
    /*projectile_with_npdfs=    */(g_is_aa && g_use_npdfs),
    /*target_with_npdfs=        */(!g_is_pp && g_use_npdfs),
    /*isoscalar_projectile=     */false,
    /*isoscalar_target=         */false,
    /*npdfs_spatial=            */g_use_snpdfs,
    /*npdf_setnumber=           */1,
    /*K_factor=                 */K_factor,
    /*A=                        */(g_is_aa)? 208u : 1u, //Pb 
    /*B=                        */(g_is_pp)? 1u : 208u, //Pb
    /*ZA=                       */(g_is_aa)? 82u : 1u,  //Pb
    /*ZB=                       */(g_is_pp)? 1u : 82u   //Pb
    /*p_n_pdf=                  */
    /*rA_spatial=               */
    /*rB_spatial=               */);
    pqcd::sigma_jet_params jet_params = pqcd::sigma_jet_params(
    /*d_params=                 */diff_params,
    /*scale_choice=             */pqcd::scaled_from_kt,
    /*scalar=                   */1.0,
    /*use_ses=                  */false);

    auto
    [
        dijet_norm,
        power_law,
        sigma_jet,
        mand_s_array,
        envelope_maximum,
        sqrt_s_array
    ] = 
    calcs::prepare_sigma_jets
    (
        reduce_nucleon_energies,
        read_sigmajets_from_file,
        p_pdf, 
        mand_s,
        sqrt_s,
        kt02, 
        kt0,
        jet_params
    );
    
    if (sigma_inel_from_sigma_jet)
    {
        auto dummy = std::get<double>(sigma_jet) / (4 * M_PI * proton_width_2);
        dummy = dummy / 10; //mb -> fm²

        sigma_inel_for_glauber = (4 * M_PI * proton_width_2) * (M_EULER + std::log(dummy) + gsl_sf_expint_E1(dummy)) * 10; //1 mb = 0.1 fm^2
        std::cout<<"sigma_inel = "<<sigma_inel_for_glauber<<std::endl;
    }

    AA_collision_params coll_params
    {
    /*mc_glauber_mode=          */g_is_mc_glauber,
    /*pp_scattering=            */g_is_pp,
    /*pA_scattering=            */g_is_pa,
    /*spatial_pdfs=             */g_use_snpdfs,
    /*calculate_end_state=      */calculate_end_state,
    /*reduce_nucleon_energies=  */reduce_nucleon_energies,
    /*sigma_inel_for_glauber=   */sigma_inel_for_glauber,
    /*Tpp=                      */Tpp,
    /*normalize_to=             */B2_normalization_mode::inelastic,
    /*sqrt_s=                   */sqrt_s,
    /*energy_threshold=         */kt0
    };

    if (verbose) std::cout<<"Done!"<<std::endl;

    auto cmpLambda = [](const io::Coll &lhs, const io::Coll &rhs) { return io::compET(lhs, rhs); };
    std::map<io::Coll, std::vector<dijet_with_coords>, decltype(cmpLambda)> colls_scatterings(cmpLambda);

//    OBSERVABLES TO BE SAVED
//
//    histo_2d jets{kt_bins, y_bins};
//    histo_2d dijets{kt_bins, y_bins};
//    std::mutex jets_mutex; 
//    std::mutex dijets_mutex;
//
//    histo_1d dETdy{y_bins};
//    histo_1d dEdy{y_bins};
//    histo_1d dNdy{y_bins};
//    histo_1d dNdET{et_bins};
//    std::mutex dETdy_mutex; 
//    std::mutex dEdy_mutex;
//    std::mutex dNdy_mutex;
//    std::mutex dNdET_mutex;
//
//    histo_1d dETdeta{y_bins};
//    histo_1d dEdeta{y_bins};
//    std::mutex dETdeta_mutex; 
//    std::mutex dEdeta_mutex;
// 
//    histo_1d dETdb{b_bins};
//    histo_1d dEdb{b_bins};
//    std::mutex dETdb_mutex; 
//    std::mutex dEdb_mutex;


    std::ofstream total_energy;
    std::mutex total_energy_mutex; 
 
    total_energy.open("total_energies_"+name_postfix+".dat");
    total_energy << "///Sum E_T Sum E" << std::endl;
    
    std::ofstream event_file;
    std::mutex event_file_mutex; 

    if (save_events_plaintext)
    {
        event_file.open("event_log_"+name_postfix+".dat");

        event_file << "///AA b: ["<<b_min<<", "<<b_max<<"] fm"<<std::endl;
        event_file << "///"<<std::endl;
        event_file << "///Nucleus params:"<<std::endl;
        event_file << "///Nucleus r: ["<<rad_min<<", "<<rad_max<<"] fm"<<std::endl;
        event_file << "///Correct overlap bias: "<<nuc_params.correct_overlap_bias<<std::endl;
        event_file << "///Nucleon min_distance: "<<nuc_params.min_distance<<" fm"<<std::endl;
        event_file << "///Shift CMS: "<<nuc_params.shift_cms<<std::endl;
        event_file << "///"<<std::endl;
        event_file << "///Collision params:"<<std::endl;
        event_file << "///Proton width: "<<sqrt(proton_width_2)<<" fm"<<std::endl;
        event_file << "///Sqrt(s): "<<sqrt_s<<" GeV"<<std::endl;
        event_file << "///kt0: "<<kt0<<" GeV"<<std::endl;
        event_file << "///Projectile with npdfs: "<<diff_params.projectile_with_npdfs<<std::endl;
        event_file << "///Target with npdfs: "<<diff_params.target_with_npdfs<<std::endl;
        event_file << "///double npdfs: "<<coll_params.spatial_pdfs<<std::endl;
        event_file << "///Projectile isoscalar: "<<diff_params.isoscalar_projectile<<std::endl;
        event_file << "///Target isoscalar: "<<diff_params.isoscalar_target<<std::endl;
        event_file << "///Normalization mode: "<<coll_params.normalize_to<<std::endl;
        event_file << "///Scale choice: "<<jet_params.scale_c<<std::endl;
        event_file << "///Scalar: "<<jet_params.scalar<<std::endl;
        event_file << "///"<<std::endl;
        event_file << "///End state params:"<<std::endl;
        event_file << "///Calculate end state: "<<coll_params.calculate_end_state<<std::endl;
        event_file << "///Power law: "<<power_law<<std::endl;
        event_file << "///"<<std::endl;
    }

    std::vector<uint_fast64_t> event_indexes(desired_N_events);
    std::iota(event_indexes.begin(), event_indexes.end(), 0); //generates the list as {0,1,2,3,...}
    std::atomic<uint_fast64_t> running_count{desired_N_events};

    std::signal(SIGINT, abort_handler);

    std::set_terminate([](){
        std::cout << std::endl << "Unhandled exception" << std::endl;
        g_bug_bool = true;
    });

    try
    {
        std::find_if
        (
            std::execution::par, 
            event_indexes.begin(), 
            event_indexes.end(), 
            [&,&b_max=b_max,
             &verbose=verbose,
             &g_is_mom_cons=g_is_mom_cons,
             &g_is_saturation=g_is_saturation,
             &g_are_ns_depleted=g_are_ns_depleted,
             &sigma_jet=sigma_jet,
             &power_law=power_law,
             &envelope_maximum=envelope_maximum,
             &K_sat=K_sat,
             &colls_scatterings=colls_scatterings
            ](const uint_fast64_t index) 
            {
                static_cast<void>(index);
                do //while (g_bug_bool)
                {
                    uint_fast32_t NColl = 0;
                    std::vector<nn_coll> binary_collisions;
                    std::vector<dijet_with_coords> filtered_scatterings;
                    double impact_parameter;

                    std::vector<nucleon> pro, tar;
                    uint_fast16_t times_discarded = 0;

                    //Demand at least one hard scattering
                    do //while (NColl<1)
                    {
                        //Keep generating nuclei until there are nucleons close enough to each other
                        //so that a collision is probable
                        bool bugged, probably_collided;
                        // "ball" diameter = distance at which two nucleons interact in MC Glauber
                        const double d2 = sigma_inel_for_glauber/(M_PI*10.0); // in fm^2
                        do //while (!probably_collided || bugged)  
                        {
                            //B^2 from a uniform distribution
                            impact_parameter = sqrt(b_min*b_min + unirand(*eng)*(b_max*b_max-b_min*b_min));
                    
                            times_discarded++;
                            if (times_discarded > 1000)
                            {
                                std::cout<<std::endl<<"Generated nuclei discarded over 1000 times. "
                                         <<"Check impact parameters and/or collsion probabilities."
                                         <<std::endl;
                                times_discarded = 0;
                            }

                            bugged = false;
                            probably_collided = false;
                            try
                            {
                                auto [pro_dummy, tar_dummy] = calcs::generate_nuclei
                                (
                                    nuc_params, 
                                    sqrt_s, 
                                    impact_parameter, 
                                    eng, 
                                    radial_sampler, 
                                    verbose
                                );
                                pro = std::move(pro_dummy);
                                tar = std::move(tar_dummy);
                            }
                            catch(const std::exception& e)
                            {
                                std::cout << e.what() << " in main, trying again"<<std::endl;
                                bugged = true;
                            }

                            for (auto A : pro)
                            {
                                for (auto B : tar)
                                {
                                    // "ball" diameter = distance at which two nucleons interact
                                    const double dij2 = A.calculate_bsquared(B);

                                    if (dij2 <= d2) //Probably at least one collision
                                    {
                                        probably_collided = true;
                                        continue;
                                    }
                                }
                                if (probably_collided)
                                {
                                    continue;
                                }
                            }

                        } while (!probably_collided || bugged);

                        binary_collisions.clear();
                        if (verbose) std::cout<<"impact_parameter: "<<impact_parameter<<std::endl;

                        calcs::collide_nuclei
                        (
                            pro, 
                            tar, 
                            binary_collisions, 
                            sigma_jet,
                            unirand, 
                            eng, 
                            coll_params, 
                            jet_params,
                            kt0,
                            p_pdf,
                            power_law,
                            envelope_maximum,
                            verbose
                        );
                        
                        NColl = static_cast<uint_fast32_t>(binary_collisions.size());
                    } while (NColl<1);
                    

                    double sum_ET = 0;
                    double sum_ET_midrap = 0;

                    if(!end_state_filtering && save_events_plaintext)
                    {
                        const std::lock_guard<std::mutex> lock(event_file_mutex);
                        io::save_event(event_file, pro, tar, impact_parameter);
                    }
                    else if (end_state_filtering)
                    {
                        double max_overlap = K_sat;
                        double sum_E = 0;
                        filter_end_state
                        (
                            binary_collisions, 
                            filtered_scatterings, 
                            eng, 
                            g_is_mom_cons, 
                            /*g_is_mom_cons_local, FEATURE NOT IMPLEMENTED*/
                            g_is_saturation, 
                            g_are_ns_depleted, 
                            sqrt_s,
                            max_overlap,
                            proton_width,
                            /*sigma_inel_for_glauber, FEATURE NOT IMPLEMENTED*/
                            calculate_tata,
                            pro,
                            tar,
                            coll_params.Tpp
                        );

                        //std::vector<std::tuple<double, double> > new_jets;
                        //std::vector<std::tuple<double, double> > new_dijets;
                        //std::vector<std::tuple<double, double> > new_ET_y;
                        //std::vector<std::tuple<double, double> > new_E_y;
                        //std::vector<std::tuple<double, double> > new_N_y;
                        //std::vector<std::tuple<double, double> > new_ET_eta;
                        //std::vector<std::tuple<double, double> > new_E_eta;
                        
                        for (auto e_co : filtered_scatterings)
                        {
                            auto e = e_co.dijet;

                            //new_jets.emplace_back(e.kt, e.y1);
                            //new_jets.emplace_back(e.kt, e.y2);
                            //
                            //new_dijets.emplace_back(e.kt, 0.5*(e.y1+e.y2));
                            //
                            //new_ET_y.emplace_back(e.y1, e.kt);
                            //new_ET_y.emplace_back(e.y2, e.kt);
                            //
                            //new_ET_eta.emplace_back(0.5*(e.y1+e.y2), 2*e.kt);
                            //
                            //new_E_y.emplace_back(e.y1, e.kt*cosh(e.y1));
                            //new_E_y.emplace_back(e.y2, e.kt*cosh(e.y2));
                            //
                            //new_N_y.emplace_back(e.y1, 1);
                            //new_N_y.emplace_back(e.y2, 1);
                            //
                            //new_E_eta.emplace_back(0.5*(e.y1+e.y2), e.kt*(cosh(e.y1) + cosh(e.y2)));
                            
                            sum_ET += 2*e.kt;
                            sum_E += e.kt*(cosh(e.y1) + cosh(e.y2));

                            if (e.y1 >= -0.5 && e.y1 <= 0.5)
                            {
                                sum_ET_midrap += e.kt;
                            }
                            if (e.y2 >= -0.5 && e.y2 <= 0.5)
                            {
                                sum_ET_midrap += e.kt;
                            }
                        }

                        //{
                        //    const std::lock_guard<std::mutex> lock(jets_mutex);
                        //    jets.add(new_jets);
                        //}
                        //
                        //{
                        //    const std::lock_guard<std::mutex> lock(dijets_mutex);
                        //    dijets.add(new_dijets);
                        //}
                        //
                        //{
                        //    const std::lock_guard<std::mutex> lock(dETdy_mutex);
                        //    dETdy.add(new_ET_y);
                        //}
                        //
                        //{
                        //    const std::lock_guard<std::mutex> lock(dEdy_mutex);
                        //    dEdy.add(new_E_y);
                        //}
                        //
                        //{
                        //    const std::lock_guard<std::mutex> lock(dNdy_mutex);
                        //    dNdy.add(new_N_y);
                        //}
                        //
                        //{
                        //    const std::lock_guard<std::mutex> lock(dETdeta_mutex);
                        //    dETdeta.add(new_ET_eta);
                        //}
                        //
                        //{
                        //    const std::lock_guard<std::mutex> lock(dEdeta_mutex);
                        //    dEdeta.add(new_E_eta);
                        //}
                        //
                        //{
                        //    const std::lock_guard<std::mutex> lock(dNdET_mutex);
                        //    dNdET.add(std::make_tuple(sum_ET, 1));
                        //}
                        //
                        //{
                        //    const std::lock_guard<std::mutex> lock(dETdb_mutex);
                        //    dETdb.add(std::make_tuple(impact_parameter, sum_ET));
                        //}
                        //
                        //{
                        //    const std::lock_guard<std::mutex> lock(dEdb_mutex);
                        //    dEdb.add(std::make_tuple(impact_parameter, sum_E));
                        //}

                        {
                            const std::lock_guard<std::mutex> lock(total_energy_mutex);
                            total_energy << sum_ET << ' ' << sum_E << std::endl;
                        }

                        if (save_events_plaintext)
                        {
                            const std::lock_guard<std::mutex> lock(event_file_mutex);
                            io::save_event(event_file, pro, tar, impact_parameter, filtered_scatterings);
                        }
                    }

                    {
                        const std::lock_guard<std::mutex> lock(AA_events_mutex);
                        AA_events_done++;
                        if (AA_events_done % 100 == 0 )
                        {
                            std::cout <<"\rA+A collisions calculated: " << AA_events_done << std::flush;
                        }
                    }

                    uint_fast32_t Npart=0;
                    for (auto &A : pro)
                    {
                        if (A.wounded)
                        {
                            Npart++;
                        }
                    }
                    for (auto &B : tar)
                    {
                        if (B.wounded)
                        {
                            Npart++;
                        }
                    }

                    {
                        const std::lock_guard<std::mutex> lock(colls_mutex);
                        for (uint_fast8_t i=0; i<4; i++)
                        {
                            io::Coll coll(NColl, Npart, 2*filtered_scatterings.size(), impact_parameter, sum_ET);
                            io::Coll coll_midrap(NColl, Npart, 2*filtered_scatterings.size(), impact_parameter, sum_ET_midrap);
                            collisions_for_reporting.push_back(coll);
                            collisions_for_reporting_midrap.push_back(coll_midrap);
                            
                            if (save_endstate_jets)
                            {
                                colls_scatterings.insert({coll, filtered_scatterings});
                            }
                        }
                    }
            
                } while (g_bug_bool);
                bool ret_value = user_aborted;
                return ret_value;
            }
        );
    }
    catch(const std::exception& e)
    {
        std::cout<<std::endl<<"Threw " << e.what() <<std::endl;
    }
    std::cout<<" ...done!"<<std::endl;
    


    //io::print_histos
    //(
    //    name_postfix,
    //    jets,
    //    dijets,
    //    dETdy,
    //    dEdy,
    //    dNdy,
    //    dNdET,
    //    dETdeta,
    //    dEdeta,
    //    dETdb,
    //    dEdb,
    //    dijet_norm,
    //    AA_events_done
    //);
    
    uint_fast8_t nBins = 18;
    double binsLow[] = {0., 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 0.0, 0.0, 0.0, 0.0};
    double binsHigh[] = {0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 0.05, 0.1, 0.8, 1.0};
    std::ofstream glauber_report_file;
    std::string g_name{"g_report_"+name_postfix+".dat"};
    std::string g_name_midrap{"g_report_midrap_"+name_postfix+".dat"};
                                              
    glauber_report_file.open(g_name, std::ios::out);
    io::mc_glauber_style_report(collisions_for_reporting, sigma_inel_for_glauber, desired_N_events, nBins, binsLow, binsHigh, glauber_report_file);
    glauber_report_file.close();
    glauber_report_file.open(g_name_midrap, std::ios::out);
    io::mc_glauber_style_report(collisions_for_reporting_midrap, sigma_inel_for_glauber, desired_N_events, nBins, binsLow, binsHigh, glauber_report_file);
    glauber_report_file.close();


    if (save_endstate_jets)
    {
        const std::array<std::tuple<double,double>, 3> centBins{std::tuple<double,double>{0.0, 0.05},
                                                                std::tuple<double,double>{0.25, 0.3},
                                                                std::tuple<double,double>{0.6, 0.8}};

        std::string name_pfs{name_postfix+".dat"};

        std::ofstream jet_file;

        for (auto [centLow, centHigh] : centBins)
        {
            std::stringstream jetsname{""};
            jetsname<<"jets_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            jet_file.open(jetsname.str(), std::ios::out | std::ios::binary);

            histo_1d dETdy_by_cent{y_bins};
            histo_1d dEdy_by_cent{y_bins};
            std::vector<std::tuple<double, double> > new_ET_y;
            std::vector<std::tuple<double, double> > new_E_y;

            uint_fast64_t N_evts_tot = colls_scatterings.size();
            // Make sure that no rounding downwards.
            double eps = 0.1/static_cast<double>(N_evts_tot);

            uint_fast64_t lower_ind = static_cast<uint_fast64_t>(centLow*static_cast<double>(N_evts_tot)+eps);
            uint_fast64_t upper_ind = static_cast<uint_fast64_t>(centHigh*static_cast<double>(N_evts_tot)+eps);
            uint_fast64_t n_in_bin = upper_ind - lower_ind;

            //total number of events in this bin
            jet_file.write(reinterpret_cast<char*>(&n_in_bin), sizeof n_in_bin);

            auto it = colls_scatterings.crbegin();
            std::advance(it, lower_ind);

            for (uint_fast64_t ii = 0; ii<n_in_bin; it++, ii++)
            {
                for (auto e_co : it->second)
                {
                    auto e = e_co.dijet;

                    new_ET_y.emplace_back(e.y1, e.kt);
                    new_ET_y.emplace_back(e.y2, e.kt);

                    new_E_y.emplace_back(e.y1, e.kt*cosh(e.y1));
                    new_E_y.emplace_back(e.y2, e.kt*cosh(e.y2));
                }
                io::append_single_coll_binary(jet_file, it->second, unirand, eng);
            }
            jet_file.close();
            std::cout<<n_in_bin<<std::endl;

            dETdy_by_cent.add(new_ET_y);
            dEdy_by_cent.add(new_E_y);

            std::stringstream outname{""};

            outname<<"dEdy_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            io::print_1d_histo
            (
                dEdy_by_cent, 
                outname.str(), 
                1.0/ static_cast<double>(n_in_bin),
                false
            );
            outname.seekp(0);
            outname<<"dETdy_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            io::print_1d_histo
            (
                dETdy_by_cent,
                outname.str(), 
                1.0/ static_cast<double>(n_in_bin),
                false
            );
        }
    }
    
    return 0;
}