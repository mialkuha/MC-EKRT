//Copyright (c) 2022 Mikko Kuha 

#include <algorithm>
#include <atomic>
#include <execution>
#include <csignal>
#include <iostream>
#include <mutex>
#include <random>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "ars.hpp"
#include "generic_helpers.hpp"
#include "high_level_calcs.hpp"
#include "histo.hpp"
#include "io_helpers.hpp"
#include "LHAPDF/GridPDF.h"
#include "nucleus_generator.hpp"
#include "pqcd.hpp"
#include "typedefs.hpp"

auto find_sigma_jet_cutoff
(
    momentum &kt02, 
    const momentum &mand_s, 
    const xsectval &target, 
    std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
    const pqcd::sigma_jet_params &jet_params, 
    const bool &verbose=true
) noexcept -> void
{
    xsectval sigma_jet=0.0;
    kt02 = 2.0;

    auto difference_to_target = [&](const momentum &_kt02)
    {
        return pqcd::calculate_sigma_jet(p_p_pdf, &mand_s, &_kt02, jet_params) - target;
    };

    helpers::secant_method(&kt02, difference_to_target, 1e-3, &sigma_jet);

    if (verbose) std::cout<<kt02<<' '<<sigma_jet+target<<' '<<target<<std::endl;
    
    return;
}

auto check_and_place_circle_among_others_with_overlaps
(
    std::tuple<double, double, double, uint16_t> cand_circle, 
    std::vector<std::tuple<double, double, double, uint16_t> > &final_circles, 
    const double &maximum_overlap,
    std::uniform_real_distribution<double> &unirand,
    std::shared_ptr<std::mt19937> random_generator
) noexcept -> bool
{
    auto & [cand_x, cand_y, cand_r, cand_overlap] = cand_circle;
    std::vector<uint16_t*> overlaps_with;
    
    for (auto & [circ_x, circ_y, circ_r, circ_overlap] : final_circles)
    {
        if ( pow((circ_x-cand_x),2) + pow((circ_y-cand_y),2) < pow((circ_r+cand_r),2) )
        {
            cand_overlap += 1;
            if (cand_overlap > maximum_overlap)
            {
                if
                (
                    ((cand_overlap - maximum_overlap) >= 1) || 
                    ((cand_overlap - maximum_overlap) > unirand(*random_generator))
                )
                {
                    return false;
                }
            }
            if ((circ_overlap + 1) > maximum_overlap)
            {
                if 
                (
                    (circ_overlap >= maximum_overlap) || 
                    ((circ_overlap + 1 - maximum_overlap) > unirand(*random_generator))
                )
                {
                    return false;
                }
            }
            overlaps_with.push_back(&circ_overlap);
        }
    }
    
    for (auto & c : overlaps_with)
    {
        (*c)++;
    }

    final_circles.emplace_back(std::move(cand_circle));
        
    return true;
}

auto check_and_place_circle_among_others_with_disks
(
    std::tuple<double, double, double, uint16_t> cand_circle, 
    std::vector<std::tuple<double, double, double, uint16_t> > &final_circles
) noexcept -> bool
{
    auto & [cand_x, cand_y, cand_r, cand_overlap] = cand_circle;
    std::vector<uint16_t*> overlaps_with;
    
    for (auto & [circ_x, circ_y, circ_r, circ_overlap] : final_circles)
    {
        if ( pow((circ_x-cand_x),2) + pow((circ_y-cand_y),2) < pow((circ_r+cand_r),2) )
        {
            return false;
        }
    }

    final_circles.emplace_back(std::move(cand_circle));
        
    return true;
}

struct dijet_with_ns
{
    dijet_specs dijet;
    nucleon * pro_nucleon;
    nucleon * tar_nucleon;
    dijet_with_ns(dijet_specs dijet_, nucleon * pro_nucleon_, nucleon * tar_nucleon_)
        : dijet(std::move(dijet_)), pro_nucleon(pro_nucleon_), tar_nucleon(tar_nucleon_) { }
};

auto throw_location_for_dijet
(
    const dijet_with_ns &cand,
    const spatial &proton_width,
    std::normal_distribution<double> &normal_dist, 
    std::shared_ptr<std::mt19937> random_generator
) noexcept -> std::tuple<double, double>
{
    auto param = std::normal_distribution<>::param_type{0., proton_width};
    auto dx = normal_dist(random_generator,param);
    auto dy = normal_dist(random_generator,param);

    auto retx = 0.5*(cand.pro_nucleon->co.x + cand.tar_nucleon->co.x + M_SQRT2*dx);
    auto rety = 0.5*(cand.pro_nucleon->co.y + cand.tar_nucleon->co.y + M_SQRT2*dy);

    return std::make_tuple(retx, rety);
}

//Checks whether the saturation criterium allows the candidate to be added TODO
auto filter_collisions_saturation
(
    std::vector<dijet_with_ns> &candidates, 
    std::vector<dijet_specs> &final_candidates,
    const double &maximum_overlap,
    const spatial &proton_width,
    std::uniform_real_distribution<double> &unirand,
    std::shared_ptr<std::mt19937> random_generator,
    const bool with_overlaps
) noexcept -> void
{
    std::vector<std::tuple<double, double, double, uint16_t> > final_circles;
    final_circles.reserve(candidates.size());
    final_candidates.reserve(final_candidates.size()+candidates.size());

    std::normal_distribution<double> normal_dist(0,0);

    std::sort(candidates.begin(), candidates.end(), //Sort the candidate events so that the one with the biggest kt is first
        [](dijet_with_ns &s1, dijet_with_ns &s2) { return (s1.dijet.kt > s2.dijet.kt); });

    if (!with_overlaps)
    {
        for (auto & cand : candidates)
        {
            auto [cand_x, cand_y] = throw_location_for_dijet(cand, proton_width, normal_dist, random_generator);
            auto cand_circle = std::make_tuple(cand_x, cand_y, 1/(cand.dijet.kt*maximum_overlap), 1);
            if (check_and_place_circle_among_others_with_disks(std::move(cand_circle), final_circles))
            {
                final_candidates.emplace_back(std::move(cand.dijet));
            }
        }
    }
    else
    {
        for (auto & cand : candidates)
        {
            auto [cand_x, cand_y] = throw_location_for_dijet(cand, proton_width, normal_dist, random_generator);
            auto cand_circle = std::make_tuple(cand_x, cand_y, 1/cand.dijet.kt, 1);
            if (check_and_place_circle_among_others_with_overlaps(std::move(cand_circle), final_circles, maximum_overlap, unirand, random_generator))
            {
                final_candidates.emplace_back(std::move(cand.dijet));
            }
        }
    }
    
    final_candidates.shrink_to_fit();
}

struct colls_with_ns
{
    momentum kt;
    rapidity y1;
    rapidity y2;
    nucleon * pro_nucleon;
    nucleon * tar_nucleon;
    dijet_specs * dijet;
    colls_with_ns(momentum kt_, rapidity y1_, rapidity y2_, nucleon * pro_nucleon_, nucleon * tar_nucleon_, dijet_specs * dijet_)
        : kt(kt_), y1(y1_), y2(y2_), pro_nucleon(pro_nucleon_), tar_nucleon(tar_nucleon_), dijet(dijet_) { }
};
//Empties all of the binary_collisions
auto filter_collisions_MC 
(
    std::vector<nn_coll> &binary_collisions,
    std::vector<dijet_with_ns> &final_candidates,
    const momentum sqrt_s,
    const bool deplete_nucleons
) noexcept -> void
{
    std::unordered_map<nucleon*, double> x1s;
    std::unordered_map<nucleon*, double> x2s;
    std::unordered_set<nucleon*> depleted_pro;
    std::unordered_set<nucleon*> depleted_tar;

    std::vector<colls_with_ns> collision_candidates;
    collision_candidates.reserve(binary_collisions.size()*10);

    for (auto &col : binary_collisions)
    {
        for (auto &dij : col.dijets)
        {
            collision_candidates.emplace_back(dij.kt, dij.y1, dij.y2, col.projectile, col.target, &dij);
        }
    }
    collision_candidates.shrink_to_fit();

    std::sort(collision_candidates.begin(), collision_candidates.end(), //Sort the candidates so that the one with the biggest kt is first
              [](colls_with_ns &s1, colls_with_ns &s2) { return (s1.kt > s2.kt); });

    for (auto & cand : collision_candidates)
    {
        if (deplete_nucleons && (depleted_pro.contains(cand.pro_nucleon) || depleted_tar.contains(cand.tar_nucleon)))
        {
            continue;
        }

        bool discard = false;
        auto x1 = (cand.kt / sqrt_s) * (exp(cand.y1) + exp(cand.y2));
        auto x2 = (cand.kt / sqrt_s) * (exp(-cand.y1) + exp(-cand.y2));

        auto i_x1_sum = x1s.find(cand.pro_nucleon);
        if (i_x1_sum != x1s.end())
        {
            i_x1_sum->second += x1;

            if (i_x1_sum->second > 1.0)
            {
                discard = true;
                if (deplete_nucleons)
                {
                    depleted_pro.insert(cand.pro_nucleon);
                }
            }
        }
        else
        {
            x1s.insert({cand.pro_nucleon, x1});
        }

        auto i_x2_sum = x2s.find(cand.tar_nucleon);
        if (i_x2_sum != x2s.end())
        {
            i_x2_sum->second += x2;

            if (i_x2_sum->second > 1.0)
            {
                discard = true;
                if (deplete_nucleons)
                {
                    depleted_tar.insert(cand.tar_nucleon);
                }
            }
        }
        else
        {
            x2s.insert({cand.tar_nucleon, x2});
        }

        if (!discard)
        {
            final_candidates.emplace_back(std::move(*(cand.dijet)), cand.pro_nucleon, cand.tar_nucleon);
        }
    }

    //std::cout << "candidates filtered: " << collision_candidates.size() - final_candidates.size() << " out of " << collision_candidates.size() << std::endl;

    binary_collisions.clear();
    //for (auto &col : binary_collisions)
    //{
    //    col.dijets.erase(col.dijets.begin(), col.dijets.end());
    //}
}

auto filter_end_state
(
    std::vector<nn_coll> &binary_collisions, 
    std::vector<dijet_specs> &filtered_scatterings,
    std::uniform_real_distribution<double> &unirand,
    std::shared_ptr<std::mt19937> random_generator,
    const bool mom_cons = false,
    const bool saturation = false,
    const bool deplete_nucleons = false,
    const bool saturation_with_overlap = false,
    const momentum sqrt_s = 0,
    const double &maximum_overlap = 2.0,
    const spatial &proton_width = 1.0
) noexcept -> void
{
    std::vector<dijet_with_ns> candidates;
    candidates.reserve(binary_collisions.size()*10); //10 events on average is just an overhead guess
    
    if (mom_cons) 
    {
        filter_collisions_MC(binary_collisions, candidates, sqrt_s, deplete_nucleons);
    }
    else
    {
        for (auto &col : binary_collisions) //all of the end states of the binary events are now candidate events
        {
            for (auto & dij : col.dijets)
            {
                candidates.emplace_back(std::move(dij), col.projectile, col.target);
            }
        }
        binary_collisions.clear();
    }
    candidates.shrink_to_fit();
    

    if (saturation)
    {
        filter_collisions_saturation(candidates, filtered_scatterings, maximum_overlap, proton_width, unirand, random_generator, saturation_with_overlap);
    }
    else
    {
        for (auto &candidate : candidates) //all of the candidates are passed to filtered
        {
            filtered_scatterings.push_back(candidate.dijet);
        }
    }
    
    //std::cout << "filtered_scatterings.size() =  " << filtered_scatterings.size() << std::endl;
}

volatile std::atomic_bool user_aborted = false;
volatile std::atomic_bool g_bug_bool = false;
void abort_handler(int num)
{
    user_aborted = true;
}
void term_handler(int num)
{
    std::cout<<std::endl<<"Terminate called "<<num<<std::endl;
}

int main(int argc, char** argv)
{ 
    /*bool g_is_pp;
    bool g_is_pa;
    bool g_is_aa;
    bool g_use_npdfs;
    bool g_use_snpdfs;
    bool g_is_mom_cons;
    bool g_is_mom_cons_new;
    bool g_are_ns_depleted;
    bool g_is_saturation;
    std::string name_postfix;*/

    auto 
    [
        name_postfix,
        g_is_pp,
        g_is_pa,
        g_is_aa,
        g_use_npdfs,
        g_use_snpdfs,
        g_is_mom_cons,
        g_is_mom_cons_new,
        g_are_ns_depleted,
        g_is_saturation,
        g_is_sat_overlap,
        desired_N_events,
        b_max
    ] = io::read_conf(std::string(argv[1]));

    /*switch (argv[1][0])
    {
    case '0':
        name_postfix = "_pA_2500k_mb_PDF";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = false;
        g_use_snpdfs = false;
        g_is_mom_cons = false;
        g_is_mom_cons_new = false;
        g_are_ns_depleted = false;
        g_is_saturation = false;
        break;
    case '1':
        name_postfix = "_pA_2500k_mb_nPDF";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = true;
        g_use_snpdfs = false;
        g_is_mom_cons = false;
        g_is_mom_cons_new = false;
        g_are_ns_depleted = false;
        g_is_saturation = false;
        break;
    case '2':
        name_postfix = "_pA_2500k_mb_snPDF";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = true;
        g_use_snpdfs = true;
        g_is_mom_cons = false;
        g_is_mom_cons_new = false;
        g_are_ns_depleted = false;
        g_is_saturation = false;
        break;
    case '3':
        name_postfix = "_pA_2500k_mb_PDF_MC";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = false;
        g_use_snpdfs = false;
        g_is_mom_cons = false;
        g_is_mom_cons_new = true;
        g_are_ns_depleted = false;
        g_is_saturation = false;
        break;
    case '4':
        name_postfix = "_pA_2500k_mb_nPDF_MC";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = true;
        g_use_snpdfs = false;
        g_is_mom_cons = false;
        g_is_mom_cons_new = true;
        g_are_ns_depleted = false;
        g_is_saturation = false;
        break;
    case '5':
        name_postfix = "_pA_2500k_mb_snPDF_MC";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = true;
        g_use_snpdfs = true;
        g_is_mom_cons = false;
        g_is_mom_cons_new = true;
        g_are_ns_depleted = false;
        g_is_saturation = false;
        break;
    case '6':
        name_postfix = "_pA_2500k_mb_PDF_MC_ND";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = false;
        g_use_snpdfs = false;
        g_is_mom_cons = false;
        g_is_mom_cons_new = true;
        g_are_ns_depleted = true;
        g_is_saturation = false;
        break;
    case '7':
        name_postfix = "_pA_2500k_mb_nPDF_MC_ND";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = true;
        g_use_snpdfs = false;
        g_is_mom_cons = false;
        g_is_mom_cons_new = true;
        g_are_ns_depleted = true;
        g_is_saturation = false;
        break;
    case '8':
        name_postfix = "_pA_2500k_mb_snPDF_MC_ND";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = true;
        g_use_snpdfs = true;
        g_is_mom_cons = false;
        g_is_mom_cons_new = true;
        g_are_ns_depleted = true;
        g_is_saturation = false;
        break;
    default:
        break;
    }*/

    //A lot of printing
    bool verbose = false;

    if (verbose) std::cout<<"Initializing..."<<std::flush;
    
    //General parameters for the simulation
    const bool    read_nuclei_from_file    = false, 
                  read_events_from_file    = false,
                  read_sigmajets_from_file = true,
                  end_state_filtering      = true, 
                  save_events              = false/*, 
                  average_spatial_taas     = false*/;
    std::string   event_file_name = "event_log_"+name_postfix+".dat";
    uint32_t      AA_events             = 0;
    const spatial b_min                 = 0;
    std::mutex AA_events_mutex; 
    std::cout<<std::endl<<std::endl<<"Doing the run "<<name_postfix<<std::endl;
    //auto eng = std::make_shared<std::mt19937>(static_cast<ulong>(1));
    auto eng = std::make_shared<std::mt19937>(static_cast<ulong>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<double> unirand{0.0, 1.0};

    //Parameters for the nuclei
    //uint8_t samplersN = 32;
    const spatial rad_min=0,
                  rad_max=30;
    const std::function<double(const double&)> rad_pdf{[](const double & x)
    {
        return x*x/(1+exp((x-6.624)/0.549));
    }};
    auto radial_sampler{std::make_shared<ars>(rad_pdf, rad_min, rad_max)};
    do //while (radial_sampler->is_adaptive())
    {
        radial_sampler->throw_one(*eng);
    } while (radial_sampler->is_adaptive());
    //std::vector<std::shared_ptr<ars> > sampler_pointers;
    //std::mutex radial_sampler_mutex;
    //for (uint8_t i=0; i<samplersN; i++)
    //{
    //    sampler_pointers.push_back(std::make_shared<ars>(rad_pdf, rad_min, rad_max));
    //}

    nucleus_generator::nucleus_params nuc_params = 
    {
        /* .NA=                   */(g_is_aa)? 208u : 1u, //Pb 
        /* .ZA=                   */(g_is_aa)? 82u : 1u, //Pb 
        /* .NB=                   */(g_is_pp)? 1u : 208u, //Pb 
        /* .ZB=                   */(g_is_pp)? 1u : 82u, //Pb 
        /* .min_distance=         */0.4, 
        /* .shift_cms=            */true, 
        /* .correct_overlap_bias= */true
    };
    
    //Parameters for the hard collisions
    const spatial proton_width = 0.573;
    const spatial proton_width_2 = pow(0.573, 2);
    const std::function<spatial(const spatial&)> Tpp{[&proton_width_2](const spatial &bsquared)
    {
        return exp(-bsquared / (4 * proton_width_2)) / (40 * M_PI * proton_width_2); // 1/fm² = mb/fm² * 1/mb = 0.1 * 1/mb
    }}; 
    const xsectval sigma_inel_for_glauber = 41.5;//mb
    const momentum sqrt_s                 = 5020;//GeV
    const momentum mand_s                 = pow(sqrt_s, 2);//GeV^2
    momentum kt0                          = 2.728321;//GeV
    momentum kt02                         = pow(kt0, 2);//GeV^2
    rapidity ylim                         = static_cast<rapidity>(log(sqrt_s / kt0));
    auto p_pdf = std::make_shared<LHAPDF::GridPDF>("CT14lo", 0);
    const AA_collision_params coll_params
    {
    /*mc_glauber_mode=          */false,
    /*pp_scattering=            */g_is_pp,
    /*pA_scattering=            */g_is_pa,
    /*spatial_pdfs=             */g_use_snpdfs,
    /*spatial_averaging=        */false,
    /*calculate_end_state=      */true,
    /*reduce_nucleon_energies=  */g_is_mom_cons,
    /*sigma_inel_for_glauber=   */sigma_inel_for_glauber,
    /*Tpp=                      */Tpp,
    /*normalize_to=             */B2_normalization_mode::inelastic,
    /*sqrt_s=                   */sqrt_s,
    /*energy_threshold=         */kt0
    };
    //std::vector<Coll> collisions_for_reporting;
    std::vector<double> kt_bins{helpers::loglinspace(kt0, 100, 16)};
    auto kt_v_dummy = helpers::loglinspace(200, sqrt_s/2.0, 5);
    for (auto ktdummy : kt_v_dummy)
    {
        kt_bins.push_back(ktdummy);
    }
    const std::vector<double> y_bins{helpers::linspace(-ylim, ylim, 40)};
    const std::vector<double> b_bins{helpers::linspace(b_min, b_max, 21)};

    //sigma_jet parameters
    /*const bool read_sigmajets_from_file = false;*/
    pqcd::diff_sigma::params diff_params = pqcd::diff_sigma::params(
    /*projectile_with_npdfs=    */(g_is_aa && g_use_npdfs),
    /*target_with_npdfs=        */(!g_is_pp && g_use_npdfs),
    /*isoscalar_projectile=     */false,
    /*isoscalar_target=         */false,
    /*npdfs_spatial=            */coll_params.spatial_pdfs,
    /*npdf_setnumber=           */1,
    /*A=                        */(g_is_aa)? 208 : 1, //Pb 
    /*B=                        */(g_is_pp)? 1 : 208, //Pb
    /*ZA=                       */(g_is_aa)? 82 : 1,  //Pb
    /*ZB=                       */(g_is_pp)? 1 : 82   //Pb
    /*p_n_pdf=                  */
    /*rA_spatial=               */
    /*rB_spatial=               */);
    pqcd::sigma_jet_params jet_params = pqcd::sigma_jet_params(
    /*d_params=                 */diff_params,
    /*scale_choice=             */pqcd::scaled_from_kt,
    /*scalar=                   */1.0,
    /*use_ses=                  */false);

    //calculate_sigma_1jet_analytical(mand_s, kt_bins, y_bins, p_pdf, &jet_params, name_postfix);
    //calculate_sigma_jet_analytical(mand_s, kt_bins, y_bins, p_pdf, &jet_params);
    //calculate_sigma_dijet_analytical(mand_s, kt_bins, y_bins, p_pdf, &jet_params, name_postfix);
//std::vector<double> mands = {2., 5., 10., 20., 50.,100.,200.,500.,1000.,2000.,5000.,};
//for (auto & m : mands) m = m*m;
//std::vector<double> kt02ss = {1., 1., 1., 1., 1.,1.,1.,1.,1.,1.,1.,};
//std::vector<double> inels = {1., 1., 1., 1., 1.,1.,1.,1.,1.,1.,1.,};
//for (int i=0; i<11 ; i++) find_sigma_jet_cutoff(kt02ss[i], mands[i], 124.6635, p_pdf, jet_params, true);
//std::cout<<kt02<<std::endl;
//return 0;
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
        (diff_params.projectile_with_npdfs || diff_params.target_with_npdfs),
        coll_params.spatial_pdfs,
        coll_params.reduce_nucleon_energies,
        read_sigmajets_from_file,
        p_pdf, 
        mand_s,
        sqrt_s,
        kt02, 
        kt0, 
        jet_params
    );
    
    //find_sigma_jet_cutoff(kt02, mand_s, 124.6635, p_pdf, jet_params, true);
    //std::cout<<kt02<<std::endl;
    
    //if (!read_sigmajets_from_file)
    //{
        //double tolerance=0.05, upper_tAA_0_limit=46.0, lower_tAA_0_limit = 30.0, upper_sumTpp_limit=0.61, lower_sumTpp_limit=0.03411;
      //  double tolerance=0.05, upper_tAA_0_limit=42.0, lower_tAA_0_limit = 33.0, upper_sumTpp_limit=0.5, lower_sumTpp_limit=0.03411;
      //  std::variant<InterpMultilinear<4, xsectval>, linear_interpolator, xsectval> 
      //      sigma_jets = calculate_spatial_sigma_jets_mf(tolerance, p_pdf, mand_s, kt02, jet_params, upper_tAA_0_limit, lower_tAA_0_limit, upper_sumTpp_limit, lower_sumTpp_limit);
        ////array<spatial,4> args{0.1, 0.3, 35.0, 40.0};
        ////std::cout<<sigma_jets.interp(args.begin())<<std::endl;
        ////std::cout<<pqcd::calculate_spatial_sigma_jet_mf(p_pdf, p_pdf, &mand_s, &kt02, &jet_params, &args[0], &args[1], &args[2], &args[3])<<std::endl;

        //double tolerance=0.25;//, upper_tAA_0_limit=46.0, lower_tAA_0_limit = 30.0, upper_sumTpp_limit=0.61, lower_sumTpp_limit=0.03411;
        //                                            1.35134e-186, 1.4788e-177, 1.02829e-215       4.24149, 3.609, 10.7237
        //const std::array<const double, 4>  lower_limits = {0.0, 0.0, 0.0, 30.0}, upper_limits = {4.5, 4.5, 11.0, 46.0};
        //const std::array<const double, 4>  lower_limits = {0.0, 0.0, 0.0, 28.4}, upper_limits = {3.6, 3.6, 8.8, 41.2};
        //const std::array<const double, 4>  lower_limits = {0.0, 0.0, 0.0, 30.5}, upper_limits = {3.6, 3.6, 8.8, 30.5};
        //double tolerance=0.05, upper_tAA_0_limit=42.0, lower_tAA_0_limit = 33.0, upper_sumTpp_limit=0.5, lower_sumTpp_limit=0.03411;
        //InterpMultilinear<5, xsectval> sigma_jets = calculate_spatial_sigma_jets_full(tolerance, p_pdf, p_pdf, mand_s, kt02, jet_params, lower_limits, upper_limits);
        //array<spatial,5> args{3.5, 4.3, 8.0, 39.0, 42.0}; //-304.782
        //std::cout<<sigma_jets.interp(args.begin())<<std::endl;

        //std::cout<<pqcd::calculate_spatial_sigma_jet_mf(p_pdf, p_pdf, &mand_s, &kt02, &jet_params, &args[0], &args[1], &args[2], &args[3])<<std::endl;
    //}    

    if (verbose) std::cout<<"Done!"<<std::endl;
    
    //std::ofstream log_file;
    //log_file.open("log2.dat", std::ios::out);

    histo_2d jets{kt_bins, y_bins};
    histo_2d dijets{kt_bins, y_bins};
    std::mutex jets_mutex; 
    std::mutex dijets_mutex;

    histo_1d dETdy {y_bins};
    histo_1d dEdy  {y_bins};
    histo_1d dNdy  {y_bins};
    std::mutex dETdy_mutex; 
    std::mutex dEdy_mutex;
    std::mutex dNdy_mutex;

    histo_1d dETdeta {y_bins};
    histo_1d dEdeta  {y_bins};
    std::mutex dETdeta_mutex; 
    std::mutex dEdeta_mutex;

    std::ofstream total_energy;
    std::mutex total_energy_mutex; 
    auto max_energy = 0.5*(diff_params.A+diff_params.B)*sqrt_s;
    int64_t max_energy_broken = 0;
 
    total_energy.open("total_energies_"+name_postfix+".dat");
    total_energy << "///Sum E_T Sum E" << std::endl;
 
 
    histo_1d dETdb {b_bins};
    histo_1d dEdb  {b_bins};
    std::mutex dETdb_mutex; 
    std::mutex dEdb_mutex;

    if (!read_events_from_file)
    {
        std::ofstream event_file;
        std::mutex event_file_mutex; 

        if (save_events)
        {
            event_file.open(event_file_name);

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
            event_file << "///Spatial npdfs: "<<coll_params.spatial_pdfs<<std::endl;
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

        std::vector<uint64_t> event_indexes(desired_N_events);
        std::iota(event_indexes.begin(), event_indexes.end(), 0); //generates the list as {0,1,2,3,...}
        std::atomic<uint64_t> running_count{desired_N_events};

        std::signal(SIGINT, abort_handler);
        //std::signal(SIGTERM, term_handler);
        //std::signal(SIGSEGV, term_handler);

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
                 &g_use_snpdfs=g_use_snpdfs,
                 &g_is_mom_cons_new=g_is_mom_cons_new,
                 &g_is_saturation=g_is_saturation,
                 &g_is_sat_overlap=g_is_sat_overlap,
                 &g_are_ns_depleted=g_are_ns_depleted,
                 &sigma_jet=sigma_jet,
                 &power_law=power_law,
                 &envelope_maximum=envelope_maximum
                ](const uint64_t index) 
                {
                    do //while (g_bug_bool)
                    {
                        uint32_t NColl = 0;
                        std::vector<nn_coll> binary_collisions;
                        std::vector<dijet_specs> filtered_scatterings;

                        //B^2 from a uniform distribution
                        spatial impact_parameter = sqrt(b_min*b_min + unirand(*eng)*(b_max*b_max-b_min*b_min)); 

                        //std::shared_ptr<ars> radial_sampler;
                        //bool sampler_found = false;
                        //uint8_t s_index = index % samplersN;

                        //{
                        //    std::unique_lock<std::mutex> lock_rad(radial_sampler_mutex);
                        //    do //while sampler_found
                        //    {
                        //        if (!sampler_pointers[s_index]->locked())
                        //        {
                        //            radial_sampler = sampler_pointers[s_index];
                        //            sampler_found = true;
                        //        }
                        //        s_index = (s_index + 1) % samplersN;
                        //        //std::this_thread::sleep_for(std::chrono::milliseconds(10));
                        //    }
                        //    while (!sampler_found);
                        //    //std::cout<<std::endl<<s_index<<std::endl;
                        //}
                        //std::unique_lock<std::mutex> lock_rad(radial_sampler_mutex);
                        //auto radial_sampler = sampler_pointers[0];
                        std::vector<nucleon> pro, tar;

                        bool bugged;
                        do //while (!bugged)
                        {
                            bugged = false;
                            try
                            {
                                auto [pro_dummy, tar_dummy] = calcs::generate_nuclei
                                (
                                    nuc_params, 
                                    sqrt_s, 
                                    impact_parameter, 
                                    eng, 
                                    radial_sampler, 
                                    read_nuclei_from_file, 
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
                        } while (bugged);
                        
                        do
                        {
                            binary_collisions.erase(binary_collisions.begin(), binary_collisions.end());
                            if (verbose) std::cout<<"impact_parameter: "<<impact_parameter<<std::endl;

                            calcs::collide_nuclei
                            (
                                g_use_snpdfs,
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

                            /*if (average_spatial_taas)
                            {
                                //collide_nuclei_with_spatial_pdfs_averaging(pro, tar, binary_collisions, sigma_jets, unirand, eng, coll_params, verbose, Tpp);
                            }
                            else
                            {
                                //log_file << impact_parameter << std::endl;
                                //collide_nuclei_with_spatial_pdfs_factored(pro, tar, binary_collisions, sigma_jets, unirand, eng, coll_params, verbose, Tpp, proton_width_2, {impact_parameter,0,0}, log_file);
                                //collide_nuclei_with_spatial_pdfs_full(pro, tar, binary_collisions, sigma_jets, unirand, eng, coll_params, verbose, Tpp, proton_width_2, {impact_parameter,0,0}, log_file);
                                //log_file << std::endl;
                            }*/
                            NColl = static_cast<uint32_t>(binary_collisions.size());
                            if (NColl<1)
                            {
                                impact_parameter = sqrt(b_min*b_min + unirand(*eng)*(b_max*b_max-b_min*b_min));
                                //const std::lock_guard<std::mutex> lock(radial_sampler_mutex);
                                std::tie(pro, tar) = calcs::generate_nuclei
                                (
                                    nuc_params, 
                                    sqrt_s, 
                                    impact_parameter, 
                                    eng, 
                                    radial_sampler, 
                                    read_nuclei_from_file, 
                                    verbose
                                );
                            }
                        } while (NColl<1);

                        if(!end_state_filtering && save_events)
                        {
                            const std::lock_guard<std::mutex> lock(event_file_mutex);
                            io::save_event(event_file, pro, tar, impact_parameter);
                        }
                        else if (end_state_filtering)
                        {
                            momentum ET=0, E=0;
                            filter_end_state
                            (
                                binary_collisions, 
                                filtered_scatterings, 
                                unirand, 
                                eng, 
                                g_is_mom_cons_new, 
                                g_is_saturation, 
                                g_are_ns_depleted, 
                                g_is_sat_overlap, 
                                sqrt_s, 
                                proton_width
                            );

                            std::vector<std::tuple<double, double> > new_jets;
                            std::vector<std::tuple<double, double> > new_dijets;
                            std::vector<std::tuple<double, double> > new_ET_y;
                            std::vector<std::tuple<double, double> > new_E_y;
                            std::vector<std::tuple<double, double> > new_N_y;
                            std::vector<std::tuple<double, double> > new_ET_eta;
                            std::vector<std::tuple<double, double> > new_E_eta;

                            for (auto e : filtered_scatterings)
                            {
                                new_jets.emplace_back(e.kt, e.y1);
                                new_jets.emplace_back(e.kt, e.y2);
                                
                                new_dijets.emplace_back(e.kt, 0.5*(e.y1+e.y2));

                                new_ET_y.emplace_back(e.y1, e.kt);
                                new_ET_y.emplace_back(e.y2, e.kt);
                                
                                new_ET_eta.emplace_back(0.5*(e.y1+e.y2), 2*e.kt);
                                
                                new_E_y.emplace_back(e.y1, e.kt*cosh(e.y1));
                                new_E_y.emplace_back(e.y2, e.kt*cosh(e.y2));
                                
                                new_N_y.emplace_back(e.y1, 1);
                                new_N_y.emplace_back(e.y2, 1);
                                
                                new_E_eta.emplace_back(0.5*(e.y1+e.y2), e.kt*(cosh(e.y1) + cosh(e.y2)));
                                
                                ET += 2*e.kt;
                                E += e.kt*(cosh(e.y1) + cosh(e.y2));
                            }

                            {
                                const std::lock_guard<std::mutex> lock(jets_mutex);
                                jets.add(new_jets);
                            }

                            {
                                const std::lock_guard<std::mutex> lock(dijets_mutex);
                                dijets.add(new_dijets);
                            }

                            {
                                const std::lock_guard<std::mutex> lock(dETdy_mutex);
                                dETdy.add(new_ET_y);
                            }

                            {
                                const std::lock_guard<std::mutex> lock(dEdy_mutex);
                                dEdy.add(new_E_y);
                            }

                            {
                                const std::lock_guard<std::mutex> lock(dNdy_mutex);
                                dNdy.add(new_N_y);
                            }

                            {
                                const std::lock_guard<std::mutex> lock(dETdeta_mutex);                    
                                dETdeta.add(new_ET_eta);
                            }

                            {
                                const std::lock_guard<std::mutex> lock(dEdeta_mutex);
                                dEdeta.add(new_E_eta);
                            }
                            
                            double sum_ET = 0, sum_E = 0;
                            for (auto e : new_ET_y)
                            {
                                sum_ET += std::get<1>(e);
                            }
                            for (auto e : new_E_y)
                            {
                                sum_E += std::get<1>(e);
                            }

                            {
                                const std::lock_guard<std::mutex> lock(total_energy_mutex);
                                total_energy << sum_ET << ' ' << sum_E << std::endl;
                            }
                            
                            if (sum_E > max_energy)
                            {
                                std::cout << std::endl 
                                        << "Energy conservation violated! Total "<< ++max_energy_broken
                                        << " times this far." << std::endl
                                        << "E_T = " << sum_ET << ", E = " << sum_E  << std::endl;
                            }
                            
                            {
                                const std::lock_guard<std::mutex> lock(dETdb_mutex);
                                dETdb.add(std::make_tuple(impact_parameter, ET));
                            }
                            {
                                const std::lock_guard<std::mutex> lock(dEdb_mutex);
                                dEdb.add(std::make_tuple(impact_parameter, E));
                            }

                            if (save_events)
                            {
                                const std::lock_guard<std::mutex> lock(event_file_mutex);
                                io::save_event(event_file, pro, tar, impact_parameter, filtered_scatterings);
                            }

                            //filtered_scatterings.erase(filtered_scatterings.begin(), filtered_scatterings.end());
                        }
                        else
                        {
                            //binary_collisions.erase(binary_collisions.begin(), binary_collisions.end());
                        }

                        /*if (verbose || (nof_collisions%100)==0)
                        {
                            std::cout << std::endl << "A+A collided thus far: " << nof_collisions << ", of which events thus far: " << AA_events << std::endl << std::endl;
                        }*/

                        {
                            const std::lock_guard<std::mutex> lock(AA_events_mutex);
                            AA_events++;
                            if (AA_events % 5000 == 0 )
                            {
                                std::cout <<"\rA+A collisions calculated: " << AA_events << std::flush;
                            }
                        }
                        //PrintThread{} <<"\rA+A collisions calculated: " << AA_events << std::flush;

                        /*uint32_t Npart=0;
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

                        Coll coll(NColl, Npart, impact_parameter, eng);
                        collisions_for_reporting.push_back(coll);*/
                        
                        //filtered_scatterings.clear();
                
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
    }/*
    else
    {
        std::ifstream event_file;
        event_file.open(event_file_name);

        if (event_file.is_open())
        {
            std::string line;
            std::getline(event_file, line);
            //Skip comments at the start
            while (line[0] == '/')
            {
                std::getline(event_file, line);
            }

            do
            {
                std::getline(event_file, line);
                std::istringstream line_stream(line);
                spatial impact_parameter;
                line_stream >> impact_parameter;

                std::vector<nucleon> pro;
                std::getline(event_file, line);
                
                //Read projectile
                do 
                {
                    std::getline(event_file, line);
                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    coords co;
                    double x, y, z;

                    line_stream.ignore(256,'{');
                    line_stream >> x;
                    line_stream.ignore(256,',');
                    line_stream >> y;
                    line_stream.ignore(256,',');
                    line_stream >> z;
                    co = coords{x, y, z};

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    momentum mom;
                    line_stream >> mom;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    bool wounded;
                    line_stream >> wounded;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    bool is_neutron;
                    line_stream >> is_neutron;

                    pro.push_back(nucleon{co, mom, wounded, is_neutron});

                    std::getline(event_file, line);
                } while (line.back() != '}');
                pro.shrink_to_fit();

                std::getline(event_file, line);
                std::vector<nucleon> tar;
                std::getline(event_file, line);

                //Read target
                do
                {
                    std::getline(event_file, line);
                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    coords co;
                    double x, y, z;

                    line_stream.ignore(256,'{');
                    line_stream >> x;
                    line_stream.ignore(256,',');
                    line_stream >> y;
                    line_stream.ignore(256,',');
                    line_stream >> z;
                    co = coords{x, y, z};

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    momentum mom;
                    line_stream >> mom;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    bool wounded;
                    line_stream >> wounded;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    bool is_neutron;
                    line_stream >> is_neutron;

                    tar.push_back(nucleon{co, mom, wounded, is_neutron});

                    std::getline(event_file, line);
                } while (line.back() != '}');
                tar.shrink_to_fit();
                
                std::getline(event_file, line);
                std::getline(event_file, line);

                //Read filtered scatterings
                do
                {
                    std::getline(event_file, line);
                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    momentum kt;
                    line_stream >> kt;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    rapidity y1;
                    line_stream >> y1;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    rapidity y2;
                    line_stream >> y2;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    particle_id init1;
                    line_stream >> init1;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    particle_id init2;
                    line_stream >> init2;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    particle_id final1;
                    line_stream >> final1;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    particle_id final2;
                    line_stream >> final2;

                    filtered_scatterings.push_back(dijet_specs{kt, y1, y2, init1, init2, final1, final2});

                    std::getline(event_file, line);
                } while (line.back() != '}');
                filtered_scatterings.shrink_to_fit();

                std::getline(event_file, line);
                std::getline(event_file, line);

                nof_collisions++;
                AA_events++;
                momentum ET=0, E=0;
                std::vector<std::tuple<double, double> > new_jets;
                std::vector<std::tuple<double, double> > new_dijets;
                std::vector<std::tuple<double, double> > new_ET_y;
                std::vector<std::tuple<double, double> > new_E_y;
                std::vector<std::tuple<double, double> > new_ET_eta;
                std::vector<std::tuple<double, double> > new_E_eta;

                for (auto e : filtered_scatterings)
                {
                    new_jets.emplace_back(e.kt, e.y1);
                    new_jets.emplace_back(e.kt, e.y2);

                    new_dijets.emplace_back(e.kt, 0.5*(e.y1+e.y2));

                    new_ET_y.emplace_back(e.y1, e.kt);
                    new_ET_y.emplace_back(e.y2, e.kt);

                    new_ET_eta.emplace_back(0.5*(e.y1+e.y2), 2*e.kt);

                    new_E_y.emplace_back(e.y1, e.kt*cosh(e.y1));
                    new_E_y.emplace_back(e.y2, e.kt*cosh(e.y2));

                    new_E_eta.emplace_back(0.5*(e.y1+e.y2), e.kt*(cosh(e.y1) + cosh(e.y2)));

                    ET += 2*e.kt;
                    E += e.kt*(cosh(e.y1) + cosh(e.y2));
                }
 //               jets.add(new_jets);
 //               dijets.add(new_dijets);

                //dETdy.add(new_ET_y);
                //dEdy.add(new_E_y);

                //dETdeta.add(new_ET_eta);
                //dEdeta.add(new_E_eta);

 //               dETdb.add(std::make_tuple(impact_parameter, ET));
 //               dEdb.add(std::make_tuple(impact_parameter, E));

                filtered_scatterings.erase(filtered_scatterings.begin(), filtered_scatterings.end());

                std::cout<<"\rA+A collisions read thus far: " << nof_collisions << std::flush;

                std::getline(event_file, line);
            } while (line[0] == '{');
        }
        else
        {
            std::cout<<"ERROR READING EVENTS"<<std::endl;
        }
    }*/

    io::print_2d_histo
    (
        jets, 
        "sigma1jet_sim_"+name_postfix+".dat", 
        2.0 * dijet_norm
    );

    io::print_2d_histo
    (
        jets, 
        "dNdpTdy_sim_"+name_postfix+".dat", 
        1.0,
        false
    );

    io::print_2d_histo
    (
        dijets, 
        "sigmadijet_sim_"+name_postfix+".dat", 
        dijet_norm
    );

    io::print_1d_histo
    (
        dETdy, 
        "dETdy_sim_"+name_postfix+".dat", 
        1.0 / AA_events,
        false
    );
    
    io::print_1d_histo
    (
        dEdy, 
        "dEdy_sim_"+name_postfix+".dat", 
        1.0 / AA_events,
        false
    );
    
    io::print_1d_histo
    (
        dNdy, 
        "dNdy_sim_"+name_postfix+".dat", 
        1.0 / AA_events,
        false
    );
    
    io::print_1d_histo
    (
        dETdeta, 
        "dETdeta_sim_"+name_postfix+".dat", 
        1.0 / AA_events,
        false
    );
    
    io::print_1d_histo
    (
        dEdeta, 
        "dEdeta_sim_"+name_postfix+".dat", 
        1.0 / AA_events,
        false
    );

    io::print_1d_histo
    (
        dETdb, 
        "dETdb_sim_"+name_postfix+".dat", 
        1.0,
        true
    );
    
    io::print_1d_histo
    (
        dEdb, 
        "dEdb_sim_"+name_postfix+".dat",  
        1.0,
        true
    );

/*
    std::cout << collisions_for_reporting.size() << " collisions generated" << std::endl;

    //log_file.close();

    uint nBins = 16;
    double binsLow[] = {0., 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 0.0, 0.0};
    double binsHigh[] = {0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 0.8, 1.0};
    mc_glauber_style_report(collisions_for_reporting, sigma_inel_for_glauber, desired_N_events, nBins, binsLow, binsHigh);*/

    

    return 0;
}
