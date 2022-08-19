//Copyright (c) 2022 Mikko Kuha
//TODO parallellize A+A colls
//TODO sigma_jet with full spatial nPDF calculation
//TODO sigma_inel calculation???
//TODO kuinka hyvin pätee ave(sum(T_pp(si-si~)))==T_AA(0) ??
//TODO voidaanko laskea <T_AA> ilman sigma_inel ??
//TODO <T_AA> = <<sigma_inel>_event> * <N_bin> 

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <cstring>
#include <execution>
#include <fstream>
#include <functional>
#include <future>
#include <gsl/gsl_multimin.h>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <sstream>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <vector>

#include "ars.hpp"
#include "histo.hpp"
#include "LHAPDF/GridPDF.h"
#include "linear_interpolator.hpp"
#include "linterp.h"
#include "nn_coll.hpp"
#include "nucleon.hpp"
#include "nucleus_generator.hpp"
#include "pqcd.hpp"
#include "typedefs.hpp"

using variant_sigma_jet = std::variant<InterpMultilinear<3, xsectval>, InterpMultilinear<2, xsectval>, linear_interpolator, xsectval>;
using variant_envelope_max = std::variant<linear_interpolator, xsectval>;

class PrintThread: public std::ostringstream
{
public:
    PrintThread() = default;

    ~PrintThread()
    {
        std::lock_guard<std::mutex> guard(_mutexPrint);
        std::cout << this->str();
    }

private:
    static std::mutex _mutexPrint;
};

std::mutex PrintThread::_mutexPrint{};

class Coll {

private:
  double b, x, nTwo, sumET;
  ulong nColl, nPart;
  bool alice = true;

public:
  Coll() { b = 0.; x = 0.; nTwo = 0.; nColl = 0; nPart = 0; sumET = 0.;}
  Coll(ulong nCollIn, ulong nPartIn, double bIn, std::shared_ptr<std::mt19937> random_generator){
    // ATLAS definitions.
    if (!alice) {
      x = 0.09; b = bIn; nColl = nCollIn; nPart = nPartIn; sumET = 0.;
      nTwo = 0.5 * (1. - x) * static_cast<double>(nPart) + x * static_cast<double>(nColl); sampleET(random_generator);
    // ALICE definitions.
  } else {
      x = 0.801; b = bIn; nColl = nCollIn; nPart = nPartIn; sumET = 0.;
      nTwo = x * static_cast<double>(nPart) + (1. - x) * static_cast<double>(nColl); sampleET(random_generator);
  }
}
  double getB() const {return b;}
  ulong getNcoll() const {return nColl;}
  ulong getNpart() const {return nPart;}
  double getNtwo() const {return nTwo;}
  double getET() const {return sumET;}
  void sampleET(std::shared_ptr<std::mt19937> random_generator) {
    double mu = 46.4;
    double k = 1.5;
    double p = (mu/k) / (mu/k + 1.);
    uint nAncestors = static_cast<uint>(floor(nTwo));
    uint mult = 0;
       //std::cout << (*random_generator)() << ' ' << std::flush;
    std::negative_binomial_distribution<uint> distr(static_cast<uint>(k),1.-p);
    for (uint i = 0; i < nAncestors; ++i) {
      mult += distr(*random_generator);
      //std::cout << mult << ' ' << std::flush;
    }
    sumET = mult;
    //std::cout << std::endl << nAncestors << ' ' << sumET/nAncestors << std::endl << std::endl;
    // cout << endl << nAncestors << ' ' << sumET/nAncestors << endl << endl;
  }
};
bool compNpart(Coll c1, Coll c2){
  return c1.getNpart() < c2.getNpart();
}
bool compNcoll(Coll c1, Coll c2){
  return c1.getNcoll() < c2.getNcoll();
}
bool compNtwo(Coll c1, Coll c2){
  return c1.getNtwo() < c2.getNtwo();
}
bool compET(Coll c1, Coll c2){
  return c1.getET() < c2.getET();
}
template <typename T>
double calc_ave(std::vector<Coll> &collisions, T (Coll::*func)()const ){
  std::vector<Coll>::iterator it;
  T tot = 0;
  for (it = collisions.begin(); it != collisions.end(); ++it){
    tot += ((*it).*func)();
  }
  return static_cast<double>(tot)/static_cast<double>(collisions.size());  
}

//Finds a zero of the function f
template<typename F, typename Ret_type, typename Arg_type>
auto secant_method
(
    Arg_type *const x, 
    F f, 
    const Ret_type error_tolerance, 
    Ret_type *const last_fx
) -> void
{
    Arg_type x_n = *x, x_np1 = 2.0*x_n, x_np2 = 0.0;
    Ret_type fx_n = f(x_n);

    if (abs(fx_n) < error_tolerance)
    {
        *last_fx = fx_n;
        return;
    }

    Ret_type fx_np1 = f(x_np1);

    while ( abs(fx_np1) > error_tolerance )
    {
        x_np2 = (x_n*fx_np1 - x_np1*fx_n)/(fx_np1 - fx_n);

        if (x_np2 <= 0.0) x_np2=1.0/std::numeric_limits<Arg_type>::max();

        fx_n = fx_np1;
        x_n = x_np1;
        x_np1 = x_np2;

        fx_np1 = f(x_np1);

        if (x_n == x_np1)
        {
            std::cout<<"Doesn't converge!!!"<<std::endl;
            std::cout<<x_n<<' '<<x_np1<<' '<<fx_np1<<std::endl;
            return;
        }
    }

    *last_fx = fx_np1;
    *x = x_np1;
    return;
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
double fdsigma
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

auto find_max_dsigma
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

    secant_method(&kt02, difference_to_target, 1e-3, &sigma_jet);

    if (verbose) std::cout<<kt02<<' '<<sigma_jet+target<<' '<<target<<std::endl;
    
    return;
}

auto read_nucleon_configs_from_file() noexcept
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

auto generate_nuclei
(
    const nucleus_generator::nucleus_params &nuc_params,
    const momentum &sqrt_s,
    const spatial &impact_parameter,
    std::shared_ptr<std::mt19937> eng,
    std::shared_ptr<ars> radial_sampler,
    const bool &read_nuclei_from_file,
    const bool &verbose
) noexcept
{
    std::vector<nucleon> pro, tar;
    if (read_nuclei_from_file)
    {
        if (verbose) std::cout<<"Reading nuclei..."<<std::flush;
        auto [pro_coords, tar_coords] = read_nucleon_configs_from_file();
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
        pro = nucleus_generator::generate_nucleus
              (
                  nuc_params, 
                  sqrt_s/2.0, 
                  -impact_parameter/2., 
                  eng, 
                  radial_sampler
              );
        tar = nucleus_generator::generate_nucleus
              (
                  nuc_params, 
                  sqrt_s/2.0, 
                  impact_parameter/2., 
                  eng, 
                  radial_sampler
              );
        if (verbose) std::cout<<"Done!"<<std::endl;
    }

    //if (verbose) std::cout<<"Shuffling nuclei..."<<std::flush;
    //std::shuffle(pro.begin(), pro.end(), *eng);
    //std::shuffle(tar.begin(), tar.end(), *eng);
    //if (verbose) std::cout<<"Done!"<<std::endl;

    return std::make_tuple(pro, tar);
}

//Checks whether the saturation criterium allows the candidate to be added TODO
auto fits_into_ps
(
    const std::vector<dijet_specs> &olds, 
    const dijet_specs &candidate
) noexcept -> bool
{
    if (olds.size() + static_cast<uint>(candidate.kt) < 0)
    {
        std::cout << "very strange" << std::endl;
    }
    return true;
}

struct colls_with_ns
{
    momentum kt;
    rapidity y1;
    rapidity y2;
    nucleon * pro_nucleon;
    nucleon * tar_nucleon;
    dijet_specs * dijet;
};

auto filter_collisions_MC //TODO discuss with Kari
(
    std::vector<nn_coll> &binary_collisions, 
    std::vector<dijet_specs> &final_candidates, 
    const momentum sqrt_s
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
            collision_candidates.push_back(colls_with_ns({dij.kt, dij.y1, dij.y2, col.projectile, col.target, &dij}));
        }
    }

    std::sort(collision_candidates.begin(), collision_candidates.end(), //Sort the candidates so that the one with the biggest kt is first
              [](colls_with_ns &s1, colls_with_ns &s2) { return (s1.kt > s2.kt); });

    for (auto & cand : collision_candidates)
    {
        if (depleted_pro.contains(cand.pro_nucleon) || depleted_tar.contains(cand.tar_nucleon))
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
                depleted_pro.insert(cand.pro_nucleon);
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
                depleted_tar.insert(cand.tar_nucleon);
            }
        }
        else
        {
            x2s.insert({cand.tar_nucleon, x2});
        }

        if (!discard)
        {
            final_candidates.push_back(*(cand.dijet));
        }
    }

    //std::cout << "candidates filtered: " << collision_candidates.size() - final_candidates.size() << " out of " << collision_candidates.size() << std::endl;

    //col.dijets.erase ?
}

auto filter_end_state
(
    std::vector<nn_coll> &binary_collisions, 
    std::vector<dijet_specs> &filtered_scatterings,
    const bool mom_cons = false,
    const momentum sqrt_s = 0
) noexcept -> void
{
    std::vector<dijet_specs> candidates;
    candidates.reserve(binary_collisions.size()*10); //10 events on average is just an overhead guess

    if (mom_cons) //TODO discuss with Kari
    {
        filter_collisions_MC(binary_collisions, candidates, sqrt_s);
    }
    else
    {
        for (auto &col : binary_collisions) //all of the end states of the binary events are now candidate events
        {
            candidates.insert(candidates.end(), std::make_move_iterator(col.dijets.begin()),
                            std::make_move_iterator(col.dijets.end()));
            col.dijets.erase(col.dijets.begin(), col.dijets.end());
        }
    }

    candidates.shrink_to_fit();

    std::sort(candidates.begin(), candidates.end(), //Sort the candidate events so that the one with the biggest kt is first
              [](dijet_specs &s1, dijet_specs &s2) { return (s1.kt > s2.kt); });
    //std::cout << "candidates.size() =  " << candidates.size() << std::endl;

    for (auto &candidate : candidates)
    {
        if (fits_into_ps(filtered_scatterings, candidate))
        {
            filtered_scatterings.push_back(candidate);
        }
    }
    //std::cout << "filtered_scatterings.size() =  " << filtered_scatterings.size() << std::endl;
}

auto mc_glauber_style_report
(
    std::vector<Coll> &collisions, 
    const xsectval &sigma_inel, 
    const uint &N_events, 
    const uint &nBins, 
    const double *const binsLow, 
    const double *const binsHigh,
    std::ostream out_stream
) noexcept -> void
{
    // Make sure that no rounding downwards.
    double eps = 0.1/N_events;

    // int comp = 1;
    for (uint8_t comp = 1; comp <= 4; ++comp) {

        // Print out the header.
        out_stream << std::endl;
        out_stream << "# sigma_inel = " << sigma_inel << " mb" << std::endl;

        switch (comp)
        {
            case 1:
                std::sort(collisions.begin(), collisions.end(), compNpart);
                out_stream << "Using Npart for centrality determination" 
                           << std::endl << std::endl
                           << "# cent%\t nPmax\t nPmin\t <b>\t <nPart> <nColl> <T_AA>"
                           << std::endl;
                break;
            case 2:
                std::sort(collisions.begin(), collisions.end(), compNcoll);
                out_stream << "Using Ncoll for centrality determination" 
                           << std::endl << std::endl
                           << "# cent%\t nCmax\t nCmin\t <b>\t <nPart> <nColl> <T_AA>"
                           << std::endl;
                break;
            case 3:
                std::sort(collisions.begin(), collisions.end(), compNtwo);
                out_stream << "Using two-component (ancestors) model for centrality determination"
                           << std::endl << std::endl
                           << "# cent%\t nAmax\t nAmin\t <b>\t <nPart> <nColl> <T_AA>"
                           << std::endl;
                break;
            default: //case 4
                std::sort(collisions.begin(), collisions.end(), compET);
                out_stream << "Using sumET model for centrality determination" 
                           << std::endl << std::endl
                           << "# cent%\t ETmax\t ETmin\t <b>\t <nPart> <nColl> <T_AA>"
                           << std::endl;
        }
        std::reverse(collisions.begin(), collisions.end());

        std::vector<Coll> centrality;
        // Number of collisions within centrality class i.
        for (uint i = 0; i < nBins; i++)
        {
            uint lower = static_cast<uint>(binsLow[i]*N_events+eps);
            uint upper = static_cast<uint>(binsHigh[i]*N_events+eps);
            centrality.clear();
            while (lower<upper)
            {
                centrality.push_back(collisions.at(lower++));
            }

            // Print the centrality selection.
            out_stream << binsLow[i]*100 << " " << binsHigh[i]*100 << '\t';

            switch (comp)
            {
                case 1:
                    out_stream << centrality[0].getNpart() << '\t' 
                               << centrality[centrality.size() - 1].getNpart() << '\t';
                    break;
                case 2:
                    out_stream << centrality[0].getNcoll() << '\t' 
                               << centrality[centrality.size() - 1].getNcoll() << '\t';
                    break;
                case 3:
                    out_stream << centrality[0].getNtwo() << '\t' 
                               << centrality[centrality.size() - 1].getNtwo() << '\t';
                    break;
                default : //case 4
                    out_stream << centrality[0].getET() << '\t' 
                               << centrality[centrality.size() - 1].getET() << '\t';
                    break;
            }

            out_stream << calc_ave<double>(centrality, &Coll::getB) << '\t' 
                       << calc_ave<ulong>(centrality, &Coll::getNpart) << '\t' 
                       << calc_ave<ulong>(centrality, &Coll::getNcoll) << '\t' 
                       << calc_ave<ulong>(centrality, &Coll::getNcoll)/sigma_inel << '\t' 
                       // << calc_ave<double>(centrality, &Coll::getET) << '\t' 
                       << std::endl;      
        }
    }
}

#include <type_traits>

template<bool temp_bool, typename std::enable_if <!temp_bool> :: type* = nullptr>
auto collide_nuclei
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
                newpair.push_end_states_to_collider_frame();

                if (AA_params.reduce_nucleon_energies)
                {
                    newpair.reduce_energy();
                }
            }                    
            newpair.wound();
            binary_collisions.push_back(newpair);
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
                newpair.push_end_states_to_collider_frame();

                if (AA_params.reduce_nucleon_energies)
                {
                    newpair.reduce_energy();
                }
            }
            newpair.wound();
            binary_collisions.push_back(newpair);
        }
    }

    if (verbose)
    {
        std::cout << "Bruteforced " << n_pairs << " pairs, got " << binary_collisions.size()+nof_softs << " collisions, of which softs "<< nof_softs<< " and hards "<< binary_collisions.size()<<" , momentum threshold broke " << mombroke << " times" << std::endl;
    }
}

auto calculate_tA
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

auto calculate_tAB
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

auto calculate_sum_tpp
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

template<bool temp_bool, typename std::enable_if <temp_bool> :: type* = nullptr>
auto collide_nuclei
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

    spatial tAA_0 = 29.5494;//30.5//calculate_tAB({0,0,0}, pro, pro, AA_params.Tpp);
    spatial tBB_0 = 29.5494;//30.5//calculate_tAB({0,0,0}, tar, tar, AA_params.Tpp);
    
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
            newpair.push_end_states_to_collider_frame();

            if (AA_params.reduce_nucleon_energies)
            {
                newpair.reduce_energy();
            }
        }
        newpair.wound();
        binary_collisions.push_back(newpair);
    }

    if (verbose)
    {
        std::cout << "Bruteforced " << n_pairs << " pairs, got " << binary_collisions.size()+nof_softs << " collisions, of which softs "<< nof_softs<< " and hards "<< binary_collisions.size()<<" , momentum threshold broke " << mombroke << " times, skipped "<< skipped << " pairs that were too far apart" << std::endl;
    }
}

// return an evenly spaced 1d grid of doubles.
auto linspace
(
    const double &first, 
    const double &last, 
    const uint16_t &len
) noexcept -> std::vector<double>
{
    std::vector<double> result(len);
    double step = (last-first) / (len - 1);
    for (uint16_t i=0; i<len; i++) { result[i] = first + i*step; }
    return result;
}

// return an  1d grid of doubles distributed evenly in logarithm
auto loglinspace
(
    const double &first, 
    const double &last, 
    const uint16_t &len
) noexcept -> std::vector<double>
{
    std::vector<double> result(len);
    double step = pow(last/first, 1.0/(len-1.0));
    result[0] = first;
    for (uint16_t i=1; i<len; i++) { result[i] = result[i-1]*step; }
    return result;
}

auto calculate_spatial_sigma_jets_mf_MC
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
    std::vector<spatial> grid1 = linspace(lower_sumTpp_limit, upper_sumTpp_limit, dim_Ns[0]);
    std::vector<spatial> grid2 = linspace(lower_sumTpp_limit, upper_sumTpp_limit, dim_Ns[1]);
    std::vector<spatial> grid3 = linspace(kt02, mand_s, dim_Ns[2]);
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
            PrintThread{} <<'\r'<<--running_count<<" left of "
                          <<num_elements<<" grid points to be calculated";
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

auto calculate_spatial_sigma_jets_mf
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
    std::vector<spatial> grid1 = linspace(lower_sumTpp_limit, upper_sumTpp_limit, dim_Ns[0]);
    std::vector<spatial> grid2 = linspace(lower_sumTpp_limit, upper_sumTpp_limit, dim_Ns[1]);
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
            PrintThread{} <<'\r'<<--running_count<<" left of "
                          <<num_elements<<" grid points to be calculated";
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

auto calculate_spatial_sigma_jets_full
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
    std::vector<spatial> grid1 = linspace(lower_limits[0], upper_limits[0], dim_Ns[0]);
    std::vector<spatial> grid2 = linspace(lower_limits[1], upper_limits[1], dim_Ns[1]);
    std::vector<spatial> grid3 = linspace(lower_limits[2], upper_limits[2], dim_Ns[2]);
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
            PrintThread{} <<'\r'<<--running_count<<" left of "<<num_elements<<" grid points to be calculated";
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

auto read_sigma_jets_mf_MC
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

auto read_sigma_jets_mf
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

auto read_sigma_jets_full
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

auto calculate_sigma_1jet_analytical
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

auto calculate_sigma_jet_analytical
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

auto calculate_sigma_dijet_analytical
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

auto print_3d_histo
(
    const histo_3d &histo_, 
    const std::string &filename,
    const double normalization = 1.0,
    const bool divide_tot = true
) noexcept -> void
{
    std::ofstream file;

    file.open(filename);

    auto [ histo, uf, of, xs, tot ] = histo_.project_2d(1);
    auto [ xs0, xs1 ] = xs;

    file << "///Total count: "<<tot<<std::endl;
    file << "///underflow: "<<uf<<std::endl;
    file << "///overflow: "<<of<<std::endl;

    file << "///y bin walls:  ";
    for (auto y : xs1)
    {
        file<<y<<' ';
    }
    file << std::endl;

    auto norm = (divide_tot)? normalization : normalization * static_cast<double>(tot);

    ulong n_x_bins = xs0.size()-1;
    ulong n_y_bins = xs1.size()-1;

    for (ulong i = 0; i < n_x_bins; i++)
    {
        file << xs0[i] << ' ' << xs0[i+1] << ' ';

        for (ulong j = 0; j < n_y_bins; j++)
        {   
            file << histo[i][j] * norm << ' ';
        }
        file << std::endl;
    }

    auto [ histo2, uf2, of2, xs2, tot2 ] = histo_.project_2d(2);
    auto [ xs02, xs12 ] = xs2;
    file << "///y bin walls:  ";
    for (auto y : xs12)
    {
        file<<y<<' ';
    }
    file << std::endl;

    norm = (divide_tot)? normalization : normalization * static_cast<double>(tot2);

    n_x_bins = xs02.size()-1;
    n_y_bins = xs12.size()-1;

    for (ulong i = 0; i < n_x_bins; i++)
    {
        file << xs02[i] << ' ' << xs02[i+1] << ' ';

        for (ulong j = 0; j < n_y_bins; j++)
        {   
            file << histo2[i][j] * norm << ' ';
        }
        file << std::endl;
    }
    std::cout<<"printed to "<<filename<<std::endl;
    file.close();
}

auto print_2d_histo
(
    const histo_2d &histo_, 
    const std::string &filename,
    const double normalization = 1.0,
    const bool divide_tot = true
) noexcept -> void
{
    std::ofstream file;

    file.open(filename);

    auto [ histo, uf, of, xs, tot ] = histo_.get_histo();
    auto [ uf0, uf1 ] = uf;
    auto [ of0, of1 ] = of;
    auto [ xs0, xs1 ] = xs;

    file << "///Total count: "<<tot<<std::endl;
    file << "///underflow: "<<uf0<<' '<<uf1<<std::endl;
    file << "///overflow: "<<of0<<' '<<of1<<std::endl;
    file << "///normalization: "<<normalization<<std::endl;

    file << "///y bin walls:  ";
    for (auto y : xs1)
    {
        file<<y<<' ';
    }
    file << std::endl;

    auto norm = (divide_tot)? normalization : normalization * static_cast<double>(tot);

    auto n_x_bins = histo.size();
    auto n_y_bins = histo[0].size();

    for (uint8_t i = 0; i < n_x_bins; i++)
    {
        file << xs0[i] << ' ' << xs0[i+1] << ' ';

        for (uint8_t j = 0; j < n_y_bins; j++)
        {   
            file << histo[i][j] * norm << ' ';
        }
        file << std::endl;
    }
    std::cout<<"printed to "<<filename<<std::endl;
    file.close();
}

auto save_event
(
    std::ofstream &event_file, 
    const std::vector<nucleon> &pro, 
    const std::vector<nucleon> &tar, 
    const spatial &impact_parameter
) noexcept -> void
{
    event_file << "{"<<std::endl;
    event_file << "    "<<impact_parameter<<','<<std::endl;
    event_file << "    {"<<std::endl;
    for (const auto & n : pro)
    {
        event_file << "        {"<<std::endl;
        event_file << "            "<<n.co<<','<<std::endl;
        event_file << "            "<<n.mom<<','<<std::endl;
        event_file << "            "<<n.wounded<<','<<std::endl;
        event_file << "            "<<n.is_neutron<<std::endl;
        if (&n != &pro.back())
        {
            event_file << "        },"<<std::endl;
        }
        else
        {
            event_file << "        }"<<std::endl;
        }
    }
    event_file << "    },"<<std::endl;
    event_file << "    {"<<std::endl;
    for (auto n : tar)
    {
        event_file << "        {"<<std::endl;
        event_file << "            "<<n.co<<std::endl;
        event_file << "            "<<n.mom<<std::endl;
        event_file << "            "<<n.wounded<<std::endl;
        event_file << "            "<<n.is_neutron<<std::endl;
        if (&n != &tar.back())
        {
            event_file << "        },"<<std::endl;
        }
        else
        {
            event_file << "        }"<<std::endl;
        }
    }
    event_file << "    }"<<std::endl;
    event_file << "},"<<std::endl;
}

auto save_event
(
    std::ofstream &event_file, 
    const std::vector<nucleon> &pro, 
    const std::vector<nucleon> &tar, 
    const spatial &impact_parameter,
    const std::vector<dijet_specs> &filtered_scatterings
) noexcept -> void
{
    event_file << "{"<<std::endl;
    event_file << "    "<<impact_parameter<<','<<std::endl;
    event_file << "    {"<<std::endl;
    for (const auto & n : pro)
    {
        event_file << "        {"<<std::endl;
        event_file << "            "<<n.co<<','<<std::endl;
        event_file << "            "<<n.mom<<','<<std::endl;
        event_file << "            "<<n.wounded<<','<<std::endl;
        event_file << "            "<<n.is_neutron<<std::endl;
        if (&n != &pro.back())
        {
            event_file << "        },"<<std::endl;
        }
        else
        {
            event_file << "        }"<<std::endl;
        }
    }
    event_file << "    },"<<std::endl;
    event_file << "    {"<<std::endl;
    for (const auto & n : tar)
    {
        event_file << "        {"<<std::endl;
        event_file << "            "<<n.co<<','<<std::endl;
        event_file << "            "<<n.mom<<','<<std::endl;
        event_file << "            "<<n.wounded<<','<<std::endl;
        event_file << "            "<<n.is_neutron<<std::endl;
        if (&n != &tar.back())
        {
            event_file << "        },"<<std::endl;
        }
        else
        {
            event_file << "        }"<<std::endl;
        }
    }
    event_file << "    },"<<std::endl;
    event_file << "    {"<<std::endl;
    for (const auto & sc : filtered_scatterings)
    {
        event_file << "        {"<<std::endl;
        event_file << "            "<<sc.kt<<','<<std::endl;
        event_file << "            "<<sc.y1<<','<<std::endl;
        event_file << "            "<<sc.y2<<','<<std::endl;
        event_file << "            "<<sc.init1<<','<<std::endl;
        event_file << "            "<<sc.init2<<','<<std::endl;
        event_file << "            "<<sc.final1<<','<<std::endl;
        event_file << "            "<<sc.final2<<std::endl;
        if (&sc != &filtered_scatterings.back())
        {
            event_file << "        },"<<std::endl;
        }
        else
        {
            event_file << "        }"<<std::endl;
        }
    }
    event_file << "    }"<<std::endl;
    event_file << "},"<<std::endl;
}

auto print_1d_histo
(
    const histo_1d &histo_, 
    const std::string &filename,
    const double normalization = 1.0,
    const bool divide_tot = true
) noexcept -> void
{
    std::ofstream file;

    file.open(filename);

    auto [ histo, uf, of, xs, tot ] = histo_.get_histo();

    file << "///Total count: "<<tot<<std::endl;
    file << "///underflow: "<<uf<<std::endl;
    file << "///overflow: "<<of<<std::endl;

    file << "///x bin walls:  ";
    for (auto x : xs)
    {
        file<<x<<' ';
    }
    file << std::endl;

    auto norm = (divide_tot)? normalization : normalization * static_cast<double>(tot);

    for (auto e : histo)
    {
        file<< e * norm <<' ';
    }
    file << std::endl;

    std::cout<<"printed to "<<filename<<std::endl;
    file.close();
}

auto calculate_sigma_jets_for_MC
(
    std::shared_ptr<LHAPDF::GridPDF> p_pdf, 
    const momentum &mand_s, 
    const momentum &kt02,  
    pqcd::sigma_jet_params params,
    const uint16_t n_divisions = 100
)
{
    const std::vector<double> mand_ss{linspace(kt02, mand_s, n_divisions)};
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

auto calculate_max_dsigmas_for_MC
(
    const momentum &kt0,
    const momentum &sqrt_s,
    std::shared_ptr<LHAPDF::GridPDF> p_pdf,    
    pqcd::sigma_jet_params params,
    const uint16_t n_divisions = 100
)
{
    const std::vector<double> sqrt_ss{linspace(kt0, sqrt_s, n_divisions)};
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
            auto max_dsigma = find_max_dsigma(kt0, ss, p_pdf, params);
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

auto prepare_sigma_jets
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
            auto [max_dsigma, err] = find_max_dsigma(kt0, sqrt_s, p_pdf, jet_params);
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
            auto [mand_ss, sigmas] = calculate_sigma_jets_for_MC(p_pdf, mand_s, kt02, jet_params);
            auto sigma_jet = linear_interpolator(mand_ss, sigmas);
            std::cout<<"done!"<<std::endl;

            dijet_norm = sigmas.back();

            std::cout<<"Calculating envelope..."<<std::flush;
            auto [sqrt_ss, max_dsigmas] = calculate_max_dsigmas_for_MC(kt0, sqrt_s, p_pdf, jet_params);
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
            auto [max_dsigma, err] = find_max_dsigma(kt0, sqrt_s, p_pdf, jet_params);
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
            auto [mand_ss, sigmas] = calculate_sigma_jets_for_MC(p_pdf, mand_s, kt02, jet_params);
            auto sigma_jet = linear_interpolator(mand_ss, sigmas);
            std::cout<<"done!"<<std::endl;

            dijet_norm = sigmas.back();

            std::cout<<"Calculating envelope..."<<std::flush;
            auto [sqrt_ss, max_dsigmas] = calculate_max_dsigmas_for_MC(kt0, sqrt_s, p_pdf, jet_params);
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
            variant_sigma_jet sigma_jet = read_sigma_jets_mf("sigma_jet_grid_mf.dat");
            std::cout<<"done!"<<std::endl;

            dijet_norm = 93.7604; //sigma_jet with ave nPDF

            std::cout<<"Calculating envelope..."<<std::flush;
            auto [max_dsigma, err] = find_max_dsigma(kt0, sqrt_s, p_pdf, jet_params);
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
                = calculate_spatial_sigma_jets_mf
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
            auto [max_dsigma, err] = find_max_dsigma(kt0, sqrt_s, p_pdf, jet_params);
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
            variant_sigma_jet sigma_jet = read_sigma_jets_mf_MC("sigma_jet_grid_mf_MC.dat");
            std::cout<<"done!"<<std::endl;

            dijet_norm = 93.7604; //sigma_jet with ave nPDF

            std::cout<<"Calculating envelope..."<<std::flush;
            auto [sqrt_ss, max_dsigmas] = calculate_max_dsigmas_for_MC(kt0, sqrt_s, p_pdf, jet_params);
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
                = calculate_spatial_sigma_jets_mf_MC
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
            auto [sqrt_ss, max_dsigmas] = calculate_max_dsigmas_for_MC(kt0, sqrt_s, p_pdf, jet_params);
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

volatile std::atomic_bool user_aborted = false;
void abort_handler(int num)
{
    user_aborted = true;
}

void calculate_and_save_nuclei_TAs_TAAs
(
    const nucleus_generator::nucleus_params &nuc_params,
    std::shared_ptr<std::mt19937> eng,
    std::shared_ptr<ars> radial_sampler
)
{
    auto num_nuclei = 10000;
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

    generate_nuclei
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

void calculate_and_save_average_nuclei_TAs_TAAs
(
    const nucleus_generator::nucleus_params &nuc_params,
    std::shared_ptr<std::mt19937> eng,
    std::shared_ptr<ars> radial_sampler
)
{
    auto num_nuclei = 10000;
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

    generate_nuclei
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
            lock_rad.~unique_lock();

            std::array<std::array<double,201>,201> TA_grid_dummy;
            for (uint16_t i=0; i<201; i++)
            {
                TA_grid_dummy[i] = std::array<double,201>();
                for (uint16_t j=0; j<201; j++)
                {
                    TA_grid_dummy[i][j] = calculate_tA(grid[i][j], pro, Tp);
                    TA_grid_dummy[i][j] += calculate_tA(grid[i][j], tar, Tp);
                    TA_grid_dummy[i][j] /= num_nuclei;
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

#ifndef IS_AA
#define IS_AA false
#endif
#ifndef SPATIAL_NPDFS
#define SPATIAL_NPDFS false
#endif
#ifndef IS_MOM_CONS
#define IS_MOM_CONS false
#endif
#ifndef IS_MOM_CONS2
#define IS_MOM_CONS2 true
#endif

int main()
{ 
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
    std::string   name_postfix = "_pp_100k_B=0_UUSMC_fr",
                  event_file_name = "event_log"+name_postfix+".dat";
    uint32_t      desired_N_events      = 100000,
                  AA_events             = 0;
    const spatial b_min                 = 0,
                  b_max                 = 0;
    std::mutex AA_events_mutex; 
    std::cout<<"Doing the run "<<name_postfix<<std::endl;
    //auto eng = std::make_shared<std::mt19937>(static_cast<ulong>(1));
    auto eng = std::make_shared<std::mt19937>(static_cast<ulong>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<double> unirand{0.0, 1.0};

    //Parameters for the nuclei
    const spatial rad_min=0,
                  rad_max=20;
    const std::function<double(const double&)> rad_pdf{[](const double & x)
    {
        return x*x/(1+exp((x-6.624)/0.549));
    }};
    auto radial_sampler = std::make_shared<ars>(rad_pdf, rad_min, rad_max);
    std::mutex radial_sampler_mutex; 
    nucleus_generator::nucleus_params nuc_params = 
    {
        /* .N=                    */208, 
        /* .Z=                    */82, 
        /* .min_distance=         */0.4, 
        /* .shift_cms=            */true, 
        /* .correct_overlap_bias= */true
    };
    
    //Parameters for the hard collisions
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
    /*spatial_pdfs=             */SPATIAL_NPDFS,
    /*spatial_averaging=        */false,
    /*calculate_end_state=      */true,
    /*reduce_nucleon_energies=  */IS_MOM_CONS,
    /*sigma_inel_for_glauber=   */sigma_inel_for_glauber,
    /*Tpp=                      */Tpp,
    /*normalize_to=             */B2_normalization_mode::inelastic,
    /*sqrt_s=                   */sqrt_s,
    /*energy_threshold=         */kt0
    };
    //std::vector<Coll> collisions_for_reporting;
    const std::vector<double> kt_bins{loglinspace(kt0, sqrt_s/2.0, 21)};
    const std::vector<double> y_bins{linspace(-ylim, ylim, 40)};
    const std::vector<double> b_bins{linspace(b_min, b_max, 21)};

    //sigma_jet parameters
    /*const bool read_sigmajets_from_file = false;*/
    pqcd::diff_sigma::params diff_params = pqcd::diff_sigma::params(
    /*projectile_with_npdfs=    */IS_AA,
    /*target_with_npdfs=        */IS_AA,
    /*isoscalar_projectile=     */false,
    /*isoscalar_target=         */false,
    /*npdfs_spatial=            */coll_params.spatial_pdfs,
    /*npdf_setnumber=           */1,
    /*A=                        */208, //Pb 
    /*B=                        */208, //Pb
    /*ZA=                       */82,  //Pb
    /*ZB=                       */82   //Pb
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

    auto
    [
        dijet_norm,
        power_law,
        sigma_jet,
        mand_s_array,
        envelope_maximum,
        sqrt_s_array
    ] = 
    prepare_sigma_jets
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
 
    total_energy.open("total_energies"+name_postfix+".dat");
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
    
        std::find_if
        (
            std::execution::par, 
            event_indexes.begin(), 
            event_indexes.end(), 
            [&,&sigma_jet=sigma_jet,&power_law=power_law,&envelope_maximum=envelope_maximum](const uint64_t index) 
            {
                uint32_t NColl = 0;
                std::vector<nn_coll> binary_collisions;
                std::vector<dijet_specs> filtered_scatterings;

                //B^2 from a uniform distribution
                spatial impact_parameter = sqrt(b_min*b_min + unirand(*eng)*(b_max*b_max-b_min*b_min)); 

                std::unique_lock<std::mutex> lock_rad(radial_sampler_mutex);
                auto [pro, tar] = generate_nuclei
                (
                    nuc_params, 
                    sqrt_s, 
                    impact_parameter, 
                    eng, 
                    radial_sampler, 
                    read_nuclei_from_file, 
                    verbose
                );
                lock_rad.~unique_lock();

                do
                {
                    binary_collisions.erase(binary_collisions.begin(), binary_collisions.end());
                    if (verbose) std::cout<<"impact_parameter: "<<impact_parameter<<std::endl;

                    collide_nuclei<SPATIAL_NPDFS>
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
                        const std::lock_guard<std::mutex> lock(radial_sampler_mutex);
                        std::tie(pro, tar) = generate_nuclei
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
                    save_event(event_file, pro, tar, impact_parameter);
                }
                else if (end_state_filtering)
                {
                    momentum ET=0, E=0;
                    filter_end_state(binary_collisions, filtered_scatterings, IS_MOM_CONS2, sqrt_s);
                    binary_collisions.erase(binary_collisions.begin(), binary_collisions.end());
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
                        save_event(event_file, pro, tar, impact_parameter, filtered_scatterings);
                    }

//                    filtered_scatterings.erase(filtered_scatterings.begin(), filtered_scatterings.end());
                }
                else
                {
//                    binary_collisions.erase(binary_collisions.begin(), binary_collisions.end());
                }

                /*if (verbose || (nof_collisions%100)==0)
                {
                    std::cout << std::endl << "A+A collided thus far: " << nof_collisions << ", of which events thus far: " << AA_events << std::endl << std::endl;
                }*/

                {
                    const std::lock_guard<std::mutex> lock(AA_events_mutex);
                    AA_events++;
                    std::cout <<"\rA+A collisions calculated: " << AA_events << std::flush;
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
                bool ret_value = user_aborted;
                return ret_value;
            }
        );
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

    print_2d_histo
    (
        jets, 
        "sigma1jet_sim"+name_postfix+".dat", 
        2.0 * dijet_norm
    );

    print_2d_histo
    (
        jets, 
        "dNdpTdy_sim"+name_postfix+".dat", 
        1.0,
        false
    );

    print_2d_histo
    (
        dijets, 
        "sigmadijet_sim"+name_postfix+".dat", 
        dijet_norm
    );

    print_1d_histo
    (
        dETdy, 
        "dETdy_sim"+name_postfix+".dat", 
        1.0 / AA_events,
        false
    );
    
    print_1d_histo
    (
        dEdy, 
        "dEdy_sim"+name_postfix+".dat", 
        1.0 / AA_events,
        false
    );
    
    print_1d_histo
    (
        dNdy, 
        "dNdy_sim"+name_postfix+".dat", 
        1.0 / AA_events,
        false
    );
    
    print_1d_histo
    (
        dETdeta, 
        "dETdeta_sim"+name_postfix+".dat", 
        1.0 / AA_events,
        false
    );
    
    print_1d_histo
    (
        dEdeta, 
        "dEdeta_sim"+name_postfix+".dat", 
        1.0 / AA_events,
        false
    );

    print_1d_histo
    (
        dETdb, 
        "dETdb_sim"+name_postfix+".dat", 
        1.0,
        true
    );
    
    print_1d_histo
    (
        dEdb, 
        "dEdb_sim"+name_postfix+".dat",  
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
