//Copyright (c) 2021 Mikko Kuha

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <future>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <stdexcept>
#include <string>
#include <sstream>
#include <thread>
#include <tuple>
#include <variant>
#include <vector>

#include "ars.hpp"
#include "LHAPDF/GridPDF.h"
#include "linear_interpolator.hpp"
#include "linterp.h"
#include "nn_coll.hpp"
#include "nucleon.hpp"
#include "nucleus_generator.hpp"
#include "pqcd.hpp"
#include "typedefs.hpp"

double largest_taa=0, largest_sumtpp=0;
double smallest_taa=1000, smallest_sumtpp=1000; 

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

auto read_nucleon_configs_from_file() noexcept
{
    std::ifstream input("mc_initial.dat");

    auto pro = std::make_unique<std::vector<coords>>();
    auto tar = std::make_unique<std::vector<coords>>();

    if (input.is_open())
    {
        std::string line;

        for (size_t i = 0; i < 8; i++)
        {
            std::getline(input, line); //Skip the 8 unimportant rows
        }

        while (true) //Loop for first nucleus
        {
            std::getline(input, line);

            if (line.empty()) //The nuclei are separated by empty line, followed by 4 comment lines
            {
                for (size_t i = 0; i < 4; i++)
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

auto generate_nuclei(const nucleus_generator::nucleus_params &nuc_params, const momentum &sqrt_s, const spatial &impact_parameter, std::shared_ptr<std::mt19937> eng, std::shared_ptr<ars> radial_sampler, const bool &read_nuclei_from_file, const bool &verbose) noexcept
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
        pro = nucleus_generator::generate_nucleus(nuc_params, sqrt_s/2.0, -impact_parameter/2., eng, radial_sampler);
        tar = nucleus_generator::generate_nucleus(nuc_params, sqrt_s/2.0, impact_parameter/2., eng, radial_sampler);
        if (verbose) std::cout<<"Done!"<<std::endl;
    }

    if (verbose) std::cout<<"Shuffling nuclei..."<<std::flush;
    std::shuffle(pro.begin(), pro.end(), *eng);
    std::shuffle(tar.begin(), tar.end(), *eng);
    if (verbose) std::cout<<"Done!"<<std::endl;

    return std::make_tuple(std::move(pro), std::move(tar));
}

//Checks whether the saturation criterium allows the candidate to be added TODO
bool fits_into_ps(const std::vector<dijet_specs> &olds, const dijet_specs &candidate) noexcept
{
    if (olds.size() + static_cast<uint>(candidate.kt) < 0)
    {
        std::cout << "very strange" << std::endl;
    }
    return true;
}

void filter_end_state(std::vector<nn_coll> &binary_collisions, std::vector<dijet_specs> &filtered_scatterings)
{
    std::vector<dijet_specs> candidates;

    candidates.reserve(binary_collisions.size()*10); //10 events on average is just an overhead guess

    for (auto &col : binary_collisions) //all of the end states of the binary events are now candidate events
    {
        candidates.insert(candidates.end(), std::make_move_iterator(col.dijets.begin()),
                          std::make_move_iterator(col.dijets.end()));
        col.dijets.erase(col.dijets.begin(), col.dijets.end());
    }

    candidates.shrink_to_fit();

    std::sort(candidates.begin(), candidates.end(), //Sort the candidate events so that the one with the biggest kt is first
              [](dijet_specs &s1, dijet_specs &s2) { return (s1.kt > s2.kt); });
    std::cout << "candidates.size() =  " << candidates.size() << std::endl;

    for (auto &candidate : candidates)
    {
        if (fits_into_ps(filtered_scatterings, candidate))
        {
            filtered_scatterings.push_back(candidate);
        }
    }
    std::cout << "filtered_scatterings.size() =  " << filtered_scatterings.size() << std::endl;
}

void mc_glauber_style_report(std::vector<Coll> &collisions, const xsectval &sigma_inel, const uint &N_events, const uint &nBins, const double *const binsLow, const double *const binsHigh) noexcept
{
    // Make sure that no rounding downwards.
    double eps = 0.1/N_events;

    // int comp = 1;
    for (int comp = 1; comp <= 4; ++comp) {

        if (comp == 1) std::sort(collisions.begin(), collisions.end(), compNpart);
        else if (comp == 2) std::sort(collisions.begin(), collisions.end(), compNcoll);
        else if (comp == 3) std::sort(collisions.begin(), collisions.end(), compNtwo);
        else if (comp == 4) std::sort(collisions.begin(), collisions.end(), compET);
        std::reverse(collisions.begin(), collisions.end());

        // Print out the header.
        std::cout  << std::endl;
        std::cout << "# sigma_inel = " << sigma_inel << " mb" << std::endl;
        if (comp == 1) 
        {
            std::cout << "Using Npart for centrality determination" << std::endl << std::endl;;
            std::cout << "# cent%\t nPmax\t nPmin\t <b>\t <nPart> <nColl> <T_AA>" << std::endl;
        } 
        else if (comp == 2) 
        {
            std::cout << "Using Ncoll for centrality determination" << std::endl << std::endl;;
            std::cout << "# cent%\t nCmax\t nCmin\t <b>\t <nPart> <nColl> <T_AA>" << std::endl;
        } 
        else if (comp == 3) 
        {
            std::cout << "Using two-component (ancestors) model for centrality"
                << " determination" << std::endl << std::endl;;
            std::cout << "# cent%\t nAmax\t nAmin\t <b>\t <nPart> <nColl> <T_AA>" << std::endl;
        } 
        else if (comp == 4) 
        {
            std::cout << "Using sumET model for centrality determination" 
                << std::endl << std::endl;;
            std::cout << "# cent%\t ETmax\t ETmin\t <b>\t <nPart> <nColl> <T_AA>" << std::endl;
        }

        std::vector<Coll> centrality;
        // Number of collisions within centrality class i.
        for (uint i = 0; i < nBins; i++)
        {
            uint lower = static_cast<uint>(binsLow[i]*N_events+eps);
            uint upper = static_cast<uint>(binsHigh[i]*N_events+eps);
            centrality.clear();
            for (;lower<upper;lower++)
            {
                centrality.push_back(collisions.at(lower));
            }
      

            // Print the centrality selection.
            std::cout << binsLow[i]*100 << " " << binsHigh[i]*100 << '\t';
            if(comp == 1)
            {
                std::cout << centrality[0].getNpart() << '\t' 
                    << centrality[centrality.size() - 1].getNpart() << '\t';
            }
            else if(comp == 2)
            {
                std::cout << centrality[0].getNcoll() << '\t' 
                    << centrality[centrality.size() - 1].getNcoll() << '\t';
            } 
            else if(comp == 3)
            {
                std::cout << centrality[0].getNtwo() << '\t' 
                    << centrality[centrality.size() - 1].getNtwo() << '\t';
            } 
            else if(comp == 4)
            {
                std::cout << centrality[0].getET() << '\t' 
                    << centrality[centrality.size() - 1].getET() << '\t';
            } 

            std::cout << calc_ave<double>(centrality, &Coll::getB) << '\t' 
                    << calc_ave<ulong>(centrality, &Coll::getNpart) << '\t' 
                    << calc_ave<ulong>(centrality, &Coll::getNcoll) << '\t' 
                    << calc_ave<ulong>(centrality, &Coll::getNcoll)/sigma_inel << '\t' 
                    // << calc_ave<double>(centrality, &Coll::getET) << '\t' 
                    << std::endl;      
        }
    }
} 

void collide_nuclei(std::vector<nucleon> &pro, std::vector<nucleon> &tar, std::vector<nn_coll> &binary_collisions, 
                    std::variant<linear_interpolator, xsectval> sigma_jets, std::uniform_real_distribution<double> unirand, 
                    std::shared_ptr<std::mt19937> eng, const AA_collision_params &params, const bool &verbose) noexcept
{

    uint n_pairs = 0, mombroke = 0, nof_softs = 0;

    for (auto &A : pro)
    {
        for (auto &B : tar)
        {
            n_pairs++;
            if ((A.mom < params.energy_threshold) || (B.mom < params.energy_threshold))
            {
                mombroke++;
                continue;
            }

            nn_coll newpair(&A, &B, 2 * sqrt(A.mom * B.mom));
            if (params.mc_glauber_mode)
            {
                // "ball" diameter = distance at which two nucleons interact
                const spatial d2 = params.sigma_inel_for_glauber/(M_PI*10); // in fm^2
                const spatial dij2 = newpair.getcr_bsquared();
                
                if (dij2 > d2) //no collision
                {
                    continue;
                }
                //collision                    
                newpair.wound();
                binary_collisions.push_back(newpair);
            }
            else
            {
                xsectval sigma_jet;
                if (params.reduce_nucleon_energies)
                {
                    sigma_jet = std::get<linear_interpolator>(sigma_jets).value_at(pow(newpair.getcr_sqrt_s(), 2));
                }
                else //Single sigma_jet
                {
                    sigma_jet = std::get<xsectval>(sigma_jets);
                }
                
                newpair.calculate_xsects(sigma_jet, params.Tpp, newpair.getcr_bsquared(), params.normalize_to);
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
                if (params.calculate_end_state)
                {
                    pqcd::generate_bin_NN_coll(&newpair);
                    newpair.push_end_states_to_collider_frame();
                    newpair.reduce_energy();
                }
                newpair.wound();
                binary_collisions.push_back(newpair);
            }
        }
    }

    if (verbose)
    {
        std::cout << "Bruteforced " << n_pairs << " pairs, got " << binary_collisions.size()+nof_softs << " collisions, of which softs "<< nof_softs<< " and hards "<< binary_collisions.size()<<" , momentum threshold broke " << mombroke << " times" << std::endl;
    }

}

spatial calculate_tAB(const coords &b, const std::vector<nucleon> &pro, const std::vector<nucleon> &tar, const std::function<spatial(const spatial&)> Tpp) noexcept
{
    spatial tAB=0.0; //sum(T_pp(b_ij + b))
    size_t A=pro.size(), B=tar.size();

    for (size_t i=0; i<A; i++)
    {
        for (size_t j=0; j<B; j++)
        {
            tAB += Tpp((pro.at(i).co - tar.at(j).co + b).magt2());
        }   
    }
    
    return tAB;
}

spatial calculate_sum_tpp(const nucleon &nuc, const std::vector<nucleon> &nucleus, const std::function<spatial(const spatial&)> Tpp) noexcept
{
    spatial sum_tpp=0.0; //sum(T_pp(b_ii'))
    size_t A=nucleus.size();

    for (size_t i=0; i<A; i++)
    {
        sum_tpp += Tpp((nuc.co - nucleus.at(i).co).magt2());
    }
    return sum_tpp;
}

void collide_nuclei_with_spatial_pdfs(std::vector<nucleon> &pro, std::vector<nucleon> &tar, std::vector<nn_coll> &binary_collisions, const InterpMultilinear<4, xsectval> &sigma_jets,
                                      std::uniform_real_distribution<double> unirand, std::shared_ptr<std::mt19937> eng, 
                                      const AA_collision_params &params, const bool &verbose, const std::function<spatial(const spatial&)> Tpp) noexcept
{

    uint n_pairs = 0, mombroke = 0, nof_softs = 0;
    spatial sum_tppa=0, sum_tppb=0;

    spatial tAA_0 = calculate_tAB({0,0,0}, pro, pro, Tpp);
    spatial tBB_0 = calculate_tAB({0,0,0}, tar, tar, Tpp);
    
    if (verbose)
    {
        std::cout << "T_AA(0)= " << tAA_0 << ", T_BB(0)= " << tBB_0 << std::endl;
    }

    for (auto &A : pro)
    {
        sum_tppa = calculate_sum_tpp(A, pro, Tpp);
        
        for (auto &B : tar)
        {
            sum_tppb = calculate_sum_tpp(B, tar, Tpp);

            n_pairs++;
            if ((A.mom < params.energy_threshold) || (B.mom < params.energy_threshold))
            {
                mombroke++;
                continue;
            }

            nn_coll newpair(&A, &B, 2 * sqrt(A.mom * B.mom));
            xsectval sigma_jet;
            if (params.reduce_nucleon_energies) 
            {
                //TODO
            }
            else
            {
                sigma_jet = 0;//pqcd::calculate_spatial_sigma_jet_mf(p_p_pdf, p_n_pdf, &mand_s, &kt02, &jet_params, &sum_tppa, &sum_tppb, &tAA_0, &tBB_0);
            }
            
            newpair.calculate_xsects(sigma_jet, params.Tpp, newpair.getcr_bsquared(), params.normalize_to);

            if (verbose)
            {
                std::cout << "<T_pp>_i= " << sum_tppa << ", <T_pp>_j= " << sum_tppb << ", sigma_jet= " 
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
            if (params.calculate_end_state)
            {
                pqcd::generate_bin_NN_coll(&newpair);
                newpair.push_end_states_to_collider_frame();
                newpair.reduce_energy();
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

// return an evenly spaced 1-d grid of doubles.
std::vector<double> linspace(const double &first, const double &last, const uint16_t &len) noexcept
{
    std::vector<double> result(len);
    double step = (last-first) / (len - 1);
    for (uint16_t i=0; i<len; i++) { result[i] = first + i*step; }
    return result;
}

InterpMultilinear<4, xsectval> calculate_spatial_sigma_jets(const double &tolerance, std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
        std::shared_ptr<LHAPDF::GridPDF> p_n_pdf, const momentum &mand_s, const momentum &kt02, const pqcd::sigma_jet_params &jet_params, 
        const double &upper_tAA_0_limit, const double &lower_tAA_0_limit, const double &upper_sumTpp_limit, const double &lower_sumTpp_limit) noexcept
{
    const double marginal = 1.2; //20% more divisions than the tolerance gives us on the edges
    uint32_t max_threads = std::thread::hardware_concurrency(); //How many threads will be used in the calculation
    std::array<uint16_t,4> dim_Ns; //How many points to calculate in each dimension
    std::array<xsectval,16> corners;

    //std::mutex value_lock;
    std::function<void(std::promise<xsectval>, const spatial *const,const spatial *const,const spatial *const,const spatial *const)> sigma_jet_function
        = [=](std::promise<xsectval> prom, const spatial *const sum_tppa, const spatial *const sum_tppb, const spatial *const tAA_0, const spatial *const tBB_0)
    {
        //std::cout<<"HEP 1"<<std::endl;
        //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        //std::unique_lock<std::mutex> lock(value_lock);
        //xsectval value = pqcd::calculate_spatial_sigma_jet_mf(p_p_pdf, p_n_pdf, &mand_s, &kt02, &jet_params, sum_tppa, sum_tppb, tAA_0, tBB_0); 
        //std::cout<<"HEP 2"<<std::endl;
        prom.set_value(pqcd::calculate_spatial_sigma_jet_mf(p_p_pdf, p_n_pdf, &mand_s, &kt02, &jet_params, sum_tppa, sum_tppb, tAA_0, tBB_0));
        //lock.unlock();
        //std::cout<<"HEP 2"<<std::endl;
    };

    std::vector<std::promise<xsectval> > corner_promises(16);
    std::vector<std::future<xsectval> > corner_futures(16);
    std::vector<std::thread> corner_threads(16);
    for (uint8_t i=0; i<16; i++)
    {
        corner_futures[i] = corner_promises[i].get_future();
    }

    corner_threads[0]  = std::thread(sigma_jet_function, std::move(corner_promises[0]),  &upper_sumTpp_limit, &upper_sumTpp_limit, &upper_tAA_0_limit, &upper_tAA_0_limit);
    corner_threads[1]  = std::thread(sigma_jet_function, std::move(corner_promises[1]),  &upper_sumTpp_limit, &upper_sumTpp_limit, &upper_tAA_0_limit, &lower_tAA_0_limit);
    corner_threads[2]  = std::thread(sigma_jet_function, std::move(corner_promises[2]),  &upper_sumTpp_limit, &upper_sumTpp_limit, &lower_tAA_0_limit, &upper_tAA_0_limit);
    corner_threads[3]  = std::thread(sigma_jet_function, std::move(corner_promises[3]),  &upper_sumTpp_limit, &upper_sumTpp_limit, &lower_tAA_0_limit, &lower_tAA_0_limit);
    corner_threads[4]  = std::thread(sigma_jet_function, std::move(corner_promises[4]),  &upper_sumTpp_limit, &lower_sumTpp_limit, &upper_tAA_0_limit, &upper_tAA_0_limit);
    corner_threads[5]  = std::thread(sigma_jet_function, std::move(corner_promises[5]),  &upper_sumTpp_limit, &lower_sumTpp_limit, &upper_tAA_0_limit, &lower_tAA_0_limit);
    corner_threads[6]  = std::thread(sigma_jet_function, std::move(corner_promises[6]),  &upper_sumTpp_limit, &lower_sumTpp_limit, &lower_tAA_0_limit, &upper_tAA_0_limit);
    corner_threads[7]  = std::thread(sigma_jet_function, std::move(corner_promises[7]),  &upper_sumTpp_limit, &lower_sumTpp_limit, &lower_tAA_0_limit, &lower_tAA_0_limit);
    corner_threads[8]  = std::thread(sigma_jet_function, std::move(corner_promises[8]),  &lower_sumTpp_limit, &upper_sumTpp_limit, &upper_tAA_0_limit, &upper_tAA_0_limit);
    corner_threads[9]  = std::thread(sigma_jet_function, std::move(corner_promises[9]),  &lower_sumTpp_limit, &upper_sumTpp_limit, &upper_tAA_0_limit, &lower_tAA_0_limit);
    corner_threads[10] = std::thread(sigma_jet_function, std::move(corner_promises[10]), &lower_sumTpp_limit, &upper_sumTpp_limit, &lower_tAA_0_limit, &upper_tAA_0_limit);
    corner_threads[11] = std::thread(sigma_jet_function, std::move(corner_promises[11]), &lower_sumTpp_limit, &upper_sumTpp_limit, &lower_tAA_0_limit, &lower_tAA_0_limit);
    corner_threads[12] = std::thread(sigma_jet_function, std::move(corner_promises[12]), &lower_sumTpp_limit, &lower_sumTpp_limit, &upper_tAA_0_limit, &upper_tAA_0_limit);
    corner_threads[13] = std::thread(sigma_jet_function, std::move(corner_promises[13]), &lower_sumTpp_limit, &lower_sumTpp_limit, &upper_tAA_0_limit, &lower_tAA_0_limit);
    corner_threads[14] = std::thread(sigma_jet_function, std::move(corner_promises[14]), &lower_sumTpp_limit, &lower_sumTpp_limit, &lower_tAA_0_limit, &upper_tAA_0_limit);
    corner_threads[15] = std::thread(sigma_jet_function, std::move(corner_promises[15]), &lower_sumTpp_limit, &lower_sumTpp_limit, &lower_tAA_0_limit, &lower_tAA_0_limit);

    for (uint8_t i=0; i<16; i++)
    {
        corners[i] = corner_futures[i].get();
    }
    for (uint8_t i=0; i<16; i++)
    {
        corner_threads[i].join();
    }

    xsectval max_corner = *std::max_element(corners.begin(), corners.end());

    //sum_Tpp_A
    std::array<xsectval,8> differences{abs(corners[0]-corners[8]), abs(corners[1]-corners[9]), abs(corners[2]-corners[10]), abs(corners[3]-corners[11]), 
                                       abs(corners[4]-corners[12]), abs(corners[5]-corners[13]), abs(corners[6]-corners[14]), abs(corners[7]-corners[15])};
    xsectval max_diff = *std::max_element(differences.begin(), differences.end());
    dim_Ns[0] = static_cast<uint16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

    //sum_Tpp_B
    differences = std::array<xsectval,8>({abs(corners[0]-corners[4]), abs(corners[1]-corners[5]), abs(corners[2]-corners[6]), abs(corners[3]-corners[7]), 
                                       abs(corners[8]-corners[12]), abs(corners[9]-corners[13]), abs(corners[10]-corners[14]), abs(corners[11]-corners[15])});
    max_diff = *std::max_element(differences.begin(), differences.end());
    dim_Ns[1] = static_cast<uint16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

    //TAA(0)
    differences = std::array<xsectval,8>({abs(corners[0]-corners[2]), abs(corners[1]-corners[3]), abs(corners[4]-corners[6]), abs(corners[5]-corners[7]), 
                                       abs(corners[8]-corners[10]), abs(corners[9]-corners[11]), abs(corners[12]-corners[14]), abs(corners[13]-corners[15])});
    max_diff = *std::max_element(differences.begin(), differences.end());
    dim_Ns[2] = static_cast<uint16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));
    
    //TBB(0)
    differences = std::array<xsectval,8>({abs(corners[0]-corners[1]), abs(corners[2]-corners[3]), abs(corners[4]-corners[5]), abs(corners[6]-corners[7]), 
                                       abs(corners[8]-corners[9]), abs(corners[10]-corners[11]), abs(corners[12]-corners[13]), abs(corners[14]-corners[15])});
    max_diff = *std::max_element(differences.begin(), differences.end());
    dim_Ns[3] = static_cast<uint16_t>(ceil( ((max_diff/max_corner) / tolerance) * marginal ));

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
    std::vector<spatial> grid3 = linspace(lower_tAA_0_limit, upper_tAA_0_limit, dim_Ns[2]);
    std::vector<spatial> grid4 = linspace(lower_tAA_0_limit, upper_tAA_0_limit, dim_Ns[3]);
    std::vector< std::vector<spatial>::iterator > grid_iter_list;
    grid_iter_list.push_back(grid1.begin());
    grid_iter_list.push_back(grid2.begin());
    grid_iter_list.push_back(grid3.begin());
    grid_iter_list.push_back(grid4.begin());
  
    // total number of elements
    size_t num_elements = dim_Ns[0] * dim_Ns[1] * dim_Ns[2] * dim_Ns[3]; 

    std::ofstream sigma_jet_grid_file;
    sigma_jet_grid_file.open("sigma_jet_grid.dat", std::ios::out);
    sigma_jet_grid_file << "%mand_s=" << mand_s << " kt02=" << kt02 << " p_pdf=" << p_p_pdf->info().get_entry("SetIndex") << std::endl;
    sigma_jet_grid_file << "%num_elements=" << num_elements << std::endl;
    sigma_jet_grid_file << "%num_sum_T_pp_A=" << dim_Ns[0] << std::endl;
    sigma_jet_grid_file << "%num_sum_T_pp_B=" << dim_Ns[1] << std::endl;
    sigma_jet_grid_file << "%num_T_AA(0)=" << dim_Ns[2] << std::endl;
    sigma_jet_grid_file << "%num_T_BB(0)=" << dim_Ns[3] << std::endl;
    sigma_jet_grid_file << std::endl;
    sigma_jet_grid_file << "%sum_T_pp_A" << std::endl;
    for (auto g : grid1) sigma_jet_grid_file << g << ' ';
    sigma_jet_grid_file << std::endl << std::endl;
    sigma_jet_grid_file << "%sum_T_pp_B" << std::endl;
    for (auto g : grid2) sigma_jet_grid_file << g << ' ';
    sigma_jet_grid_file << std::endl << std::endl;
    sigma_jet_grid_file << "%T_AA(0)" << std::endl;
    for (auto g : grid3) sigma_jet_grid_file << g << ' ';
    sigma_jet_grid_file << std::endl << std::endl;
    sigma_jet_grid_file << "%T_BB(0)" << std::endl;
    for (auto g : grid4) sigma_jet_grid_file << g << ' ';
    sigma_jet_grid_file << std::endl << std::endl;
    sigma_jet_grid_file << "%sigma_jet_values" << std::endl;

    // fill in the values of f(x) at the gridpoints. 
    // we will pass in a contiguous sequence, values are assumed to be laid out C-style
    if (max_threads<dim_Ns[2]*dim_Ns[3])
    {
        std::cout<<"Too few possible threads!!!!"<<std::endl;
        std::cout<<"Need "<<dim_Ns[2]*dim_Ns[3]<<", have "<<max_threads<<std::endl;
        max_threads=dim_Ns[2]*dim_Ns[3];
    }

    std::vector<std::promise<xsectval> > grid_promises(dim_Ns[2]*dim_Ns[3]);
    std::vector<std::future<xsectval> > grid_futures(max_threads);
    std::vector<std::thread> grid_threads(max_threads);

    std::vector<xsectval> f_values(num_elements);
    size_t running_count=0;
    xsectval dummy;
    std::mutex dummy_lock;

    for (uint16_t i=14; i<dim_Ns[0]; i++)
    {
        for (uint16_t j=14; j<dim_Ns[1]; j++)
        {
            for (uint16_t ii=0; ii<max_threads; ii++)
            {
                grid_promises[ii] = std::promise<xsectval>();
                grid_futures[ii] = grid_promises[ii].get_future();
            }
            uint16_t gi=0;

    	    for (uint16_t k=0; k<dim_Ns[2]; k++)
            {
                for (uint16_t l=0; l<dim_Ns[3]; l++,gi++)
                {
                    grid_threads[gi] = std::thread(sigma_jet_function, std::move(grid_promises[gi]), &grid1[i], &grid2[j], &grid3[k], &grid4[l]);
            	}
            }

            gi=0;

    	    for (uint16_t k=0; k<dim_Ns[2]; k++)
            {
                for (uint16_t l=0; l<dim_Ns[3]; l++,gi++)
                {

                    std::unique_lock<std::mutex> lock(dummy_lock);
                    dummy = grid_futures[gi].get();
            	    f_values[i*dim_Ns[1]*dim_Ns[2]*dim_Ns[3] + j*dim_Ns[2]*dim_Ns[3] + k*dim_Ns[3] + l] = dummy;
                    sigma_jet_grid_file << dummy << ' ' <<std::flush;
                    dummy_lock.unlock();
                    std::cout<<"\rCalculated "<<++running_count<<" of "<<num_elements<<" grid points";
            	}
                sigma_jet_grid_file << std::endl;
            }
            sigma_jet_grid_file << std::endl;
            for (uint16_t ii=0; ii<max_threads; ii++)
            {
                grid_threads[ii].join();
            }
        }
        sigma_jet_grid_file << std::endl;
    }
    sigma_jet_grid_file.close();
    std::cout<<std::endl;

    return InterpMultilinear<4, xsectval>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data() + num_elements);
}

InterpMultilinear<4, xsectval> read_sigma_jets(const std::string &filename) noexcept
{

    std::ifstream input(filename);

    std::array<size_t,4> dim_Ns;
    std::vector<spatial> grid1, grid2, grid3, grid4;
    std::vector< std::vector<spatial>::iterator > grid_iter_list;
    std::vector<xsectval> f_values;

    if (input.is_open())
    {
        std::string line;
        std::getline(input, line); //#1 Don't need anything from here

        std::getline(input, line); //#2
        std::istringstream line_stream(line);
        size_t num_elements;
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

        std::getline(input, line); //#7 empty
        std::getline(input, line); //#8 Don't need anything from here
        std::getline(input, line); //#9
        line_stream = std::istringstream(line);
        spatial num;
        while (line_stream >> num)
        {
            grid1.push_back(num);
        }
        std::getline(input, line); //#10 empty
        std::getline(input, line); //#11 Don't need anything from here
        std::getline(input, line); //#12
        line_stream = std::istringstream(line);
        while (line_stream >> num)
        {
            grid2.push_back(num);
        }
        std::getline(input, line); //#13 empty
        std::getline(input, line); //#14 Don't need anything from here
        std::getline(input, line); //#15
        line_stream = std::istringstream(line);
        while (line_stream >> num)
        {
            grid3.push_back(num);
        }
        std::getline(input, line); //#16 empty
        std::getline(input, line); //#17 Don't need anything from here
        std::getline(input, line); //#18
        line_stream = std::istringstream(line);
        while (line_stream >> num)
        {
            grid4.push_back(num);
        }
        grid_iter_list.push_back(grid1.begin());
        grid_iter_list.push_back(grid2.begin());
        grid_iter_list.push_back(grid3.begin());
        grid_iter_list.push_back(grid4.begin());
        
        std::getline(input, line); //#19 empty
        std::getline(input, line); //#20 Don't need anything from here
        f_values.reserve(num_elements);
        uint16_t l=0;
        xsectval sigma_jet;

        for (uint16_t i=0; i<dim_Ns[0]; i++)
        {
            for (uint16_t j=0; j<dim_Ns[1]; j++)
            {
                for (uint16_t k=0; k<dim_Ns[2]; k++)
                {
                    std::getline(input, line);
                    line_stream = std::istringstream(line);
                    while (line_stream >> sigma_jet)
                    {
                        f_values[i*dim_Ns[1]*dim_Ns[2]*dim_Ns[3] + j*dim_Ns[2]*dim_Ns[3] + k*dim_Ns[3] + l] = sigma_jet;
                        l++;
                    }
                    l=0;
                }
                std::getline(input, line); //empty
            }
            std::getline(input, line); //empty
        }

        return InterpMultilinear<4, xsectval>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data() + num_elements);
    }

    std::cout<<"ERROR READING SIGMA_JETS"<<std::endl;
    
    return InterpMultilinear<4, xsectval>(grid_iter_list.begin(), dim_Ns.begin(), f_values.data(), f_values.data());
}

int main()
{
    //A lot of printing
    bool verbose = false;

    if (verbose) std::cout<<"Initializing..."<<std::flush;
    
    //General parameters for the simulation
    const bool read_nuclei_from_file = false, end_state_filtering = false;
    uint desired_N_events = 10000, AA_events = 0, nof_collisions = 0;
    const spatial b_min=0, b_max=20;
    auto eng = std::make_shared<std::mt19937>(static_cast<ulong>(1000));
    //auto eng = std::make_shared<std::mt19937>(static_cast<ulong>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<double> unirand(0.0, 1.0);

    //Parameters for the nuclei
    const spatial rad_min=0, rad_max=20;
    const std::function<double(const double&)> rad_pdf{[](const double & x){return x*x/(1+exp((x-6.62)/0.546));}};
    auto radial_sampler = std::make_shared<ars>(rad_pdf, rad_min, rad_max);
    nucleus_generator::nucleus_params nuc_params = {/*.N=*/208, /*.Z=*/82, /*.min_distance=*/0.4, /*.shift_cms=*/true, /*.correct_overlap_bias=*/true};
    
    //Parameters for the hard collisions
    const spatial proton_width_2 = pow(0.483, 2);
    const std::function<spatial(const spatial&)> Tpp{[&proton_width_2](const spatial &bsquared) { return exp(-bsquared / (4 * proton_width_2)) / (40 * M_PI * proton_width_2); }}; // 1/fm² = mb/fm² * 1/mb = 0.1 * 1/mb
    const xsectval sigma_inel_for_glauber = 70;//mb
    const momentum sqrt_s = 2760;//GeV
    const momentum mand_s = pow(sqrt_s, 2);//GeV^2
    momentum kt0 = 2.0;//GeV
    momentum kt02 = pow(kt0, 2);//GeV^2
    //rapidity ycut = 10.0;
    std::shared_ptr<LHAPDF::GridPDF> p_pdf(new LHAPDF::GridPDF("CT14lo", 0));
    const AA_collision_params coll_params{
    /*mc_glauber_mode=          */false,
    /*spatial_pdfs=             */true,
    /*calculate_end_state=      */false,
    /*reduce_nucleon_energies=  */false,
    /*sigma_inel_for_glauber=   */sigma_inel_for_glauber,
    /*Tpp=                      */Tpp,
    /*normalize_to=             */B2_normalization_mode::inelastic,
    /*sqrt_s=                   */sqrt_s,
    /*energy_threshold=         */2.0/*GeV*/};
    
    std::vector<nn_coll> binary_collisions;
    std::vector<dijet_specs> filtered_scatterings;
    std::vector<Coll> collisions_for_reporting;

    //sigma_jet stuff
    const bool read_sigmajets_from_file = true;
    const pqcd::sigma_jet_params jet_params = pqcd::sigma_jet_params(
    /*scale_choice=             */pqcd::scaled_from_kt,
    /*scalar=                   */1.0,
    /*projectile_with_npdfs=    */true,
    /*target_with_npdfs=        */true,
    /*isoscalar_projectile=     */false,
    /*isoscalar_target=         */false,
    /*npdf_setnumber=           */1,
    /*use_ses=                  */false);
    ////Pre-calculated sigma_jets
    //std::vector<double> xs({8, 124.235, 182.353, 211.412, 240.47, 269.529, 298.588, 327.647, 356.706, 385.764, 414.823, 443.882, 472.941, 502, 531.059, 560.117, 589.176, 618.235, 647.294, 676.353, 705.411, 734.47, 763.529, 792.588, 821.647, 850.705, 879.764, 908.823, 937.882, 966.941, 995.999, 1025.06, 1054.12, 1083.18, 1112.23, 1141.29, 1170.35, 1228.47, 1286.59, 1344.71, 1402.82, 1460.94, 1519.06, 1577.18, 1635.29, 1693.41, 1751.53, 1809.65, 1867.76, 1925.88, 1984, 2042.12, 2100.23, 2158.35, 2216.47, 2274.59, 2332.7, 2390.82, 2448.94, 2507.06, 2565.18, 2623.29, 2681.41, 2739.53, 2797.65, 2855.76, 2913.88, 2972, 3030.12, 3088.23, 3146.35, 3204.47, 3262.59, 3320.7, 3378.82, 3436.94, 3495.06, 3553.17, 3611.29, 3669.41, 3727.53, 3785.64, 3843.76, 3901.88, 3960, 4018.12, 4076.23, 4134.35, 4192.47, 4250.59, 4308.7, 4366.82, 4424.94, 4483.06, 4541.17, 4599.29, 4657.41, 4715.53, 4773.64, 4831.76, 4889.88, 4948, 5006.11, 5064.23, 5122.35, 5180.47, 5238.59, 5296.7, 5354.82, 5412.94, 5471.06, 5529.17, 5587.29, 5645.41, 5703.53, 5761.64, 5819.76, 5877.88, 5936, 5994.11, 6052.23, 6110.35, 6168.47, 6226.58, 6284.7, 6342.82, 6400.94, 6459.06, 6517.17, 6575.29, 6633.41, 6691.53, 6749.64, 6807.76, 6865.88, 6924, 6982.11, 7040.23, 7098.35, 7156.47, 7214.58, 7272.7, 7330.82, 7388.94, 7447.05, 7505.17, 7563.29, 7621.41, 7679.53, 7737.64, 7795.76, 7853.88, 7912, 7970.11, 8028.23, 8086.35, 8144.47, 8202.58, 8260.7, 8318.82, 8376.94, 8435.05, 8493.17, 8551.29, 8609.41, 8725.64, 8841.88, 8958.11, 9074.35, 9190.58, 9306.82, 9423.05, 9539.29, 9655.52, 9771.76, 9887.99, 10004.2, 10120.5, 10236.7, 10352.9, 10469.2, 10585.4, 10701.6, 10817.9, 10934.1, 11050.3, 11166.6, 11282.8, 11399.1, 11515.3, 11631.5, 11747.8, 11864, 11980.2, 12096.5, 12212.7, 12328.9, 12445.2, 12561.4, 12677.6, 12793.9, 12910.1, 13026.3, 13142.6, 13258.8, 13375.1, 13491.3, 13607.5, 13723.8, 13840, 13956.2, 14072.5, 14188.7, 14304.9, 14421.2, 14537.4, 14653.6, 14769.9, 14886.1, 15002.3, 15118.6, 15234.8, 15351.1, 15467.3, 15583.5, 15699.8, 15816, 16048.5, 16280.9, 16513.4, 16745.9, 16978.3, 17210.8, 17443.3, 17675.8, 17908.2, 18140.7, 18373.2, 18605.6, 18838.1, 19070.6, 19303, 19535.5, 19768, 20000.5, 20232.9, 20465.4, 20697.9, 20930.3, 21162.8, 21395.3, 21627.8, 21860.2, 22092.7, 22325.2, 22557.6, 22790.1, 23022.6, 23255, 23487.5, 23720, 23952.5, 24184.9, 24417.4, 24649.9, 24882.3, 25114.8, 25347.3, 25579.8, 25812.2, 26044.7, 26277.2, 26509.6, 26742.1, 26974.6, 27207, 27439.5, 27672, 27904.5, 28136.9, 28369.4, 28601.9, 28834.3, 29066.8, 29299.3, 29764.2, 30229.2, 30694.1, 31159, 31624, 32088.9, 32553.9, 33018.8, 33483.7, 33948.7, 34413.6, 34878.6, 35343.5, 35808.5, 36273.4, 36738.3, 37203.3, 37668.2, 38133.2, 38598.1, 39063, 39528, 39992.9, 40457.9, 40922.8, 41387.7, 41852.7, 42317.6, 42782.6, 43247.5, 43712.4, 44177.4, 44642.3, 45107.3, 45572.2, 46037.2, 46502.1, 46967, 47432, 47896.9, 48361.9, 48826.8, 49291.7, 49756.7, 50221.6, 50686.6, 51151.5, 51616.4, 52081.4, 52546.3, 53011.3, 53476.2, 53941.1, 54406.1, 54871, 55800.9, 56730.8, 57660.7, 58590.6, 59520.4, 60450.3, 61380.2, 62310.1, 63240, 64169.8, 65099.7, 66029.6, 66959.5, 67889.4, 68819.3, 69749.1, 70679, 71608.9, 72538.8, 73468.7, 74398.5, 75328.4, 76258.3, 77188.2, 78118.1, 79048, 79977.8, 80907.7, 81837.6, 82767.5, 83697.4, 84627.2, 85557.1, 86487, 87416.9, 88346.8, 89276.7, 90206.5, 91136.4, 92066.3, 92996.2, 93926.1, 94855.9, 95785.8, 96715.7, 97645.6, 98575.5, 99505.4, 100435, 101365, 102295, 104155, 106015, 107874, 109734, 111594, 113454, 115313, 117173, 119033, 120893, 122752, 124612, 126472, 128332, 130191, 132051, 133911, 135771, 137631, 139490, 141350, 143210, 145070, 146929, 148789, 150649, 152509, 154368, 156228, 158088, 159948, 161807, 163667, 165527, 167387, 169246, 171106, 172966, 174826, 176686, 178545, 180405, 182265, 184125, 185984, 187844, 189704, 191564, 193423, 197143, 200862, 204582, 208302, 212021, 215741, 219460, 223180, 226899, 230619, 234338, 238058, 241777, 245497, 249216, 252936, 256655, 260375, 264094, 267814, 271533, 275253, 278973, 282692, 286412, 290131, 293851, 297570, 301290, 305009, 308729, 312448, 316168, 319887, 323607, 327326, 331046, 334765, 338485, 342205, 345924, 349644, 353363, 357083, 360802, 364522, 371961, 379400, 386839, 394278, 401717, 409156, 416595, 424034, 431473, 438912, 446351, 453790, 461229, 468668, 476108, 483547, 490986, 498425, 505864, 513303, 520742, 528181, 535620, 543059, 550498, 557937, 565376, 572815, 580254, 587693, 595132, 602571, 610010, 617450, 624889, 632328, 639767, 647206, 654645, 662084, 669523, 676962, 684401, 699279, 714157, 729035, 743913, 758792, 773670, 788548, 803426, 818304, 833182, 848060, 862938, 877816, 892695, 907573, 922451, 937329, 952207, 967085, 981963, 996841, 1.01172e+06, 1.0266e+06, 1.04148e+06, 1.05635e+06, 1.07123e+06, 1.08611e+06, 1.10099e+06, 1.11587e+06, 1.13074e+06, 1.14562e+06, 1.1605e+06, 1.17538e+06, 1.19026e+06, 1.20513e+06, 1.22001e+06, 1.23489e+06, 1.24977e+06, 1.26465e+06, 1.27953e+06, 1.2944e+06, 1.30928e+06, 1.33904e+06, 1.36879e+06, 1.39855e+06, 1.42831e+06, 1.45806e+06, 1.48782e+06, 1.51758e+06, 1.54733e+06, 1.57709e+06, 1.60684e+06, 1.6366e+06, 1.66636e+06, 1.69611e+06, 1.72587e+06, 1.75562e+06, 1.78538e+06, 1.81514e+06, 1.84489e+06, 1.87465e+06, 1.90441e+06, 1.93416e+06, 1.96392e+06, 1.99367e+06, 2.02343e+06, 2.05319e+06, 2.08294e+06, 2.1127e+06, 2.14246e+06, 2.17221e+06, 2.20197e+06, 2.23172e+06, 2.26148e+06, 2.29124e+06, 2.32099e+06, 2.35075e+06, 2.38051e+06, 2.41026e+06, 2.44002e+06, 2.46977e+06, 2.49953e+06, 2.55904e+06, 2.61856e+06, 2.67807e+06, 2.73758e+06, 2.79709e+06, 2.8566e+06, 2.91612e+06, 2.97563e+06, 3.03514e+06, 3.09465e+06, 3.15417e+06, 3.21368e+06, 3.27319e+06, 3.3327e+06, 3.39222e+06, 3.45173e+06, 3.51124e+06, 3.57075e+06, 3.63027e+06, 3.68978e+06, 3.74929e+06, 3.8088e+06, 3.86832e+06, 3.92783e+06, 3.98734e+06, 4.04685e+06, 4.10637e+06, 4.16588e+06, 4.22539e+06, 4.2849e+06, 4.34442e+06, 4.40393e+06, 4.46344e+06, 4.52295e+06, 4.58247e+06, 4.64198e+06, 4.70149e+06, 4.761e+06, 4.88003e+06, 4.99905e+06, 5.11808e+06, 5.2371e+06, 5.35613e+06, 5.47515e+06, 5.59418e+06, 5.7132e+06, 5.83223e+06, 5.95125e+06, 6.07028e+06, 6.1893e+06, 6.30833e+06, 6.42735e+06, 6.54638e+06, 6.6654e+06, 6.78443e+06, 6.90345e+06, 7.02248e+06, 7.1415e+06, 7.26053e+06, 7.37955e+06, 7.49858e+06, 7.6176e+06, 7.7534e+06, 7.8892e+06, 8.025e+06, 8.1608e+06, 8.2966e+06, 8.4324e+06, 8.5682e+06, 8.704e+06, 8.8398e+06, 8.9756e+06, 9.1114e+06, 9.2472e+06, 9.383e+06, 9.5188e+06, 9.6546e+06, 9.7904e+06, 9.9262e+06, 1.0062e+07, 1.01978e+07, 1.03336e+07, 1.06052e+07, 1.08768e+07, 1.11484e+07, 1.142e+07, 1.16916e+07, 1.19632e+07, 1.22348e+07, 1.25064e+07, 1.2778e+07, 1.30496e+07, 1.33212e+07, 1.35928e+07, 1.38644e+07, 1.4136e+07, 1.44076e+07, 1.46792e+07, 1.49508e+07, 1.52224e+07, 1.5494e+07, 1.57656e+07, 1.60372e+07, 1.63088e+07, 1.65804e+07, 1.6852e+07, 1.71236e+07, 1.73952e+07, 1.76668e+07, 1.79384e+07, 1.821e+07, 1.84816e+07, 1.87532e+07, 1.90248e+07, 1.92964e+07, 1.9568e+07, 1.98396e+07, 2.01112e+07, 2.06544e+07, 2.11976e+07, 2.17408e+07, 2.2284e+07, 2.28272e+07, 2.33704e+07, 2.39136e+07, 2.44568e+07, 2.5e+07});
    //std::vector<double> ys({0, 0.0935744, 0.188578, 0.239477, 0.291646, 0.344675, 0.398283, 0.452277, 0.506536, 0.560964, 0.615471, 0.670006, 0.724548, 0.778992, 0.833378, 0.887649, 0.941796, 0.995794, 1.04965, 1.10332, 1.15682, 1.21014, 1.26326, 1.31619, 1.36894, 1.42145, 1.47379, 1.52593, 1.57785, 1.62955, 1.68108, 1.73237, 1.78347, 1.83436, 1.88503, 1.93555, 1.98579, 2.08574, 2.18488, 2.28324, 2.38082, 2.47764, 2.57371, 2.66905, 2.76367, 2.85758, 2.95079, 3.04331, 3.13519, 3.22638, 3.31701, 3.40697, 3.49631, 3.58507, 3.67318, 3.76075, 3.84775, 3.93419, 4.02008, 4.10543, 4.19026, 4.27459, 4.35839, 4.44168, 4.5245, 4.60684, 4.68873, 4.77011, 4.85103, 4.93156, 5.01159, 5.09122, 5.17041, 5.24917, 5.32755, 5.40549, 5.48306, 5.56022, 5.63701, 5.71337, 5.78938, 5.86503, 5.94031, 6.01522, 6.08979, 6.16402, 6.23788, 6.3114, 6.38459, 6.45747, 6.53004, 6.60226, 6.67419, 6.74577, 6.8171, 6.88808, 6.95879, 7.02916, 7.09927, 7.16904, 7.23856, 7.30779, 7.37676, 7.4454, 7.5138, 7.58196, 7.64985, 7.71746, 7.78485, 7.85195, 7.91882, 7.98543, 8.05181, 8.11793, 8.18377, 8.2494, 8.3148, 8.37998, 8.44489, 8.50958, 8.57407, 8.63833, 8.70236, 8.76619, 8.82976, 8.89314, 8.95633, 9.01921, 9.08198, 9.14453, 9.20687, 9.26901, 9.33096, 9.39269, 9.45423, 9.51552, 9.57669, 9.63761, 9.69837, 9.759, 9.8194, 9.8796, 9.93964, 9.99949, 10.0592, 10.1187, 10.178, 10.2372, 10.2962, 10.355, 10.4136, 10.472, 10.5303, 10.5885, 10.6465, 10.7043, 10.7619, 10.8195, 10.8768, 10.934, 10.991, 11.0479, 11.1046, 11.1612, 11.2176, 11.3299, 11.4417, 11.5528, 11.6634, 11.7735, 11.883, 11.992, 12.1004, 12.2083, 12.3157, 12.4225, 12.5289, 12.6348, 12.7401, 12.8451, 12.9494, 13.0533, 13.1567, 13.2597, 13.3623, 13.4644, 13.566, 13.6672, 13.7679, 13.8683, 13.9682, 14.0677, 14.1668, 14.2655, 14.3637, 14.4616, 14.5591, 14.6562, 14.7529, 14.8492, 14.9452, 15.0408, 15.136, 15.2309, 15.3254, 15.4196, 15.5134, 15.6064, 15.7, 15.7928, 15.8853, 15.9774, 16.0692, 16.1606, 16.2517, 16.3426, 16.4331, 16.5233, 16.6131, 16.7026, 16.7919, 16.8809, 16.9695, 17.0579, 17.146, 17.2338, 17.3213, 17.4955, 17.6685, 17.8405, 18.0115, 18.1814, 18.3503, 18.5181, 18.6849, 18.8505, 19.0154, 19.1793, 19.3423, 19.5043, 19.6635, 19.8258, 19.9853, 20.1437, 20.3014, 20.4582, 20.6144, 20.7696, 20.924, 21.0776, 21.2303, 21.3817, 21.5337, 21.6841, 21.834, 21.9831, 22.1315, 22.2792, 22.4262, 22.5725, 22.7182, 22.863, 23.0073, 23.151, 23.294, 23.4364, 23.5781, 23.719, 23.8595, 23.9994, 24.1386, 24.2776, 24.4156, 24.5531, 24.6901, 24.8266, 24.9625, 25.0978, 25.233, 25.3673, 25.501, 25.6342, 25.7668, 25.8953, 26.0306, 26.2923, 26.5522, 26.81, 27.0662, 27.3207, 27.5731, 27.8237, 28.0725, 28.3195, 28.5649, 28.8088, 29.0509, 29.2916, 29.5306, 29.7682, 30.0043, 30.2389, 30.472, 30.7036, 30.9339, 31.1627, 31.3902, 31.6163, 31.8412, 32.0648, 32.2871, 32.5082, 32.728, 32.9465, 33.164, 33.3803, 33.5953, 33.8093, 34.0219, 34.2336, 34.4442, 34.6539, 34.8623, 35.0697, 35.2761, 35.4814, 35.6857, 35.8889, 36.0913, 36.2927, 36.4931, 36.6927, 36.8912, 37.0896, 37.2855, 37.4813, 37.6763, 37.8703, 38.0635, 38.256, 38.6383, 39.0172, 39.3929, 39.7655, 40.135, 40.5015, 40.865, 41.2257, 41.5835, 41.9386, 42.2908, 42.6407, 42.9878, 43.3325, 43.6746, 44.0136, 44.3508, 44.6858, 45.0184, 45.3494, 45.6776, 46.0035, 46.3273, 46.6491, 46.9687, 47.2862, 47.6017, 47.9152, 48.2269, 48.5367, 48.845, 49.151, 49.4553, 49.7577, 50.0586, 50.3576, 50.6544, 50.9501, 51.244, 51.5364, 51.8271, 52.1151, 52.404, 52.69, 52.9746, 53.2584, 53.5392, 53.8193, 54.098, 54.3751, 54.651, 55.1987, 55.7411, 56.2781, 56.8102, 57.3372, 57.8595, 58.3771, 58.89, 59.3983, 59.9023, 60.402, 60.8974, 61.3887, 61.8761, 62.3595, 62.8392, 63.315, 63.7871, 64.2556, 64.7202, 65.1818, 65.6392, 66.094, 66.5454, 66.9944, 67.4393, 67.8811, 68.3199, 68.7558, 69.1888, 69.6188, 70.0461, 70.4706, 70.8923, 71.3115, 71.7278, 72.1415, 72.5528, 72.9614, 73.3675, 73.7679, 74.0824, 74.5719, 74.9685, 75.3629, 75.755, 76.1448, 76.5324, 76.918, 77.682, 78.439, 79.1876, 79.9279, 80.661, 81.3864, 82.1047, 82.8156, 83.5201, 84.2178, 84.9089, 85.5936, 86.2722, 86.9445, 87.6111, 88.2717, 88.9267, 89.576, 90.22, 90.8587, 91.4924, 92.1207, 92.7441, 93.3627, 93.9764, 94.5857, 95.19, 95.79, 96.3855, 96.9767, 97.5627, 98.1454, 98.7241, 99.2995, 99.8702, 100.437, 101, 101.559, 102.114, 102.666, 103.215, 103.759, 104.3, 104.838, 105.373, 105.904, 106.957, 107.997, 109.025, 110.041, 111.046, 112.04, 113.023, 113.996, 114.958, 115.909, 116.851, 117.784, 118.707, 119.622, 120.528, 121.425, 122.313, 123.194, 124.067, 124.932, 125.79, 126.639, 127.482, 128.318, 129.146, 129.968, 130.783, 131.591, 132.393, 133.189, 133.978, 134.762, 135.539, 136.311, 137.077, 137.837, 138.592, 139.342, 140.087, 140.824, 141.558, 142.287, 143.011, 144.445, 145.86, 147.256, 148.635, 149.997, 151.343, 152.673, 153.987, 155.286, 156.569, 157.839, 159.095, 160.337, 161.566, 162.783, 163.989, 165.179, 166.359, 167.528, 168.685, 169.832, 170.968, 172.093, 173.21, 174.316, 175.408, 176.494, 177.571, 178.639, 179.698, 180.747, 181.788, 182.821, 183.845, 184.861, 185.869, 186.87, 187.863, 188.849, 189.827, 190.798, 191.761, 193.668, 195.548, 197.402, 199.231, 201.037, 202.818, 204.577, 206.315, 208.03, 209.726, 211.406, 213.057, 214.695, 216.313, 217.913, 219.496, 221.064, 222.613, 224.147, 225.666, 227.169, 228.658, 230.133, 231.592, 233.038, 234.471, 235.893, 237.3, 238.691, 240.072, 241.442, 242.8, 244.147, 245.485, 246.809, 248.112, 249.421, 250.713, 251.995, 253.267, 255.782, 258.259, 260.711, 263.11, 265.474, 267.813, 270.12, 272.396, 274.643, 276.861, 279.051, 281.215, 283.354, 285.466, 287.553, 289.617, 291.657, 293.675, 295.67, 297.643, 299.597, 301.531, 303.444, 305.338, 307.215, 309.072, 310.911, 312.733, 314.537, 316.325, 318.096, 319.851, 321.591, 323.315, 325.024, 326.719, 328.399, 330.065, 333.363, 336.596, 339.784, 342.932, 346.016, 349.064, 352.067, 355.028, 357.948, 360.828, 363.67, 366.488, 369.246, 371.982, 374.684, 377.353, 379.99, 382.596, 385.173, 387.721, 390.24, 392.732, 395.199, 397.638, 400.39, 403.11, 405.798, 408.457, 411.086, 413.685, 416.257, 418.802, 421.32, 423.812, 426.279, 428.72, 431.139, 433.533, 435.905, 438.255, 440.583, 442.888, 445.173, 447.447, 451.907, 456.301, 460.622, 464.871, 469.052, 473.168, 477.222, 481.216, 485.149, 489.027, 492.87, 496.648, 500.343, 504.014, 507.636, 511.213, 514.746, 518.234, 521.683, 525.089, 528.456, 531.783, 535.074, 538.327, 541.545, 544.725, 547.875, 550.992, 554.076, 557.13, 560.153, 563.15, 566.114, 569.05, 571.957, 574.838, 580.521, 586.104, 591.588, 596.98, 602.282, 607.499, 612.632, 617.687, 622.665});
    //std::variant<linear_interpolator, xsectval> sigma_jets = linear_interpolator(xs, ys);
    //Only one sigma_jet
    //std::variant<linear_interpolator, xsectval> sigma_jets = xsectval(179.2);
    //std::variant<linear_interpolator, xsectval> sigma_jets = pqcd::calculate_sigma_jet(p_pdf, &mand_s, &kt02, &jet_params);

    InterpMultilinear<4, xsectval> sigma_jets = read_sigma_jets("sigma_jet_grid (copy).dat");
    
    if (!read_sigmajets_from_file)
    {
        double tolerance=0.01, upper_tAA_0_limit=40.0, lower_tAA_0_limit = 35.0, upper_sumTpp_limit=0.40, lower_sumTpp_limit=0.30;
        //double tolerance=0.05, upper_tAA_0_limit=42.0, lower_tAA_0_limit = 33.0, upper_sumTpp_limit=0.5, lower_sumTpp_limit=0.03411;
        InterpMultilinear<4, xsectval> sigma_jets = calculate_spatial_sigma_jets(tolerance, p_pdf, p_pdf, mand_s, kt02, jet_params, upper_tAA_0_limit, lower_tAA_0_limit, upper_sumTpp_limit, lower_sumTpp_limit);
        //array<spatial,4> args{0.1, 0.3, 35.0, 40.0};
        //std::cout<<sigma_jets.interp(args.begin())<<std::endl;
        //std::cout<<pqcd::calculate_spatial_sigma_jet_mf(p_pdf, p_pdf, &mand_s, &kt02, &jet_params, &args[0], &args[1], &args[2], &args[3])<<std::endl;
    }    

    if (verbose) std::cout<<"Done!"<<std::endl;
    
    do // while (AA_events < desired_N_events);
    {
        spatial impact_parameter = sqrt(b_min*b_min + unirand(*eng)*(b_max*b_max-b_min*b_min)); //B^2 from a uniform distribution
        if (verbose) std::cout<<"impact_parameter: "<<impact_parameter<<std::endl;

        auto [pro, tar] = generate_nuclei(nuc_params, sqrt_s, impact_parameter, eng, radial_sampler, read_nuclei_from_file, verbose);        
        //collide_nuclei(pro, tar, binary_collisions, sigma_jets, unirand, eng, coll_params, verbose);
        collide_nuclei_with_spatial_pdfs(pro, tar, binary_collisions, sigma_jets, unirand, eng, coll_params, verbose, Tpp);

        nof_collisions++;
        ulong NColl=binary_collisions.size();

        if (NColl<1) 
        {
            binary_collisions.erase(binary_collisions.begin(), binary_collisions.end());
        
            if (verbose || (nof_collisions%1000)==0)
            {
                std::cout << std::endl << "A+A collided thus far: " << nof_collisions << ", of which events thus far: " << AA_events << std::endl << std::endl;
            }
            continue;
        }
        AA_events++;

        if (end_state_filtering)
        {
            filter_end_state(binary_collisions, filtered_scatterings);
        }
        
        binary_collisions.erase(binary_collisions.begin(), binary_collisions.end());
        
        if (verbose || (nof_collisions%1000)==0)
        {
            std::cout << std::endl << "A+A collided thus far: " << nof_collisions << ", of which events thus far: " << AA_events << std::endl << std::endl;
        }

        uint Npart=0;
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
        collisions_for_reporting.push_back(coll);
        
        filtered_scatterings.clear();
        
    } while (AA_events < desired_N_events);

    std::cout << collisions_for_reporting.size() << " collisions generated" << std::endl;

    uint nBins = 16;
    double binsLow[] = {0., 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 0.0, 0.0};
    double binsHigh[] = {0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 0.8, 1.0};
    mc_glauber_style_report(collisions_for_reporting, sigma_inel_for_glauber, desired_N_events, nBins, binsLow, binsHigh);

    return 0;
}
