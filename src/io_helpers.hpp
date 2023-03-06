//Copyright (c) 2022 Mikko Kuha

#ifndef IO_HELPERS_HPP
#define IO_HELPERS_HPP

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "histo.hpp"
#include "nucleon.hpp"
#include "typedefs.hpp"

class io
{
public:
    class Coll 
    {

        private:
        double b, x, nTwo, sumET;
        ulong nColl, nPart, nEnd;
        bool alice = true;

        public:
        Coll() { b = 0.; x = 0.; nTwo = 0.; nColl = 0; nPart = 0; nEnd = 0; sumET = 0.;}
        Coll(ulong nColl_, ulong nPart_, ulong nEnd_, double b_, double sumET_){
            b = b_; nColl = nColl_; nPart = nPart_; nEnd = nEnd_; sumET = sumET_;
            // ATLAS definitions. 
            if (!alice) {
            x = 0.09;
            nTwo = 0.5 * (1. - x) * static_cast<double>(nPart) + x * static_cast<double>(nColl);
            // ALICE definitions.
        } else {
            x = 0.801;
            nTwo = x * static_cast<double>(nPart) + (1. - x) * static_cast<double>(nColl);
        }
        }
        Coll(ulong nCollIn, ulong nPartIn, ulong nEndIn, double bIn, std::shared_ptr<std::mt19937> random_generator){
            b = bIn; nColl = nCollIn; nPart = nPartIn; nEnd = nEndIn; sumET = 0.;
            // ATLAS definitions. 
            if (!alice) {
            x = 0.09; 
            nTwo = 0.5 * (1. - x) * static_cast<double>(nPart) + x * static_cast<double>(nColl); sampleET(random_generator);
            // ALICE definitions.
        } else {
            x = 0.801;
            nTwo = x * static_cast<double>(nPart) + (1. - x) * static_cast<double>(nColl); sampleET(random_generator);
        }
        }
        double getB() const {return b;}
        ulong getNcoll() const {return nColl;}
        ulong getNpart() const {return nPart;}
        ulong getNend() const {return nEnd;}
        double getNtwo() const {return nTwo;}
        double getET() const {return sumET;}
        void sampleET(std::shared_ptr<std::mt19937> random_generator) {
            double mu = 46.4;
            double k = 1.5;
            double p = (mu/k) / (mu/k + 1.);
            uint_fast32_t nAncestors = static_cast<uint_fast32_t>(floor(nTwo));
            uint_fast32_t mult = 0;
            //std::cout << (*random_generator)() << ' ' << std::flush;
            std::negative_binomial_distribution<uint_fast32_t> distr(static_cast<uint_fast32_t>(k),1.-p);
            for (uint_fast32_t i = 0; i < nAncestors; ++i) {
            mult += distr(*random_generator);
            //std::cout << mult << ' ' << std::flush;
            }
            sumET = static_cast<double>(mult);
            //std::cout << std::endl << nAncestors << ' ' << sumET/nAncestors << std::endl << std::endl;
            // cout << endl << nAncestors << ' ' << sumET/nAncestors << endl << endl;
        }
    };
    static bool compNpart(Coll c1, Coll c2){
        return c1.getNpart() < c2.getNpart();
        }
    static bool compNcoll(Coll c1, Coll c2){
        return c1.getNcoll() < c2.getNcoll();
        }
    static bool compNend(Coll c1, Coll c2){
        return c1.getNend() < c2.getNend();
        }
    static bool compNtwo(Coll c1, Coll c2){
        return c1.getNtwo() < c2.getNtwo();
        }
    static bool compET(Coll c1, Coll c2){
        return c1.getET() < c2.getET();
        }
    template <typename T>
    static double calc_ave(std::vector<Coll> &collisions, T (Coll::*func)()const ){
        std::vector<Coll>::iterator it;
        T tot = 0;
        for (it = collisions.begin(); it != collisions.end(); ++it){
            tot += ((*it).*func)();
        }
        return static_cast<double>(tot)/static_cast<double>(collisions.size());  
    }


    static auto mc_glauber_style_report
    (
        std::vector<Coll> &collisions, 
        const double &sigma_inel, 
        const uint_fast64_t &N_events, 
        const uint_fast64_t &nBins, 
        const double *const binsLow, 
        const double *const binsHigh,
        std::ostream &out_stream
    ) noexcept -> void
    {
        // Make sure that no rounding downwards.
        double eps = 0.1/static_cast<double>(N_events);
        
        // int comp = 1;
        for (uint_fast8_t comp = 1; comp <= 5; ++comp) {
            
            // Print out the header.
            out_stream << std::endl;
            out_stream << "# sigma_inel = " << sigma_inel << " mb" << std::endl;
            
            switch (comp)
            {
                case 1:
                    std::sort(collisions.begin(), collisions.end(), compNpart);
                    out_stream << "Using Npart for centrality determination" 
                            << std::endl << std::endl
                            << "#  cent     nPmax     nPmin       <b>   <nPart>   <nColl>    <nEnd>    <T_AA>     <E_T>"
                            << std::endl;
                    break;
                case 2:
                    std::sort(collisions.begin(), collisions.end(), compNcoll);
                    out_stream << "Using Ncoll for centrality determination" 
                            << std::endl << std::endl
                            << "#  cent     nCmax     nCmin       <b>   <nPart>   <nColl>    <nEnd>    <T_AA>     <E_T>"
                            << std::endl;
                    break;
                case 3:
                    std::sort(collisions.begin(), collisions.end(), compNtwo);
                    out_stream << "Using two-component (ancestors) model for centrality determination"
                            << std::endl << std::endl
                            << "#  cent     nAmax     nAmin       <b>   <nPart>   <nColl>    <nEnd>    <T_AA>     <E_T>"
                            << std::endl;
                    break;
                case 4:
                    std::sort(collisions.begin(), collisions.end(), compET);
                    out_stream << "Using sumET model for centrality determination" 
                            << std::endl << std::endl
                            << "#  cent     ETmax     ETmin       <b>   <nPart>   <nColl>    <nEnd>    <T_AA>     <E_T>"
                            << std::endl;
                    break;
                default: //case 5
                    std::sort(collisions.begin(), collisions.end(), compET);
                    out_stream << "Using the number of produced jets for centrality determination" 
                            << std::endl << std::endl
                            << "#  cent     nJmax     nJmin       <b>   <nPart>   <nColl>    <nEnd>    <T_AA>     <E_T>"
                            << std::endl;
            }
            std::reverse(collisions.begin(), collisions.end());
            
            std::vector<Coll> centrality;
            // Number of collisions within centrality class i.
            for (uint_fast8_t i = 0; i < nBins; i++)
            {
                uint_fast64_t lower = static_cast<uint_fast64_t>(binsLow[i]*static_cast<double>(N_events)+eps);
                uint_fast64_t upper = static_cast<uint_fast64_t>(binsHigh[i]*static_cast<double>(N_events)+eps);
                centrality.clear();
                while (lower<upper)
                {
                    centrality.push_back(collisions.at(lower++));
                }
                
                // Print the centrality selection.
                out_stream << std::setw(2) << binsLow[i]*100 << " " << std::setw(4) << binsHigh[i]*100;
                
                switch (comp)
                {
                    case 1:
                        out_stream << std::setw(10) << centrality[0].getNpart()
                                   << std::setw(10) << centrality[centrality.size() - 1].getNpart();
                        break;
                    case 2:
                        out_stream << std::setw(10) << centrality[0].getNcoll() 
                                   << std::setw(10) << centrality[centrality.size() - 1].getNcoll();
                        break;
                    case 3:
                        out_stream << std::setw(10) <<centrality[0].getNtwo()
                                   << std::setw(10) << centrality[centrality.size() - 1].getNtwo();
                        break;
                    case 4:
                        out_stream << std::setw(10) << centrality[0].getET() 
                                   << std::setw(10) << centrality[centrality.size() - 1].getET();
                        break;
                    default : //case 5
                        out_stream << std::setw(10) << centrality[0].getNend() 
                                   << std::setw(10) << centrality[centrality.size() - 1].getNend();
                        break;
                }
                
                out_stream << std::setw(10) << io::calc_ave<double>(centrality, &Coll::getB) 
                           << std::setw(10) << io::calc_ave<ulong>(centrality, &Coll::getNpart) 
                           << std::setw(10) << io::calc_ave<ulong>(centrality, &Coll::getNcoll)
                           << std::setw(10) << io::calc_ave<ulong>(centrality, &Coll::getNend) 
                           << std::setw(10) << io::calc_ave<ulong>(centrality, &Coll::getNcoll)/sigma_inel
                           << std::setw(10) << calc_ave<double>(centrality, &Coll::getET)
                           << std::endl;      
            }
        }
    }

    static auto print_3d_histo
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

    static auto print_2d_histo
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

        for (uint_fast8_t i = 0; i < n_x_bins; i++)
        {
            file << xs0[i] << ' ' << xs0[i+1] << ' ';

            for (uint_fast8_t j = 0; j < n_y_bins; j++)
            {   
                file << histo[i][j] * norm << ' ';
            }
            file << std::endl;
        }
        std::cout<<"printed to "<<filename<<std::endl;
        file.close();
    }

    static auto save_event
    (
        std::ofstream &event_file, 
        const std::vector<nucleon> &pro, 
        const std::vector<nucleon> &tar, 
        const double &impact_parameter
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

    static auto save_event
    (
        std::ofstream &event_file, 
        const std::vector<nucleon> &pro, 
        const std::vector<nucleon> &tar, 
        const double &impact_parameter,
        const std::vector<dijet_with_coords> &filtered_scatterings
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
        for (const auto & scc : filtered_scatterings)
        {
            const auto sc = scc.dijet;

            event_file << "        {"<<std::endl;
            event_file << "            "<<sc.kt<<','<<std::endl;
            event_file << "            "<<sc.y1<<','<<std::endl;
            event_file << "            "<<sc.y2<<','<<std::endl;
            event_file << "            "<<sc.init1<<','<<std::endl;
            event_file << "            "<<sc.init2<<','<<std::endl;
            event_file << "            "<<sc.final1<<','<<std::endl;
            event_file << "            "<<sc.final2<<std::endl;
            if (&scc != &filtered_scatterings.back())
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

/*
    static auto save_single_coll_csv
    (
        const std::array<std::vector<dijet_with_coords>, 4> &filtered_scatterings,
        std::uniform_real_distribution<double> &unirand,
        std::shared_ptr<std::mt19937> random_generator,
        bool include_TATA = false
    ) noexcept -> void
    {
        std::ofstream event_file;
        std::string name = "jets";
        std::string separator = ",";
        std::array<std::string, 4> names{name+".csv",
                                         name+"_MC.csv",
                                         name+"_SAT.csv",
                                         name+"_SAT_MC.csv"};
        
        double pt, energy, px, py, pz;
        double y1, y2;
        double theta;

        for (size_t i=0; i<4; i++)
        {
            event_file.open(names[i]);
            //Output header
            event_file << "t0" << separator;
            event_file << "x"   << separator;
            event_file << "y"   << separator;
            event_file << "z"   << separator;
            event_file << "E"   << separator;
            event_file << "p_x" << separator;
            event_file << "p_y" << separator;
            event_file << "p_z";
            if(include_TATA)
            { 
                event_file << separator;
                event_file << "p_T" << separator;
                event_file << "T_AT_A" << std::endl;
            }
            else
            {
                event_file << std::endl;
            }
            event_file << std::setprecision(std::numeric_limits<double>::digits10 + 1);

            for (auto e_co : filtered_scatterings[i])
            {
                pt = e_co.dijet.kt;
                y1 = e_co.dijet.y1;
                y2 = e_co.dijet.y2;
                // jet 1
                energy = pt*std::cosh(y1);
                theta = 2*M_PI*unirand(*random_generator);
                px = pt*std::cos(theta);
                py = pt*std::sin(theta);
                pz = pt*std::sinh(y1);
                event_file << e_co.t0  << separator; //t0
                event_file << e_co.co.x << separator; //x
                event_file << e_co.co.y << separator; //y
                event_file << e_co.co.z << separator; //z
                event_file << energy    << separator; //E
                event_file << px        << separator; //p_x
                event_file << py        << separator; //p_y
                event_file << pz                    ; //p_z
                if(include_TATA)
                { 
                    event_file << separator;
                    event_file << pt << separator;        //p_T
                    event_file << e_co.tata << std::endl; //T_A * T_A
                }
                else
                {
                    event_file << std::endl;
                }
                auto m = energy*energy-(px*px+py*py+pz*pz);
		        if (m >= 0.0001)
                {
                    std::cout<<"1: "<<m<<std::endl;
                }
		

                // jet 2
                energy = pt*std::cosh(y2);
                pz = pt*std::sinh(y2);
                event_file << e_co.t0  << separator; //t0
                event_file << e_co.co.x << separator; //x
                event_file << e_co.co.y << separator; //y
                event_file << e_co.co.z << separator; //z
                event_file << energy    << separator; //E
                event_file << -px       << separator; //p_x
                event_file << -py       << separator; //p_y
                event_file << pz                    ; //p_z
                if(include_TATA)
                { 
                    event_file << separator;
                    event_file << pt << separator;        //p_T
                    event_file << e_co.tata << std::endl; //T_A * T_A
                }
                else
                {
                    event_file << std::endl;
                }
                m = energy*energy-(px*px+py*py+pz*pz);
		        if (m >= 0.0001)
                {
                    std::cout<<"2: "<<m<<std::endl;
                }
            }
            event_file.close();
        }
    }


    static auto save_single_coll_binary
    (
        const std::array<std::vector<dijet_with_coords>, 4> &filtered_scatterings,
        std::uniform_real_distribution<double> &unirand,
        std::shared_ptr<std::mt19937> random_generator
    ) noexcept -> void
    {
        std::ofstream jet_file;
        std::string name = "jets";
        std::array<std::string, 4> names{name+".dat",
                                         name+"_MC.dat",
                                         name+"_SAT.dat",
                                         name+"_SAT_MC.dat"};
        
        double pt, y1, y2, t0, x, y, z, energy, px, py, pz, tata, theta;
        uint_fast64_t n_jets;

        for (uint_fast8_t i=0; i<4; i++)
        {
            n_jets = filtered_scatterings[i].size()*2;

            jet_file.open(names[i], std::ios::out | std::ios::binary);
            jet_file.write(reinterpret_cast<char*>(&n_jets), sizeof n_jets); //total number of jets

            for (auto e_co : filtered_scatterings[i])
            {
                pt = e_co.dijet.kt;
                y1 = e_co.dijet.y1;
                y2 = e_co.dijet.y2;
                t0 = e_co.t0;
                x = e_co.co.x;
                y = e_co.co.y;
                z = e_co.co.z;
                tata = e_co.tata;
                
                // jet 1
                energy = pt*std::cosh(y1);
                theta = 2*M_PI*unirand(*random_generator);
                px = pt*std::cos(theta);
                py = pt*std::sin(theta);
                pz = pt*std::sinh(y1);
                jet_file.write(reinterpret_cast<char*>(&t0)     , sizeof t0);     //t0
                jet_file.write(reinterpret_cast<char*>(&x)      , sizeof x);      //x
                jet_file.write(reinterpret_cast<char*>(&y)      , sizeof y);      //y
                jet_file.write(reinterpret_cast<char*>(&z)      , sizeof z);      //z
                jet_file.write(reinterpret_cast<char*>(&energy) , sizeof energy); //E
                jet_file.write(reinterpret_cast<char*>(&px)     , sizeof px);     //p_x
                jet_file.write(reinterpret_cast<char*>(&py)     , sizeof py);     //p_y
                jet_file.write(reinterpret_cast<char*>(&pz)     , sizeof pz);     //p_z
                jet_file.write(reinterpret_cast<char*>(&tata)   , sizeof tata);   //T_A * T_A

                // jet 2
                energy = pt*std::cosh(y2);
                px = -px;
                py = -py;
                pz = pt*std::sinh(y2);
                jet_file.write(reinterpret_cast<char*>(&t0)     , sizeof t0);     //t0
                jet_file.write(reinterpret_cast<char*>(&x)      , sizeof x);      //x
                jet_file.write(reinterpret_cast<char*>(&y)      , sizeof y);      //y
                jet_file.write(reinterpret_cast<char*>(&z)      , sizeof z);      //z
                jet_file.write(reinterpret_cast<char*>(&energy) , sizeof energy); //E
                jet_file.write(reinterpret_cast<char*>(&px)     , sizeof px);     //p_x
                jet_file.write(reinterpret_cast<char*>(&py)     , sizeof py);     //p_y
                jet_file.write(reinterpret_cast<char*>(&pz)     , sizeof pz);     //p_z
                jet_file.write(reinterpret_cast<char*>(&tata)   , sizeof tata);   //T_A * T_A
            }
            jet_file.close();
        }
    }
*/

    static auto append_single_coll_binary
    (
        std::ofstream &jet_file,
        const std::vector<dijet_with_coords> &filtered_scatterings,
        std::uniform_real_distribution<double> &unirand,
        std::shared_ptr<std::mt19937> random_generator
    ) noexcept -> void
    {        
        //double pt, y1, y2, t0, x, y, z, energy, px, py, pz, tata, phi;
        double pt, y1, y2, t01, t02, x, y, tata, phi;
        uint_fast64_t n_dijets;
        n_dijets = filtered_scatterings.size();
        jet_file.write(reinterpret_cast<char*>(&n_dijets), sizeof n_dijets); //total number of dijets

        for (auto e_co : filtered_scatterings)
        {
            pt = e_co.dijet.kt;
            y1 = e_co.dijet.y1;
            y2 = e_co.dijet.y2;
            t01 = e_co.t01;
            t02 = e_co.t02;
            x = e_co.co.x;
            y = e_co.co.y;
            //z = e_co.co.z;
            tata = e_co.tata;
            
            // jet 1
            //energy = pt*std::cosh(y1);
            phi = 2*M_PI*unirand(*random_generator);
            //px = pt*std::cos(theta);
            //py = pt*std::sin(theta);
            //pz = pt*std::sinh(y1);
            jet_file.write(reinterpret_cast<char*>(&t01)    , sizeof t01);    //t01
            jet_file.write(reinterpret_cast<char*>(&t02)    , sizeof t02);    //t02
            jet_file.write(reinterpret_cast<char*>(&x)      , sizeof x);      //x
            jet_file.write(reinterpret_cast<char*>(&y)      , sizeof y);      //y
            jet_file.write(reinterpret_cast<char*>(&pt)     , sizeof pt);     //p_T
            jet_file.write(reinterpret_cast<char*>(&y1)     , sizeof y1);     //y1
            jet_file.write(reinterpret_cast<char*>(&y2)     , sizeof y2);     //y2
            jet_file.write(reinterpret_cast<char*>(&phi)    , sizeof phi);    //phi
            jet_file.write(reinterpret_cast<char*>(&tata)   , sizeof tata);   //T_A * T_A

            // jet 2
            //energy = pt*std::cosh(y2);
            //px = -px;
            //py = -py;
            //pz = pt*std::sinh(y2);
            //jet_file.write(reinterpret_cast<char*>(&t0)     , sizeof t0);     //t0
            //jet_file.write(reinterpret_cast<char*>(&x)      , sizeof x);      //x
            //jet_file.write(reinterpret_cast<char*>(&y)      , sizeof y);      //y
            //jet_file.write(reinterpret_cast<char*>(&z)      , sizeof z);      //z
            //jet_file.write(reinterpret_cast<char*>(&energy) , sizeof energy); //E
            //jet_file.write(reinterpret_cast<char*>(&px)     , sizeof px);     //p_x
            //jet_file.write(reinterpret_cast<char*>(&py)     , sizeof py);     //p_y
            //jet_file.write(reinterpret_cast<char*>(&pz)     , sizeof pz);     //p_z
            //jet_file.write(reinterpret_cast<char*>(&tata)   , sizeof tata);   //T_A * T_A
        }
    }

    static auto print_1d_histo
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

    static auto read_conf
    (
        const std::string filename
    ) noexcept
    {
        std::ifstream input(filename);
        
        std::string name{"example_name"};
        std::string sigmajet_filename{"sigma_jet.dat"}; 
        uint_fast32_t n_events{500};
        double b_max{20.0};
        double b_min{0.0};
        double sqrt_s{5020.0};
        double K_factor{1.0};
        double kt0{1.0};
        double proton_width{0.573};
        double sigma_inel{70.0};
        double sigma_inel_AA{70.0};
        double T_AA_0_for_snpdfs{0.0};
        double spatial_cutoff{0.0};
        uint_fast16_t A{208};
        uint_fast16_t ZA{82};
        double M_factor{2.5};
        double nn_min_dist{0.4};
        double nuclear_RA{6.62435};
        double nuclear_d{0.5498};
        double rad_max{30.0};
        double rad_min{0.0};
        bool is_AA{true};
        bool is_pA{false};
        bool is_pp{false};
        bool is_mc_glauber{true};
        bool read_sigmajets_from_file{false};
        bool proton_width_static{false};
        bool sigma_inel_from_sigma_jet{true};
        bool AA_inel_same_as_NN{false};
        bool only_protons{false};
        bool use_npdfs{true};
        bool use_snpdfs{false}; 
        bool snpdfs_linear{false};
        bool calculate_spatial_cutoff{true};
        bool calculate_end_state{true};
        bool calculate_tata{true};
        bool save_endstate_jets{true};
        bool end_state_filtering{true};
        bool is_mom_cons{true};
        bool is_saturation{true};

        uint_fast16_t count = 0;

        if (input.is_open())
        {
            std::string line;

            for (; std::getline(input, line); )
            {
                std::istringstream line_stream(line);
                std::string param_name;
                line_stream >> param_name;

                if (param_name == "name")
                {
                    line_stream >> name;
                }
                else if (param_name == "sigmajet_filename")
                {
                    line_stream >> sigmajet_filename;
                }
                else if (param_name == "n_events")
                {
                    line_stream >> n_events;
                }
                else if (param_name == "b_max")
                {
                    line_stream >> b_max;
                }
                else if (param_name == "b_min")
                {
                    line_stream >> b_min;
                }
                else if (param_name == "sqrt_s")
                {
                    line_stream >> sqrt_s;
                }
                else if (param_name == "K_factor")
                {
                    line_stream >> K_factor;
                }
                else if (param_name == "kt0")
                {
                    line_stream >> kt0;
                }
                else if (param_name == "proton_width")
                {
                    line_stream >> proton_width;
                }
                else if (param_name == "sigma_inel")
                {
                    line_stream >> sigma_inel;
                }
                else if (param_name == "sigma_inel_AA")
                {
                    line_stream >> sigma_inel_AA;
                }
                else if (param_name == "T_AA_0_for_snpdfs")
                {
                    line_stream >> T_AA_0_for_snpdfs;
                }
                else if (param_name == "spatial_cutoff")
                {
                    line_stream >> spatial_cutoff;
                }
                else if (param_name == "A")
                {
                    line_stream >> A;
                }
                else if (param_name == "ZA")
                {
                    line_stream >> ZA;
                }
                else if (param_name == "M_factor")
                {
                    line_stream >> M_factor;
                }
                else if (param_name == "nn_min_dist")
                {
                    line_stream >> nn_min_dist;
                }
                else if (param_name == "nuclear_RA")
                {
                    line_stream >> nuclear_RA;
                }
                else if (param_name == "nuclear_d")
                {
                    line_stream >> nuclear_d;
                }
                else if (param_name == "rad_max")
                {
                    line_stream >> rad_max;
                }
                else if (param_name == "rad_min")
                {
                    line_stream >> rad_min;
                }
                else if (param_name == "is_aa")
                {
                    line_stream >> std::boolalpha >> is_AA;
                }
                else if (param_name == "is_pa")
                {
                    line_stream >> std::boolalpha >> is_pA;
                }
                else if (param_name == "is_pp")
                {
                    line_stream >> std::boolalpha >> is_pp;
                }
                else if (param_name == "is_mc_glauber")
                {
                    line_stream >> std::boolalpha >> is_mc_glauber;
                }
                else if (param_name == "read_sigmajets_from_file")
                {
                    line_stream >> std::boolalpha >> read_sigmajets_from_file;
                }
                else if (param_name == "proton_width_static")
                {
                    line_stream >> std::boolalpha >> proton_width_static;
                }
                else if (param_name == "sigma_inel_from_sigma_jet")
                {
                    line_stream >> std::boolalpha >> sigma_inel_from_sigma_jet;
                }
                else if (param_name == "AA_inel_same_as_NN")
                {
                    line_stream >> std::boolalpha >> AA_inel_same_as_NN;
                }
                else if (param_name == "only_protons")
                {
                    line_stream >> std::boolalpha >> only_protons;
                }
                else if (param_name == "use_npdfs")
                {
                    line_stream >> std::boolalpha >> use_npdfs;
                }
                else if (param_name == "use_snpdfs")
                {
                    line_stream >> std::boolalpha >> use_snpdfs;
                }
                else if (param_name == "snpdfs_linear")
                {
                    line_stream >> std::boolalpha >> snpdfs_linear;
                }
                else if (param_name == "calculate_spatial_cutoff")
                {
                    line_stream >> std::boolalpha >> calculate_spatial_cutoff;
                }
                else if (param_name == "calculate_end_state")
                {
                    line_stream >> std::boolalpha >> calculate_end_state;
                }
                else if (param_name == "calculate_tata")
                {
                    line_stream >> std::boolalpha >> calculate_tata;
                }
                else if (param_name == "save_endstate_jets")
                {
                    line_stream >> std::boolalpha >> save_endstate_jets;
                }
                else if (param_name == "end_state_filtering")
                {
                    line_stream >> std::boolalpha >> end_state_filtering;
                }
                else if (param_name == "is_mom_cons")
                {
                    line_stream >> std::boolalpha >> is_mom_cons;
                }
                else if (param_name == "is_saturation")
                {
                    line_stream >> std::boolalpha >> is_saturation;
                }
                else
                {
                    continue;
                }
                count++;
            }
        }
        else
        {
            std::cout<<"Error opening the configuration file!!"<<std::endl;
            std::terminate();
        }

        std::cout<<"Read "<< count <<" parameters from a file "<< filename <<std::endl;

        return std::make_tuple
                (
                    name,
                    sigmajet_filename,
                    n_events,
                    b_max,
                    b_min,
                    sqrt_s,
                    K_factor,
                    kt0,
                    proton_width,
                    sigma_inel,
                    sigma_inel_AA,
                    T_AA_0_for_snpdfs,
                    spatial_cutoff,
                    A,
                    ZA,
                    M_factor,
                    nn_min_dist,
                    nuclear_RA,
                    nuclear_d,
                    rad_max,
                    rad_min,
                    is_AA,
                    is_pA,
                    is_pp,
                    is_mc_glauber,
                    read_sigmajets_from_file,
                    proton_width_static,
                    sigma_inel_from_sigma_jet,
                    AA_inel_same_as_NN,
                    only_protons,
                    use_npdfs,
                    use_snpdfs,
                    snpdfs_linear,
                    calculate_spatial_cutoff,
                    calculate_end_state,
                    calculate_tata,
                    save_endstate_jets,
                    end_state_filtering,
                    is_mom_cons,
                    is_saturation
                );
    }

    static auto print_histos
    (
        const std::string &name_postfix,
        const histo_2d &jets,
        const histo_2d &dijets,
        const histo_1d &dETdy,
        const histo_1d &dEdy,
        const histo_1d &dNdy,
        const histo_1d &dNdET,
        const histo_1d &dETdeta,
        const histo_1d &dEdeta,
        const histo_1d &dETdb,
        const histo_1d &dEdb,
        const double &dijet_norm,
        const uint_fast32_t &AA_events
    ) noexcept -> void
    {
        std::string true_postfix{name_postfix+".dat"};

        //sigma1jet
        print_2d_histo
        (
            jets, 
            "sigma1jet_sim_"+true_postfix, 
            2.0 * dijet_norm
        );

        //dNdpTdy
        print_2d_histo
        (
            jets,
            "dNdpTdy_sim_"+true_postfix, 
            1.0,
            false
        );

        //sigmadijet
        print_2d_histo
        (
            dijets,
            "sigmadijet_sim_"+true_postfix, 
            dijet_norm
        );

        //dETdy
        print_1d_histo
        (
            dETdy,
            "dETdy_sim_"+true_postfix, 
            1.0/ static_cast<double>(AA_events),
            false
        );

        //dEdy
        print_1d_histo
        (
            dEdy, 
            "dEdy_sim_"+true_postfix, 
            1.0/ static_cast<double>(AA_events),
            false
        );

        //dNdy
        print_1d_histo
        (
            dNdy, 
            "dNdy_sim_"+true_postfix, 
            1.0/ static_cast<double>(AA_events),
            false
        );

        //dNdET
        print_1d_histo
        (
            dNdET, 
            "dNdET_sim_"+true_postfix, 
            1.0,
            false
        );

        //dETdeta
        print_1d_histo
        (
            dETdeta, 
            "dETdeta_sim_"+true_postfix, 
            1.0 / static_cast<double>(AA_events),
            false
        );

        //dEdeta
        print_1d_histo
        (
            dEdeta, 
            "dEdeta_sim_"+true_postfix, 
            1.0 / static_cast<double>(AA_events),
            false
        );

        //dETdb
        print_1d_histo
        (
            dETdb, 
            "dETdb_sim_"+true_postfix, 
            1.0,
            true
        );

        //dEdb
        print_1d_histo
        (
            dEdb,
            "dEdb_sim_"+true_postfix, 
            1.0,
            true
        );
        
    }
};

#endif // IO_HELPERS_HPP