#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

struct dijet
{
    double t01;
    double t02;
    double x;
    double y;
    double pt;
    double y1;
    double y2;
    double phi;
    double tata;
    int b1;
    int q1;
    int b2;
    int q2;
    dijet(const double &t01_, const double &t02_, const double &x_, const double &y_, 
            const double &pt_, const double &y1_, const double &y2_, const double &phi_, 
            const double &tata_, const int_fast8_t &f1, const int_fast8_t &f2)
    : t01(t01_), t02(t02_), x(x_), y(y_), pt(pt_), y1(y1_), y2(y2_), phi(phi_), tata(tata_)
    {
        if (f1 != 0)
        {
            b1 = static_cast<int>(f1 / abs(f1));
            if (f1 % 2 == 0)
            {
                q1 = 2*b1;
            }
            else
            {
                q1 = -1*b1;
            }
        }
        else
        {
            b1 = 0;
            q1 = 0;
        }
        if (f2 != 0)
        {
            b2 = static_cast<int>(f2 / abs(f2));
            if (f2 % 2 == 0)
            {
                q2 = 2*b2;
            }
            else
            {
                q2 = -1*b2;
            }
        }
        else
        {
            b2 = 0;
            q2 = 0;
        }
    }
};

 
int main()
{
    std::string filename = "jets_0_100_LHC_ref_flav.dat";
    std::ifstream jet_file(filename, std::ios::in | std::ios::binary);

    if (!jet_file.is_open())
    {
        std::cout << "Error opening file " << filename << std::endl;
        exit(1);
    }

    uint_fast64_t n_events, n_jets;
    double t01, t02, x, y, pt, y1, y2, phi, tata;
    int_fast8_t f1, f2;
    int quarks, charges;
    double total_B = 0;
    double total_Q = 0;
    std::vector<std::vector<dijet> > events_jets;
    
    if (sizeof n_jets != 8 || sizeof t01 != 8 )
    {
        std::cout << "Change the types! n_jets should be 64 bit unsigned integer"
                  << " and the coordinates should all be 64 bit float numbers" << std::endl;
        exit(1);
    }

    jet_file.read(reinterpret_cast<char*>(&n_events), sizeof n_events);
    std::cout << "Trying to read " << n_events << " events from the file " 
              << filename << " ..." << std::endl;
    events_jets.reserve(n_events);

    uint_fast16_t n_bins = 100;

    uint_fast64_t n_dijet_total = 0;
    for (uint_fast64_t i=0; i<n_events; i++)
    {
        jet_file.read(reinterpret_cast<char*>(&n_jets), sizeof n_jets);

        std::vector<dijet> jets;
        jets.reserve(n_jets);
        quarks = 0;
        charges = 0;

        for (uint_fast64_t ii=0; ii<n_jets; ii++, n_dijet_total++)
        {
            jet_file.read(reinterpret_cast<char*>(&t01)  , sizeof t01);  
            jet_file.read(reinterpret_cast<char*>(&t02)  , sizeof t02);  
            jet_file.read(reinterpret_cast<char*>(&x)    , sizeof x);    
            jet_file.read(reinterpret_cast<char*>(&y)    , sizeof y);    
            jet_file.read(reinterpret_cast<char*>(&pt)   , sizeof pt);   
            jet_file.read(reinterpret_cast<char*>(&y1)   , sizeof y1);   
            jet_file.read(reinterpret_cast<char*>(&y2)   , sizeof y2);   
            jet_file.read(reinterpret_cast<char*>(&phi)  , sizeof phi);  
            jet_file.read(reinterpret_cast<char*>(&tata) , sizeof tata); 
            jet_file.read(reinterpret_cast<char*>(&f1) , sizeof f1); 
            jet_file.read(reinterpret_cast<char*>(&f2) , sizeof f2); 

            jets.emplace_back(t01, t02, x, y, pt, y1, y2, phi, tata, f1, f2);

            if (f1 != 0)
            {
                quarks += static_cast<int>(f1 / abs(f1));
                if (f1 % 2 == 0)
                {
                    charges += 2*static_cast<int>(f1 / abs(f1));
                }
                else
                {
                    charges += -1*static_cast<int>(f1 / abs(f1));
                }
            }
            if (f2 != 0)
            {
                quarks += static_cast<int>(f2 / abs(f2));
                if (f2 % 2 == 0)
                {
                    charges += 2*static_cast<int>(f2 / abs(f2));
                }
                else
                {
                    charges += -1*static_cast<int>(f2 / abs(f2));
                }
            }
        }
        total_B += static_cast<double>(quarks)/3.0;
        total_Q += static_cast<double>(charges)/3.0;
        events_jets.emplace_back(std::move(jets));
    }
    jet_file.close();

    total_B /= static_cast<double>(events_jets.size());
    total_Q /= static_cast<double>(events_jets.size());

    std::cout << "Done! Read " << events_jets.size() << " events totaling "
              << n_dijet_total << " dijets, average baryon number was "
              << total_B << " and average em-charge was "
              << total_Q <<std::endl;

    double delta = (8.6 + 8.6)/n_bins;
    std::ofstream dBdy;
    std::ofstream dQdy;
    dBdy.open("dBdy.csv");
    dQdy.open("dQdy.csv");
    for (int i = 0; i < n_bins+1; i++)
    {
        dBdy << -8.6 + i*delta << ',';
        dQdy << -8.6 + i*delta << ',';
    }
    dBdy << std::endl;
    dQdy << std::endl;

    std::vector<int> B_in_bin(n_bins, 0);
    std::vector<int> Q_in_bin(n_bins, 0);

    double B_var = 0.0;
    double Q_var = 0.0;
    double B_av = 0.0;
    double Q_av = 0.0;

    for (auto event : events_jets)
    {
        double B_ev = 0.0;
        double Q_ev = 0.0;
        for (auto dijet : event)
        {
            B_ev += dijet.b1 + dijet.b2;
            Q_ev += dijet.q1 + dijet.q2;
            if (dijet.b1 != 0)
            {
                bool found =false;
                for (int i = 0; i < n_bins-1; i++)
                {
                    if ( -8.6+(i+1)*delta > dijet.y1)
                    {
                        found = true;
                        B_in_bin[i] += dijet.b1;
                        Q_in_bin[i] += dijet.q1;
                        break;
                    }
                }
                
                if (!found)
                {
                    B_in_bin[n_bins-1] += dijet.b1;
                    Q_in_bin[n_bins-1] += dijet.q1;
                }
            }
            if (dijet.b2 != 0)
            {
                bool found =false;
                for (int i = 0; i < n_bins-1; i++)
                {
                    if ( -8.6+(i+1)*delta > dijet.y2)
                    {
                        found = true;
                        B_in_bin[i] += dijet.b2;
                        Q_in_bin[i] += dijet.q2;
                        break;
                    }
                }
                
                if (!found)
                {
                    B_in_bin[n_bins-1] += dijet.b2;
                    Q_in_bin[n_bins-1] += dijet.q2;
                }
            }
        }
        B_var += pow(B_ev/3.0 - total_B, 2);
        Q_var += pow(Q_ev/3.0 - total_Q, 2);
        B_av += B_ev/3.0;
        Q_av += Q_ev/3.0;
    }

    std::cout<<"Variances were dB="<<B_var / static_cast<double>(events_jets.size())
             <<" and dQ="<<Q_var / static_cast<double>(events_jets.size())<<std::endl;
    std::cout<<"Averages again were dB="<<B_av / static_cast<double>(events_jets.size())
             <<" and dQ="<<Q_av / static_cast<double>(events_jets.size())<<std::endl;
        
    for (int i = 0; i < n_bins+1; i++)
    {
        dBdy << B_in_bin[i] / static_cast<double>(3.0*events_jets.size()*delta) << ',';
        dQdy << Q_in_bin[i] / static_cast<double>(3.0*events_jets.size()*delta) << ',';
    }
    dBdy << std::endl;
    dBdy.close();
    dQdy << std::endl;
    dQdy.close();
}


//0-100
//Trying to read 75000 events from the file jets_0_100_LHC_ref_flav.dat ...
//Done! Read 75000 events totaling 117326693 dijets, average baryon number was 48.5017 and average em-charge was 57.2184
//Variances were dB=12405.6 and dQ=17777
//     std::cout << "x = [";
//     double delta = (40.0-1.0)/n_bins;
//     for (int i = 0; i < n_bins; i++)
//     {
//         std::cout << 1.0+i*delta << ',';
//     }
//     std::cout << 1.0+n_bins*delta << ']' << std::endl;

//     std::vector<std::vector<double> > pts_in_bin
//     (
//         n_bins+1,
//         std::vector<double>()
//     );

//     for (auto event : events_jets)
//     {
//         for (auto dijet : event)
//         {
//             //if (3 < dijet.y1 && 4 > dijet.y1)
//             {
//             bool found =false;
//             for (int i = 0; i < n_bins; i++)
//             {
//                 if ( 1.0+(i+1)*delta > dijet.pt)
//                 {
//                     found = true;
//                     pts_in_bin[i].push_back(dijet.pt);
//                     break;
//                 }
//             }
//             if (!found)
//             {
//                 pts_in_bin[n_bins].push_back(dijet.pt);
//             }}
//             //if (3 < dijet.y2 && 4 > dijet.y2)
//             {
//             bool found =false;
//             for (int i = 0; i < n_bins; i++)
//             {
//                 if ( 1.0+(i+1)*delta > dijet.pt)
//                 {
//                     found = true;
//                     pts_in_bin[i].push_back(dijet.pt);
//                     break;
//                 }
//             }
//             if (!found)
//             {
//                 pts_in_bin[n_bins].push_back(dijet.pt);
//             }}
//         }
//     }

//     std::cout << "y = [";
//     for (int i = 0; i < n_bins; i++)
//     {
//         std::cout << pts_in_bin[i].size()/(events_jets.size()*delta) << ',';
//     }
//     std::cout << pts_in_bin[n_bins].size()/(events_jets.size()*delta) << ']' << std::endl;
// }