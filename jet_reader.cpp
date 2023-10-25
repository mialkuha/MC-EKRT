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
    uint_fast16_t init1;
    uint_fast16_t init2;
    uint_fast16_t final1;
    uint_fast16_t final2;
    uint_fast16_t ia;
    uint_fast16_t ib;
    dijet(const double &t01_, const double &t02_, const double &x_, const double &y_, 
            const double &pt_, const double &y1_, const double &y2_, const double &phi_, 
            const double &tata_, const uint_fast16_t &init1_, const uint_fast16_t &init2_,
            const uint_fast16_t &final1_, const uint_fast16_t &final2_, const uint_fast16_t &ia_, 
            const uint_fast16_t &ib_)
    : t01(t01_), t02(t02_), x(x_), y(y_), pt(pt_), y1(y1_), y2(y2_), phi(phi_) , tata(tata_),
      init1(init1_), init2(init2_), final1(final1_), final2(final2_), ia(ia_), ib(ib_) { }
};

 
int main()
{
    std::string filename = "jets_0_5_example_name.dat";
    std::ifstream jet_file(filename, std::ios::in | std::ios::binary);

    if (!jet_file.is_open())
    {
        std::cout << "Error opening file " << filename << std::endl;
        exit(1);
    }

    uint_fast64_t n_events, n_jets;
    double t01, t02, x, y, pt, y1, y2, phi, tata;
    uint_fast16_t init1, init2, final1, final2, ia, ib;
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

    uint_fast16_t n_bins = 20;
    double largest_tata = 0;

    uint_fast64_t n_dijet_total = 0;
    for (uint_fast64_t i=0; i<n_events; i++)
    {
        jet_file.read(reinterpret_cast<char*>(&n_jets), sizeof n_jets);

        std::vector<dijet> jets;
        jets.reserve(n_jets);

        for (uint_fast64_t ii=0; ii<n_jets; ii++, n_dijet_total++)
        {
            jet_file.read(reinterpret_cast<char*>(&t01)    , sizeof t01);    //t01
            jet_file.read(reinterpret_cast<char*>(&t02)    , sizeof t02);    //t02
            jet_file.read(reinterpret_cast<char*>(&x)      , sizeof x);      //x
            jet_file.read(reinterpret_cast<char*>(&y)      , sizeof y);      //y
            jet_file.read(reinterpret_cast<char*>(&pt)     , sizeof pt);     //pT
            jet_file.read(reinterpret_cast<char*>(&y1)     , sizeof y1);     //y1
            jet_file.read(reinterpret_cast<char*>(&y2)     , sizeof y2);     //y2
            jet_file.read(reinterpret_cast<char*>(&phi)    , sizeof phi);    //phi_1
            jet_file.read(reinterpret_cast<char*>(&tata)   , sizeof tata);   //T_A * T_A
            jet_file.read(reinterpret_cast<char*>(&init1)  , sizeof init1);  //flavour of incoming 1
            jet_file.read(reinterpret_cast<char*>(&init2)  , sizeof init2);  //flavour of incoming 2
            jet_file.read(reinterpret_cast<char*>(&final1) , sizeof final1); //flavour of outgoing 1
            jet_file.read(reinterpret_cast<char*>(&final2) , sizeof final2); //flavour of outgoing 2
            jet_file.read(reinterpret_cast<char*>(&ia)     , sizeof ia);     //index of mother a
            jet_file.read(reinterpret_cast<char*>(&ib)     , sizeof ib);     //index of mother b

            jets.emplace_back(t01, t02, x, y, pt, y1, y2, phi, tata, init1, init2, final1, final2, ia, ib);
			//std::cout<<t01<<' '<<t02<<' '<<x<<' '<<y<<' '<<pt<<' '<<y1<<' '<<y2<<' '<<phi<<' '<<tata<<std::endl;
            if (largest_tata < tata)
            {
                largest_tata = tata;
            }
        }
        events_jets.emplace_back(std::move(jets));
    }
    jet_file.close();

    std::cout << "Done! Read " << events_jets.size() << " events totaling "
              << n_dijet_total << " dijets, the largest tata was "
              << largest_tata <<std::endl;

    double delta = largest_tata*1.01/n_bins;
    std::ofstream smallest_pts;
    smallest_pts.open("smallest_pts_0_5.csv");
    for (uint i = 0; i < n_bins; i++)
    {
        smallest_pts << (i+1)*delta << ',';
    }
    smallest_pts << std::endl;

    std::vector<std::vector<double> > pts_in_event
    (
        n_bins, 
        std::vector<double>()
    );

    for (auto event : events_jets)
    {
        for (auto dijet : event)
        {
            //std::cout << dijet.tata << ' '; 
            for (uint i = 0; i < n_bins; i++)
            {
                if ( (i+1)*delta > dijet.tata)
                {
                    //std::cout << (i+1)*delta << ' ' << ;
                    pts_in_event[i].push_back(dijet.pt);
                    break;
                }
            }
        }
        
        for (uint i = 0; i < n_bins; i++)
        {
            if (pts_in_event[i].size() != 0)
            {
                std::sort(pts_in_event[i].begin(),pts_in_event[i].end());
                smallest_pts << pts_in_event[i][0] << ',';
                pts_in_event[i].clear();
            }
            else
            {
                smallest_pts << "-1" << ',';
            }
        }
        smallest_pts << std::endl;
    }
    smallest_pts.close();
}

