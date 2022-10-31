#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

struct minijet
{
    double t0;
    double x;
    double y;
    double z;
    double energy;
    double px;
    double py;
    double pz;
    double tata;
    minijet(const double &t0_, const double &x_, const double &y_, const double &z_, 
            const double &energy_, const double &px_, const double &py_, const double &pz_, 
            const double &tata_)
    : t0(t0_), x(x_), y(y_), z(z_), energy(energy_), px(px_), py(py_), pz(pz_) , tata(tata_) { }
};

 
int main()
{
    std::string filename = "jets_0_5.dat";
    std::ifstream jet_file(filename, std::ios::in | std::ios::binary);

    if (!jet_file.is_open())
    {
        std::cout << "Error opening file " << filename << std::endl;
        exit(1);
    }

    uint64_t n_events, n_jets;
    double t0, x, y, z, energy, px, py, pz, tata;
    std::vector<std::vector<minijet> > events_jets;
    
    if (sizeof n_jets != 8 || sizeof t0 != 8 )
    {
        std::cout << "Change the types! n_jets should be 64 bit unsigned integer"
                  << " and the coordinates should all be 64 bit float numbers" << std::endl;
        exit(1);
    }

    jet_file.read(reinterpret_cast<char*>(&n_events), sizeof n_events);
    std::cout << "Trying to read " << n_events << " events from the file " 
              << filename << " ..." << std::endl;
    events_jets.reserve(n_events);

    uint64_t n_jet_total = 0;
    for (uint64_t i=0; i<n_events; i++)
    {
        jet_file.read(reinterpret_cast<char*>(&n_jets), sizeof n_jets);

        events_jets[i].reserve(n_jets);

        for (uint64_t i=0; i<n_jets; i++, n_jet_total++)
        {
            jet_file.read(reinterpret_cast<char*>(&t0)     , sizeof t0);     //t0
            jet_file.read(reinterpret_cast<char*>(&x)      , sizeof x);      //x
            jet_file.read(reinterpret_cast<char*>(&y)      , sizeof y);      //y
            jet_file.read(reinterpret_cast<char*>(&z)      , sizeof z);      //z
            jet_file.read(reinterpret_cast<char*>(&energy) , sizeof energy); //E
            jet_file.read(reinterpret_cast<char*>(&px)     , sizeof px);     //p_x
            jet_file.read(reinterpret_cast<char*>(&py)     , sizeof py);     //p_y
            jet_file.read(reinterpret_cast<char*>(&pz)     , sizeof pz);     //p_z
            jet_file.read(reinterpret_cast<char*>(&tata)   , sizeof tata);   //T_A * T_A

            events_jets[i].emplace_back(t0, x, y, z, energy, px, py, pz, tata);
        }
        jet_file.close();
    }

    std::cout << "Done! Read " << events_jets.size() << " events totaling "
              << n_jet_total << " jets, the last one had "
              << energy << " GeV energy"<<std::endl;

    double pt_cumulant = 0;
    for (auto event : events_jets)
    {
        for (auto jet : event)
        {
            pt_cumulant += std::sqrt(jet.px*jet.px + jet.py*jet.py);
        }
    }
    std::cout << "In all of the events, average jet had pT=" 
              << pt_cumulant / n_jet_total << " GeV"<<std::endl;
}

