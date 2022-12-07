//Copyright (c) 2022 Mikko Kuha 
//TODO implement NN-tree 

#include "mcaa.hpp"

int main()
{
    std::vector<double> p0s{/*0.75,*/1/*,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4*/};
    std::vector<double> Ks{0.5,1,1.5,2,2.5,3,3.5,4};
    std::vector<double> Ms{1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6};

    std::ofstream file;
    
    std::string dummy{"for_sigma_jets"};
    file.open(dummy);
    file << "name " <<dummy<<std::endl;
    file << "kt0 " <<1<<std::endl;
    file << "K_factor " <<1<<std::endl;
    file << "calculate_tata false"<<std::endl;
    file << "use_snpdfs true"<<std::endl;
    file << "read_sigmajets_from_file false"<<std::endl;
    file << "n_events 10"<<std::endl;
    file.close();

    mcaa sim0(dummy);
    sim0.run();

    for (auto p0 : p0s)
    {
        for (auto K : Ks)
        {
            for (auto M : Ms)
            {
                std::stringstream name{""};
                name << p0<<'_'<<K<<'_'<<M<<"_s";
                file.open(name.str());
                file << "name " <<name.str()<<std::endl;
                file << "kt0 " <<p0<<std::endl;
                file << "K_factor " <<K<<std::endl;
                file << "M_factor " <<M<<std::endl;
                file << "calculate_tata false"<<std::endl;
                file << "use_snpdfs true"<<std::endl;
                file << "read_sigmajets_from_file true"<<std::endl;
                file.close();

                mcaa sim(name.str());
                sim.run();
            }
        }
    }

    return 0;
}