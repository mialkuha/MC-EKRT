//Copyright (c) 2021 Mikko Kuha

#include "nucleus_generator.hpp"

std::uniform_real_distribution<double> nucleus_generator::unif_dist = std::uniform_real_distribution<double>(0.0,1.0);

std::vector<nucleon> nucleus_generator::generate_nucleus(const nucleus_params & params, const bool &target, const momentum & mom, 
        const spatial & xshift, std::shared_ptr<std::mt19937> random_generator, std::shared_ptr<ars> radial_sampler)
{
    std::vector<nucleon> generated_nucleus;
    std::vector<coords> generated_coords;
    auto N = (target) ? params.NB : params.NA;
    generated_nucleus.reserve(N);
    generated_coords.reserve(N);

    if (N==1)
    {
        generated_nucleus.emplace_back(coords({xshift,0.0,0.0}), mom);
        return generated_nucleus;
    }

    coords new_coords;
    bool coords_do_fit;
    auto count=0;

    do // while (generated_coords.size() < N);
    {
        do // while (!coords_do_fit);
        {
            try
            {
                new_coords = nucleus_generator::throw_nucleon_coords(random_generator, radial_sampler);
                coords_do_fit = nucleus_generator::coords_fit(new_coords, generated_coords, params.min_distance);
                //std::cout<<"coords "<<new_coords.x<<' '<<new_coords.y<<' '<<new_coords.z<<" fit? "<<coords_do_fit<<std::endl;
                if (params.correct_overlap_bias && !coords_do_fit)
                {
                    break;
                }
            }
            catch(const std::exception& e)
            {
                std::cout << e.what() << " in nucleus_generator, trying again"<<std::endl;
                coords_do_fit = false;
            }
        } while (!coords_do_fit);

        if (params.correct_overlap_bias && !coords_do_fit)
        {
            generated_coords.clear();
            count++;
            //std::cout<<"throw! "<<count<<' '<<std::flush;
        }
        generated_coords.push_back(std::move(new_coords));

    } while (generated_coords.size() < N);
    
    //std::cout<<"threw away "<<count<<" times\n";

    coords com{0.0,0.0,0.0}; //center of mass

    if (params.shift_cms)
    {
        for (const auto& co : generated_coords)
        {
            com += co;
        }
        com /= static_cast<double>(generated_coords.size());
    }

    for (auto& co : generated_coords)
    {
        if (params.shift_cms)
        {
            co -= com; //Shift all nucleons so that CoM is (0,0,0)
        }
        co.x += xshift;
        generated_nucleus.emplace_back(co, mom);
    }

    nucleus_generator::throw_neutrons(&generated_nucleus, (target) ? params.ZB : params.ZA, random_generator);

    return generated_nucleus;
}

coords nucleus_generator::throw_nucleon_coords(std::shared_ptr<std::mt19937> random_generator, std::shared_ptr<ars> radial_sampler) noexcept
{
    double new_r=-1, new_phi, new_cos_theta, new_sin_theta, rand1,rand2;

    do
    {
        try
        {
            new_r = radial_sampler->throw_one(*random_generator);
        }
        catch(const std::exception& e)
        {
            std::cout<<std::endl<<"Threw: "<<e.what()<<std::endl;
            new_r=-1;
        }
    } while (new_r < 0);
    
    rand1 = nucleus_generator::unif_dist(*random_generator);
    rand2 = nucleus_generator::unif_dist(*random_generator);
    new_phi = 2*M_PI*rand1;
    new_cos_theta = 2*rand2-1.;
    new_sin_theta = sqrt(1 - pow(new_cos_theta,2));

    //coords ret = coords({new_r * new_sin_theta * cos(new_phi),
    //                     new_r * new_sin_theta * sin(new_phi),
    //                     new_r * new_cos_theta});

    //std::cout<<"threw: "<<ret.x<<' '<<ret.y<<' '<<ret.z<<' '<<std::endl;
    //std::cout<<"threw: "<<new_r<<' '<<rand1<<' '<<rand2<<' '<<std::endl;

    //return std::move(ret);
    return coords({new_r * new_sin_theta * cos(new_phi),
                   new_r * new_sin_theta * sin(new_phi),
                   new_r * new_cos_theta});
}

void nucleus_generator::throw_neutrons(std::vector<nucleon> *const nucleus, const uint_fast16_t & Z, std::shared_ptr<std::mt19937> random_generator) noexcept
{
    for (uint_fast16_t i=0,iz=0; i<nucleus->size(); ++i)
    {
        double frac = static_cast<double>(Z-iz) / static_cast<double>(nucleus->size()-i);
        double rn = nucleus_generator::unif_dist(*random_generator);
        if (rn<frac) 
        {
            nucleus->at(i).is_neutron = false;
            ++iz;
        }
        else
        {
            nucleus->at(i).is_neutron = true;
        }
    }

}

bool nucleus_generator::coords_fit(const coords& co, const std::vector<coords>& other_coords, const spatial& min_distance) noexcept
{
    if (min_distance<=0)
    {
        return true;
    }
    
    const spatial md2 = pow(min_distance,2);
    
    for (const auto& oco : other_coords)
    {
        if ((co - oco).mag2() < md2) 
        {
            //std::cout<<co.x<<' '<<co.y<<' '<<co.z<<" is too near to "<<oco.x<<' '<<oco.y<<' '<<oco.z<<'\n';
            return false;
        }
    }
    return true;
}
