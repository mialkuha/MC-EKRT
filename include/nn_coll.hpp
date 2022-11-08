//Copyright (c) 2022 Mikko Kuha

#ifndef NN_COLL_HPP
#define NN_COLL_HPP

#include <cmath>
#include <functional>
#include <random>

#include "nucleon.hpp"
#include "typedefs.hpp"

class nn_coll
{
public:
    explicit nn_coll(const coords &) noexcept;
    nn_coll(nucleon *const A, nucleon *const B, const double &sqrt_s_) noexcept 
        : dijets(0), target(B), projectile(A), sqrt_s(sqrt_s_), max_inel_xsect(0), effective_inel_xsect(0), max_tot_xsect(0), effective_tot_xsect(0)
    {
        this->bsquared = this->projectile->calculate_bsquared(*(this->target));
    }

    const double & getcr_sqrt_s() const noexcept {return this->sqrt_s;} 
    const double & getcr_bsquared() const noexcept {return this->bsquared;}
    const double & getcr_max_inel_xsect() const noexcept {return this->max_inel_xsect;}
    const double & getcr_effective_inel_xsect() const noexcept {return this->effective_inel_xsect;}
    const double & getcr_max_tot_xsect() const noexcept {return this->max_tot_xsect;}
    const double & getcr_effective_tot_xsect() const noexcept {return this->effective_tot_xsect;}
    
    void reduce_energy() noexcept;
    void wound() noexcept;
    void calculate_xsects(const double &sigma_jet, const std::function<double(const double&)> &Tpp, const double &bsquared, const B2_normalization_mode & mode)  noexcept;
    void push_end_states_to_collider_frame() noexcept;
    void reduce_energy_and_push_end_states_to_collider_frame() noexcept;
    double randomize_bsquared_for_this(const double & min, const double & max, std::mt19937 random_generator_) noexcept;
        //bgen = TMath::Sqrt((fBmax*fBmax-fBmin*fBmin)*gRandom->Rndm()+fBmin*fBmin);

    static double randomize_bsquared(const double & min, const double & max, std::mt19937 random_generator_) noexcept;

    std::vector<dijet_specs> dijets;
    nucleon *target;
    nucleon *projectile;
    
protected:
private:

    double sqrt_s{0};
    double bsquared{0};
    double max_inel_xsect{0};       //Inelastic xsect at b=0
    double effective_inel_xsect{0}; //Inelastic xsect at this b
    double max_tot_xsect{0};       //Total xsect at b=0
    double effective_tot_xsect{0}; //Total xsect at this b
};

#endif // NN_COLL_HPP
