//Copyright (c) 2021 Mikko Kuha

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
    nn_coll(nucleon * A, nucleon * B, const momentum & sqrt_s_) noexcept;

    const momentum & getcr_sqrt_s() const noexcept {return this->sqrt_s;} 
    const spatial & getcr_bsquared() const noexcept {return this->bsquared;}
    const xsectval & getcr_max_inel_xsect() const noexcept {return this->max_inel_xsect;}
    const xsectval & getcr_effective_inel_xsect() const noexcept {return this->effective_inel_xsect;}
    const xsectval & getcr_max_tot_xsect() const noexcept {return this->max_tot_xsect;}
    const xsectval & getcr_effective_tot_xsect() const noexcept {return this->effective_tot_xsect;}
    
    void reduce_energy() noexcept;
    void wound() noexcept;
    void calculate_xsects(const xsectval &sigma_jet, const std::function<spatial(const spatial&)> &Tpp, const spatial &bsquared, const B2_normalization_mode & mode)  noexcept;
    void push_end_states_to_collider_frame() noexcept;
    spatial randomize_bsquared_for_this(const spatial & min, const spatial & max, std::mt19937 random_generator_) noexcept;
        //bgen = TMath::Sqrt((fBmax*fBmax-fBmin*fBmin)*gRandom->Rndm()+fBmin*fBmin);

    static spatial randomize_bsquared(const spatial & min, const spatial & max, std::mt19937 random_generator_) noexcept;

    std::vector<dijet_specs> dijets;
protected:
private:

    nucleon *target;
    nucleon *projectile;
    momentum sqrt_s{0};
    spatial bsquared{0};
    xsectval max_inel_xsect{0};       //Inelastic xsect at b=0
    xsectval effective_inel_xsect{0}; //Inelastic xsect at this b
    xsectval max_tot_xsect{0};       //Total xsect at b=0
    xsectval effective_tot_xsect{0}; //Total xsect at this b
};

#endif // NN_COLL_HPP
