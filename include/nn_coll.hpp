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
    nn_coll(nucleon *const A, nucleon *const B, const momentum &sqrt_s_) noexcept 
        : dijets(0), target(B), projectile(A), sqrt_s(sqrt_s_), max_inel_xsect(0), effective_inel_xsect(0), max_tot_xsect(0), effective_tot_xsect(0)
    {
        this->bsquared = this->projectile->calculate_bsquared(*(this->target));
    }
    //nn_coll& operator=(const nn_coll &) = default;
    //nn_coll& operator=(nn_coll&&) = default;
    //nn_coll(const nn_coll &) = default;
    //nn_coll(nn_coll&& other) 
    //    : dijets(std::move(other.dijets)), sqrt_s(other.sqrt_s), bsquared(other.bsquared), max_inel_xsect(other.max_inel_xsect),
    //      effective_inel_xsect(other.effective_inel_xsect), max_tot_xsect(other.max_tot_xsect), effective_tot_xsect(other.effective_tot_xsect) 
    //{
    //    this->target = other.target;
    //    this->projectile = other.projectile;
    //    other.target = nullptr;
    //    other.projectile = nullptr;
    //}

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
    void reduce_energy_and_push_end_states_to_collider_frame() noexcept;
    spatial randomize_bsquared_for_this(const spatial & min, const spatial & max, std::mt19937 random_generator_) noexcept;
        //bgen = TMath::Sqrt((fBmax*fBmax-fBmin*fBmin)*gRandom->Rndm()+fBmin*fBmin);

    static spatial randomize_bsquared(const spatial & min, const spatial & max, std::mt19937 random_generator_) noexcept;

    std::vector<dijet_specs> dijets;
    nucleon *target;
    nucleon *projectile;
    
protected:
private:

    momentum sqrt_s{0};
    spatial bsquared{0};
    xsectval max_inel_xsect{0};       //Inelastic xsect at b=0
    xsectval effective_inel_xsect{0}; //Inelastic xsect at this b
    xsectval max_tot_xsect{0};       //Total xsect at b=0
    xsectval effective_tot_xsect{0}; //Total xsect at this b
};

#endif // NN_COLL_HPP
