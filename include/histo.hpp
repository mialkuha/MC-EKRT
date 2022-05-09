//Copyright (c) 2022 Mikko Kuha

#ifndef HISTO_H
#define HISTO_H

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <tuple>
#include <vector>

/**
 * @brief An object class to make 1D histograms
 * 
 */
class histo_1d
{
public:

    /**
     * @brief Construct a new 1D histogram object
     * 
     * @param xs_ Coordinates of the bin edges
     */
    histo_1d
    (
        const std::vector<double> &xs_
    ) noexcept : 
        xs(xs_),
        counts(xs_.size()-1, 0.0),
        total_counts(0),
        underf(0),
        overf(0)
    {}

    /**
     * @brief Construct a new 1D histogram object
     * 
     * @param xs_ Coordinates of the bin edges
     * @param ys_ A vector of points to add to the histogram
     */
    histo_1d
    (
        const std::vector<double> &xs_,
        const std::vector<double> &ys_
    ) noexcept : 
        xs(xs_),
        counts(xs_.size()-1, 0.0),
        total_counts(0),
        underf(0),
        overf(0)
    { this->add(ys_); }

    /**
     * @brief Adds a single point to the histogram
     * 
     * @param y The coordinate of the point
     */
    auto add
    (
        const double &y
    ) noexcept -> void;

    /**
     * @brief Adds multiple points to the histogram, faster than
     * adding one by one.
     * 
     * @param ys A vector of points to add to the histogram
     */
    auto add
    (
        const std::vector<double> &ys
    ) noexcept -> void;

    /**
     * @brief Adds all the points of another histogram to the histogram
     * 
     * @param other The other histogram from which the points are added
     */
    auto add
    (
        const histo_1d &other
    ) noexcept -> void;

    /**
     * @brief Gets the collected data from the histogram. The returned
     * histogram values are normalized (divided) with the #of total
     * points and bin size (for each bin separately). The overflow and
     * the underflow are normalized only with #of total points.
     * 
     * @return std::tuple
     *         <std::vector<double>&&, = Normalized histogram
     *          double,                = Underflow
     *          double,                = Overflow
     *          std::vector<double>,   = Bin edge coords
     *          uint64_t>              = Total # of points
     */
    auto get_histo() const noexcept -> std::tuple
    <   std::vector<double>&&, 
        double, 
        double, 
        std::vector<double>, 
        uint64_t >;

    friend bool operator==
    (
        const histo_1d& c1, 
        const histo_1d& c2
    ) noexcept;

protected:
private:
    const std::vector<double> xs;
    std::vector<double> counts;
    uint64_t total_counts{0};
    double underf{0.0};
    double overf{0.0};
};

/**
 * @brief An object class to make 2D histograms
 * 
 */
class histo_2d
{
public:

    /**
     * @brief Construct a new 2D histogram object
     * 
     * @param xs_ Coordinates of the bin edges, 1st dim
     * @param ys_ Coordinates of the bin edges, 2nd dim
     */
    histo_2d
    (
        const std::vector<double> &xs_,
        const std::vector<double> &ys_
    ) noexcept : 
        xs(xs_),
        ys(ys_),
        counts
        (
            xs_.size()-1, 
            std::vector<double>
            (
                ys_.size()-1, 
                0.0
            )
        ),
        total_counts(0),
        underf({0.0, 0.0}),
        overf({0.0, 0.0})
    {}

    /**
     * @brief Adds a single point to the histogram
     * 
     * @param x The coordinates of the point
     */
    auto add
    (
        const std::tuple<double, double> &x
    ) noexcept -> void;

    /**
     * @brief Adds multiple points to the histogram, faster than
     * adding one by one.
     * 
     * @param news A vector of points to add to the histogram
     */
    auto add
    (
        std::vector<std::tuple<double, double> > news
    ) noexcept -> void;

    /**
     * @brief Adds all the points of another histogram to the histogram
     * 
     * @param other The other histogram from which the points are added
     */
    auto add
    (
        const histo_2d &other
    ) noexcept -> void;

    /**
     * @brief  Gets the collected data from the histogram. The returned
     * histogram values are normalized (divided) with the #of total
     * points and bin size (for each bin separately). The overflow and
     * the underflow are normalized only with #of total points.
     * 
     * @return std::tuple
     * <   std::vector<std::vector<double> >&&, = Normalized histogram
     *     std::tuple<double, double>,          = Underflow (both dims)
     *     std::tuple<double, double>,          = Overflow (both dims)
     *     std::tuple<std::vector<double>,      
     *                std::vector<double> >,    = Bin edge coords (both dims)
     *     uint64_t >                           = Total # of points
     */
    auto get_histo() const noexcept -> std::tuple
    <   std::vector<std::vector<double> >&&, 
        std::tuple<double, double>, 
        std::tuple<double, double>, 
        std::tuple<std::vector<double>, std::vector<double> >,
        uint64_t >;

    /**
     * @brief Integrate over one of the dimensions. The returned
     * histogram values are normalized (divided) with the #of total
     * points and bin size (for each bin separately). The overflow and
     * the underflow are normalized only with #of total points.
     * 
     * @param project_ys True  = integrate over 2nd dimension
     *                   False = integrate over 1st dimension
     * 
     * @return std::tuple
     *         <std::vector<double>&&, = Normalized histogram
     *          double,                = Underflow
     *          double,                = Overflow
     *          std::vector<double>&&, = Bin edge coords
     *          uint64_t>              = Total # of points 
     */
    auto project_1d(const bool project_ys) const noexcept -> std::tuple
    <   std::vector<double>&&, 
        double, 
        double, 
        std::vector<double>&&, 
        uint64_t >;

    friend bool operator==
    (
        const histo_2d& c1, 
        const histo_2d& c2
    ) noexcept;

protected:
private:
    const std::vector<double> xs;
    const std::vector<double> ys;
    std::vector<std::vector<double> > counts;
    uint64_t total_counts{0};
    std::tuple<double, double> underf{0.0, 0.0};
    std::tuple<double, double> overf{0.0, 0.0};
};

class histo_3d
{
public:
    histo_3d
    (
        const std::vector<double> &xs_,
        const std::vector<double> &ys_,
        const std::vector<double> &zs_
    ) noexcept : 
        xs(xs_),
        ys(ys_),
        zs(zs_),
        counts
        (
            xs_.size()-1, 
            std::vector<std::vector<double> >
            (
                ys_.size()-1, 
                std::vector<double>
                (
                    zs_.size()-1, 
                    0.0
                )
            )
        ),
        total_counts(0),
        underf({0.0, 0.0, 0.0}),
        overf({0.0, 0.0, 0.0})
    {}

    auto add
    (
        const std::tuple<double, double, double> &x
    ) noexcept -> void;

    auto add
    (
        std::vector<std::tuple<double, double, double> > news
    ) noexcept -> void;

    auto add
    (
        const histo_3d &other
    ) noexcept -> void;

    auto get_histo() const noexcept;

    auto project_1d(const uint8_t dim_left) const noexcept;

    friend bool operator==
    (
        const histo_3d& c1, 
        const histo_3d& c2
    ) noexcept;

protected:
private:
    const std::vector<double> xs;
    const std::vector<double> ys;
    const std::vector<double> zs;
    std::vector<std::vector<std::vector<double> > > counts;
    uint64_t total_counts{0};
    std::tuple<double, double, double> underf{0.0, 0.0, 0.0};
    std::tuple<double, double, double> overf{0.0, 0.0, 0.0};
};

#endif // HISTO_H
