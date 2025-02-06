// Copyright (c) 2025 Mikko Kuha.
// This program is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software 
// Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU General Public License for more details. You should have received a copy
// of the GNU General Public License along with this program. If not, see
// <http://www.gnu.org/licenses/>.

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
     * @brief Adds a single point with a weight to the histogram
     * 
     * @param y_ The point and the weight
     */
    auto add
    (
        const std::tuple<double, double> &y_
    ) noexcept -> void;

    /**
     * @brief Adds multiple points to the histogram, faster than
     * adding one by one. Goes only once through each of the bins
     * 
     * @param news A vector of points to add to the histogram
     */
    auto add
    (
        const std::vector<double> &news
    ) noexcept -> void;

    /**
     * @brief Adds multiple points with varying weights to the 
     * histogram, faster than adding one by one. Goes only once 
     * through each of the bins
     * 
     * @param ys A vector of points and their weights to add to the 
     *           histogram
     */
    auto add
    (
        const std::vector<std::tuple<double, double> > &news
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
     *         <std::vector<double>,   = Normalized histogram
     *          double,                = Underflow
     *          double,                = Overflow
     *          std::vector<double>,   = Bin edge coords
     *          uint_fast64_t>              = Total # of points
     */
    auto get_histo() const noexcept -> std::tuple
    <   std::vector<double>, 
        double, 
        double, 
        std::vector<double>, 
        uint_fast64_t >;

    /**
     * @brief Checks whether the binning coords in the histograms
     * are the same. DOES NOT compare the points in the histograms
     * 
     * @param c1 One histogram
     * @param c2 Another histogram
     * @return true The bins are the same
     * @return false Not all the bins are same
     */
    friend bool operator==
    (
        const histo_1d& c1, 
        const histo_1d& c2
    ) noexcept;

protected:
private:
    const std::vector<double> xs;
    std::vector<double> counts;
    uint_fast64_t total_counts{0};
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
     * adding one by one. Goes only once through each of the bins
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
     *     uint_fast64_t >                           = Total # of points
     */
    auto get_histo() const noexcept -> std::tuple
    <   std::vector<std::vector<double> >, 
        std::tuple<double, double>, 
        std::tuple<double, double>, 
        std::tuple<std::vector<double>, std::vector<double> >,
        uint_fast64_t >;

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
     *          double,                = #of all underflow
     *          double,                = #of all overflow
     *          std::vector<double>&&, = Bin edge coords
     *          uint_fast64_t>              = Total # of points 
     */
    auto project_1d(const bool project_ys) const noexcept -> std::tuple
    <   std::vector<double>&&, 
        double, 
        double, 
        std::vector<double>&&, 
        uint_fast64_t >;

    /**
     * @brief Checks whether the binning coords in the histograms
     * are the same. DOES NOT compare the points in the histograms
     * 
     * @param c1 One histogram
     * @param c2 Another histogram
     * @return true The bins are the same
     * @return false Not all the bins are same
     */
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
    uint_fast64_t total_counts{0};
    std::tuple<double, double> underf{0.0, 0.0};
    std::tuple<double, double> overf{0.0, 0.0};
};

/**
 * @brief An object class to make 3D histograms
 * 
 */
class histo_3d
{
public:

    /**
     * @brief Construct a new 3D histogram object
     * 
     * @param xs_ Coordinates of the bin edges, 1st dim
     * @param ys_ Coordinates of the bin edges, 2nd dim
     * @param zs_ Coordinates of the bin edges, 3rd dim
     */
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

    /**
     * @brief Adds a single point to the histogram
     * 
     * @param x The coordinates of the point
     */
    auto add
    (
        const std::tuple<double, double, double> &x
    ) noexcept -> void;

    /**
     * @brief Adds multiple points to the histogram, faster than
     * adding one by one. Goes only once through each of the bins
     * 
     * @param news A vector of points to add to the histogram
     */
    auto add
    (
        std::vector<std::tuple<double, double, double> > news
    ) noexcept -> void;

    /**
     * @brief Adds all the points of another histogram to the histogram
     * 
     * @param other The other histogram from which the points are added
     */
    auto add
    (
        const histo_3d &other
    ) noexcept -> void;

    /**
     * @brief  Gets the collected data from the histogram. The returned
     * histogram values are normalized (divided) with the #of total
     * points and bin size (for each bin separately). The overflow and
     * the underflow are normalized only with #of total points.
     * 
     * @return std::tuple
     * <   std::vector<std::vector<std::vector<double>>>&&, = Normalized histogram
     *     std::tuple<double, double, double>,              = Underflow (all dims)
     *     std::tuple<double, double, double>,              = Overflow (all dims)
     *     std::tuple<std::vector<double>, 
     *                std::vector<double>, 
     *                std::vector<double> >,                = Bin edge coords (all dims)
     *     uint_fast64_t >                                       = Total # of points
     */
    auto get_histo() const noexcept -> std::tuple
    <   std::vector<std::vector<std::vector<double> > >, 
        std::tuple<double, double, double>, 
        std::tuple<double, double, double>, 
        std::tuple<std::vector<double>, std::vector<double>, std::vector<double> >,
        uint_fast64_t >;

    /**
     * @brief Integrate over two of the dimensions. The returned
     * histogram values are normalized (divided) with the #of total
     * points and bin size (for each bin separately). The overflow and
     * the underflow are normalized only with #of total points.
     * 
     * @param dim_left The dimension that is not integrated over. 0, 1 or 2.
     * 
     * @return std::tuple
     *         <std::vector<double>&&, = Normalized histogram
     *          double,                = #of all underflow
     *          double,                = #of all overflow
     *          std::vector<double>&&, = Bin edge coords
     *          uint_fast64_t>              = Total # of points 
     */
    auto project_1d(const uint_fast8_t dim_left) const noexcept -> std::tuple
    <   std::vector<double>&&, 
        double, 
        double, 
        std::vector<double>&&, 
        uint_fast64_t >;

    /**
     * @brief Integrate over one of the dimensions. The returned
     * histogram values are normalized (divided) with the #of total
     * points and bin size (for each bin separately). The overflow and
     * the underflow are normalized only with #of total points.
     * 
     * @param dim_integrated The integrated dimension. 0, 1 or 2.
     * 
     * @return std::tuple
     * <   std::vector<std::vector<double> >&&, = Normalized histogram
     *     double,                              = #of all underflow
     *     double,                              = #of all overflow
     *     std::tuple<std::vector<double>,      
     *                std::vector<double> >,    = Bin edge coords (both dims)
     *     uint_fast64_t>                            = Total # of points 
     */
    auto project_2d(const uint_fast8_t dim_integrated) const noexcept 
    -> std::tuple
    <   std::vector<std::vector<double> >, 
        double, 
        double, 
        std::tuple<std::vector<double>, std::vector<double> >,
        uint_fast64_t >;

    /**
     * @brief Checks whether the binning coords in the histograms
     * are the same. DOES NOT compare the points in the histograms
     * 
     * @param c1 One histogram
     * @param c2 Another histogram
     * @return true The bins are the same
     * @return false Not all the bins are same
     */
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
    uint_fast64_t total_counts{0};
    std::tuple<double, double, double> underf{0.0, 0.0, 0.0};
    std::tuple<double, double, double> overf{0.0, 0.0, 0.0};
};

#endif // HISTO_H
