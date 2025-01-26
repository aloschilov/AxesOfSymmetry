#include <iostream>
#include <numeric>
#include <map>
#include "matrix.h"

enum RelationToOx
{
    ABOVE,
    BELOW,
    ON
};

//////////////////////////////////////////////////////////////////////////////
/// @brief Splits points of those above ox and below given
///
/// @param relations - input - precalculated vector of relations {ABOVE, ON, BELOW, BELOW, ON, ABOVE}
///
/// @returns - pairs of indexes such as {4, 5, 0}, {1, 2, 3}
///
//////////////////////////////////////////////////////////////////////////////
std::pair<std::vector<unsigned int>, std::vector<unsigned int>>
split_above_to_below(const std::vector<RelationToOx> &relations);

//////////////////////////////////////////////////////////////////////////////
/// @brief Calculates the (d, alpha) representation of original polygon
///
/// @param points - input - vector of regular points represented as column vector
///
/// @returns - vector of pairs (d, alpha), where d - edge length, and alpha - angle of the edge with ox
///
//////////////////////////////////////////////////////////////////////////////
std::vector<std::pair<double, double>> calculate_representation(const std::vector<Matrix> &points);

//////////////////////////////////////////////////////////////////////////////
/// @brief Rotates representation of (d, alpha) on provided angle
///
/// @param representation - input - vector of pairs (d, alpha), where d - edge length, and alpha - angle of the edge with ox
/// @param angle - input - angle value in radian
///
/// @returns - vector of pairs (d, alpha), where d - edge length, and alpha - angle of edge with ox
///
//////////////////////////////////////////////////////////////////////////////
std::vector<std::pair<double, double>> rotate_representation(const std::vector<std::pair<double, double>> &representation, double angle);


//////////////////////////////////////////////////////////////////////////////
/// @brief Checks whether there is a symmetry based on representation and index of point above and below
///
/// @param above - input - vector of indexes considered above ox. See split_above_to_below
/// @param below - input - vector of indexes considered below ox. See split_above_to_below
/// @param representation - input - vector of pairs (d, alpha), where d - edge length, and alpha - angle of the edge with ox
///
/// @returns - vector of pairs (d, alpha), where d - edge length, and alpha - angle of edge with ox
///
//////////////////////////////////////////////////////////////////////////////
bool check_symmetry(
        const std::vector<unsigned int> &above,
        const std::vector<unsigned int> &below,
        const std::vector<std::pair<double, double>> &representation);

//////////////////////////////////////////////////////////////////////////////
/// @brief Gets polygon symmetry lines by polygon vertices
///
/// @param points - input - vector of regular points represented as column vector
///
/// @returns - vector of pairs with regular points represented as column vector
///
//////////////////////////////////////////////////////////////////////////////
std::vector<std::pair<size_t, double>> get_symmetry_angles(const std::vector<Matrix> &points);
//////////////////////////////////////////////////////////////////////////////
/// @brief Returns mapping of points to its angles
///
/// @param points - input - vector of regular points represented as column vector
///
/// @returns - vector of angles in radians
///
//////////////////////////////////////////////////////////////////////////////
std::vector<double> get_point_angles(const std::vector<Matrix> &points);

//////////////////////////////////////////////////////////////////////////////
/// @brief Filters symmetry angles candidates
///
/// @param point_angles - input - angle which point vector makes with ox
/// @param symmetry_angles - input - vector of pairs (idx, angle), where idx is index of point taken
///                                   and angle is a associated angle
/// @param representation - input - vector of pairs (d, alpha), where d - edge length, and alpha - angle of the edge with ox
///
/// @returns - symmetry_angles filtered by symmetry check
///
//////////////////////////////////////////////////////////////////////////////
std::vector<std::pair<size_t, double>> filter_symmetry_angles(const std::vector<double> &point_angles,
                                                              const std::vector<std::pair<size_t, double>> &symmetry_angles,
                                           const std::vector<std::pair<double, double>> &representation);
//////////////////////////////////////////////////////////////////////////////
/// @brief Add mid edges points to original representation
///
/// @param points - input - original points as vector of column vectors
///
/// @returns - points with added mid-edge points
///
//////////////////////////////////////////////////////////////////////////////
std::vector<Matrix> get_representation_with_mid_points(const std::vector<Matrix> &points);

//////////////////////////////////////////////////////////////////////////////
/// @brief Returns symmetry lines of a polygon
///
/// @param original_points - input - original points as vector of column vectors
///
/// @returns - vector of pairs with two points, which represent symmetry line
///
//////////////////////////////////////////////////////////////////////////////
std::vector<std::pair<Matrix, Matrix>> get_symmetry_lines(const std::vector<Matrix> &original_points);
