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

std::pair<std::vector<unsigned int>, std::vector<unsigned int>>
split_above_to_below(const std::vector<RelationToOx> &relations);
std::vector<std::pair<double, double>> calculate_representation(const std::vector<Matrix> &points);
std::vector<std::pair<double, double>> rotate_representation(const std::vector<std::pair<double, double>> &representation, double angle);
bool check_symmetry(
        const std::vector<unsigned int> &above,
        const std::vector<unsigned int> &below,
        const std::vector<std::pair<double, double>> &representation);
std::vector<std::pair<size_t, double>> get_symmetry_angles(const std::vector<Matrix> &points);
std::vector<double> get_point_angles(const std::vector<Matrix> &points);
std::vector<std::pair<size_t, double>> filter_symmetry_angles(const std::vector<double> &point_angles,
                                                              const std::vector<std::pair<size_t, double>> &symmetry_angles,
                                           const std::vector<std::pair<double, double>> &representation);
std::vector<Matrix> get_representation_with_mid_points(const std::vector<Matrix> &points);
std::vector<std::pair<Matrix, Matrix>> get_symmetry_lines(const std::vector<Matrix> &original_points);