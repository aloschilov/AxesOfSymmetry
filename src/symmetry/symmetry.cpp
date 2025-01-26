#include <numeric>
#include <tuple>
#include "symmetry/symmetry.h"
#include "symmetry/matrix.h"

std::pair<std::vector<unsigned int>, std::vector<unsigned int>>
split_above_to_below(const std::vector<RelationToOx> &relations) {
    std::vector<unsigned int> above_indexes;
    std::vector<unsigned int> below_indexes;

    size_t n = relations.size();
    size_t on_or_below_to_above_index = 0;
    size_t above_to_on_or_below = 0;

    for (int i = 0; i < n; ++i) {
        if ((relations[i] == ON || relations[i] == BELOW) && relations[(i + 1) % n] == ABOVE) {
            if (relations[i] == ON)
                on_or_below_to_above_index = (i) % n;
            else
                on_or_below_to_above_index = (i + 1) % n;
        }

        if (relations[i] == ABOVE && (relations[(i + 1) % n] == ON || relations[(i + 1) % n] == BELOW)) {
            above_to_on_or_below = (i + 1) % n;
        }
    }

    for (size_t i = on_or_below_to_above_index; i % n != above_to_on_or_below; ++i) {
        above_indexes.push_back(i % n);
    }

    for (size_t i = above_to_on_or_below; i % n != on_or_below_to_above_index; ++i) {
        below_indexes.push_back(i % n);
    }

    return {above_indexes, below_indexes};
}

double normalize_angle(double angle) {
    if (angle < 0) {
        return 2 * M_PI + angle;
    } else if (angle >= 2 * M_PI) {
        return angle - 2 * M_PI;
    } else {
        return angle;
    }
}

double atan(const Matrix &p) {
    return atan2(p.at(1, 0), p.at(0, 0));
}

std::vector<std::pair<double, double>> calculate_representation(const std::vector<Matrix> &points) {
    std::vector<std::pair<double, double>> res;
    const unsigned int &n = points.size();

    for (int i = 0; i < points.size(); i++) {
        Matrix edge_vector = points[(i + 1) % points.size()] - points[i % points.size()];
        Matrix next_edge_vector = points[(i + 2) % points.size()] - points[(i + 1) % points.size()];

        double &x = edge_vector[0];
        double &y = edge_vector[1];
        double alpha = normalize_angle(atan(edge_vector));

        const double length = sqrt(x * x + y * y);

        res.emplace_back(length, alpha);
    }

    return res;
}

std::vector<std::pair<double, double>>
rotate_representation(const std::vector<std::pair<double, double>> &representation,
                      double angle) {
    std::vector<std::pair<double, double>> res;

    for (auto r: representation) {
        res.emplace_back(r.first, normalize_angle(r.second - angle));
    }

    return res;
}

bool check_symmetry(
        const std::vector<unsigned int> &above,
        const std::vector<unsigned int> &below,
        const std::vector<std::pair<double, double>> &representation) {

    const unsigned int n = above.size();
    if (above.size() == below.size()) {
        for (unsigned int i = 0; i < n; ++i) {
            const std::pair<double, double> &r_above = representation[above[i]];
            const std::pair<double, double> &r_below = representation[below[n - i - 1]];

            double rep_above_d = r_above.first;
            double rep_above_a = r_above.second;
            double rep_below_d = r_below.first;
            double rep_below_a = normalize_angle(M_PI - r_below.second);

            if (!(fabs(rep_above_d - rep_below_d) < EPSILON && fabs(rep_above_a - rep_below_a) < EPSILON)) {
                return false;
            }
        }
        return true;
    } else {
        return false;
    }
}

std::vector<std::pair<size_t, double>> get_symmetry_angles(const std::vector<Matrix> &points) {
    std::vector<std::pair<size_t, double>> symmetry_angles;

    const unsigned int n = points.size();

    for (size_t i = 0; i < n; i++) {

        const Matrix &p_i = points[i];
        const double &point_symmetry_angle = normalize_angle(atan(p_i));

        if (point_symmetry_angle < M_PI) {
//            TODO
            symmetry_angles.emplace_back(i, point_symmetry_angle);
        }
    }
    return symmetry_angles;
}

//////////////////////////////////////////////////////////////////////////////
/// @brief This function tests position of point relation to ox
///
/// @param angle - input - angle w.r.t. ox
///
/// @returns - ON, ABOVE or BELOW
///
//////////////////////////////////////////////////////////////////////////////
RelationToOx angle_above_ox(double angle) {
    double normalized_angle = normalize_angle(angle);
    if (fabs(normalized_angle - 0.0) < EPSILON || fabs(normalized_angle - M_PI) < EPSILON) {
        return ON;
    } else {
        if (normalized_angle < M_PI) {
            return ABOVE;
        } else {
            return BELOW;
        }
    }
}

std::vector<double> get_point_angles(const std::vector<Matrix> &points) {
    std::vector<double> point_angles;

    for (const Matrix &point: points) {
        const double &point_angle = normalize_angle(atan(point));

        point_angles.push_back(point_angle);
    }
    return point_angles;
}

std::vector<RelationToOx> calculate_rtx(const std::vector<double> &point_angles, double symmetry_angle) {
    std::vector<RelationToOx> rtx;

    for (double point_angle: point_angles) {
        const double diff = point_angle - symmetry_angle;

        RelationToOx a_on_b = angle_above_ox(diff);
        rtx.push_back(a_on_b);
    }
    return rtx;
}

std::vector<std::pair<size_t, double>> filter_symmetry_angles(const std::vector<double> &point_angles,
                                                              const std::vector<std::pair<size_t, double>> &symmetry_angles,
                                                              const std::vector<std::pair<double, double>> &representation) {
    std::vector<std::pair<size_t, double>> res;
    for (const std::pair<size_t, double> &symmetry_angle: symmetry_angles) {

        auto rtx = calculate_rtx(point_angles, symmetry_angle.second);

        auto above_and_below = split_above_to_below(rtx);

        auto rotated_representation = rotate_representation(representation, symmetry_angle.second);
        bool is_symmetric = check_symmetry(above_and_below.first, above_and_below.second, rotated_representation);
        if (is_symmetric)
            res.push_back(symmetry_angle);
    }
    return res;
}

std::vector<Matrix> get_representation_with_mid_points(const std::vector<Matrix> &points) {
    const size_t n = points.size();
    std::vector<Matrix> res;

    for (int i = 0; i < n; i++) {
        const Matrix &p_i = points[i];
        const Matrix &p_i_1 = points[(i + 1) % n];
        const Matrix mid_edge_point = (p_i_1 + p_i) / 2.0;
        res.push_back(p_i);
        res.push_back(mid_edge_point);
    }
    return res;
}

std::vector<std::pair<Matrix, Matrix>> get_symmetry_lines(const std::vector<Matrix> &original_points) {
    std::vector<std::pair<Matrix, Matrix>> res;

    Matrix centroid = std::accumulate(original_points.begin(), original_points.end(), Vector(0, 0)) /
                      double(original_points.size());


    auto representation_with_mid_points = get_representation_with_mid_points(original_points);

    std::vector<Matrix> points;

    for (const Matrix &representation_with_mid_point: representation_with_mid_points) {
        points.push_back(representation_with_mid_point - centroid);
    }

    const size_t n = points.size();
    auto symmetry_angles = get_symmetry_angles(points);
    auto representation = calculate_representation(points);
    auto point_angles = get_point_angles(points);

    auto angles_with_symmetries = filter_symmetry_angles(
            point_angles,
            symmetry_angles,
            representation);
    for (auto angles_with_symmetry: angles_with_symmetries) {
        size_t i = angles_with_symmetry.first;
        res.emplace_back(representation_with_mid_points[i], representation_with_mid_points[(i + n / 2) % (n)]);
    }
    return res;
}