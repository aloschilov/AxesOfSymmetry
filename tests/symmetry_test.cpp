#include <symmetry/symmetry.h>
#include <symmetry/matrix.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using testing::Pointwise;
using testing::DoubleNear;


MATCHER_P2(ContainerDoubleNear, container, abs_err, "")
{
    if (container.size() != arg.size()) return false;
    for (int i = 0; i < container.size(); i++)
    {
        if (!::testing::ExplainMatchResult(::testing::DoubleNear(arg[i].first, abs_err), container[i].first, result_listener))
        {
            *result_listener << " for element at idx " << i;
            return false;
        }
        if (!::testing::ExplainMatchResult(::testing::DoubleNear(arg[i].second, abs_err), container[i].second, result_listener))
        {
            *result_listener << " for element at idx " << i;
            return false;
        }
    }
    return true;
}


TEST(SymmetryTests, SequeceWithNodesOnSymmetryLineConvex)
{
    std::vector<RelationToOx> relations = {ABOVE, ON, BELOW, BELOW, ON, ABOVE};
    auto actual = split_above_to_below(relations);

    std::vector<unsigned int> expected_above = {4, 5, 0};
    std::vector<unsigned int> expected_below = {1, 2, 3};
    std::pair<std::vector<unsigned int>, std::vector<unsigned int>> expected{expected_above, expected_below};
    ASSERT_EQ(expected, actual);
}

TEST(SymmetryTests, SequeceWithNodesOnSymmetryLineConvexNoOnNodes)
{
    std::vector<RelationToOx> relations = {BELOW, BELOW, ABOVE, ABOVE};

    auto actual = split_above_to_below(relations);
    std::vector<unsigned int> expected_above = {2, 3};
    std::vector<unsigned int> expected_below = {0, 1};
    std::pair<std::vector<unsigned int>, std::vector<unsigned int>> expected{expected_above, expected_below};
    ASSERT_EQ(expected, actual);
}

TEST(SymmetryTests, CalculateRepresentation1)
{
    const std::vector<Matrix> points = {
            Vector(0, 0),
            Vector(1, 0),
            Vector(1, 1),
            Vector(0, 1)
    };

    std::vector<std::pair<double, double>> actual = calculate_representation(points);
    std::vector<std::pair<double, double>> expected = {
            {1.0, 0}, {1.0, M_PI_2}, {1.0, M_PI}, {1.0, 3.0*M_PI/2.0}
    };
    ASSERT_EQ(expected, actual);
}

TEST(SymmetryTests, CalculateRepresentation2)
{
    const std::vector<Matrix> points = {
            Vector(1, 0),
            Vector(0, sqrt(3)),
            Vector(-1, 0)
    };

    std::vector<std::pair<double, double>> actual = calculate_representation(points);
    std::vector<std::pair<double, double>> expected = {
            {2.0, 2*M_PI/3}, {2.0, 4*M_PI/3}, {2.0, 0}
    };

    EXPECT_THAT(expected, ContainerDoubleNear(actual,EPSILON));
}

TEST(SymmetryTests, CheckSymmetry)
{
    std::vector<unsigned int> above = {0, 1};
    std::vector<unsigned int> below = {2, 3};

    std::vector<std::pair<double, double>> representation = {
            {1.0, M_PI_2}, {1.0, M_PI_2}, {1.0, M_PI_2}, {1.0, M_PI_2}
    };
    bool actual = check_symmetry(above, below, representation);
    bool expected = true;
    ASSERT_EQ(expected, actual);
}

TEST(SymmetryTests, GetSymmetryAngles)
{
    const std::vector<Matrix> points = {
            Vector(-0.5, -0.5),
            Vector(0.0, -0.5),
            Vector(0.5, -0.5),
            Vector(0.5, 0.0),
            Vector(0.5, 0.5),
            Vector(0.0, 0.5),
            Vector(-0.5, 0.5),
            Vector(-0.5, 0.0),
    };
    std::vector<std::pair<size_t,double>> actual = get_symmetry_angles(points);
    std::vector<std::pair<size_t, double>> expected = {std::pair{3, 0}, std::pair{4, M_PI_4}, std::pair{5, M_PI_2}, std::pair{6, 3*M_PI_4}};
    ASSERT_EQ(expected, actual);
}

TEST(SymmetryTests, FilterSymmetries)
{
    const std::vector<Matrix> points = {
            Vector(-0.5, -0.5),
            Vector(0.0, -0.5),
            Vector(0.5, -0.5),
            Vector(0.5, 0.0),
            Vector(0.5, 0.5),
            Vector(0.0, 0.5),
            Vector(-0.5, 0.5),
            Vector(-0.5, 0.0),
    };
    std::vector<std::pair<size_t, double>>  symmetry_angles = {std::pair{3, 0}, std::pair{4, M_PI_4}, std::pair{5, M_PI_2}, std::pair{6, 3*M_PI_4}};
    auto point_angles = get_point_angles(points);
    auto representation = calculate_representation(points);

    auto actual = filter_symmetry_angles(point_angles, symmetry_angles, representation);
    ASSERT_EQ(symmetry_angles, actual);
}

TEST(SymmetryTests, GetRepresentationWithMidPoints)
{
    const std::vector<Matrix> points = {
            Vector(0, 0),
            Vector(1, 0),
            Vector(1, 1),
            Vector(0, 1)
    };
    const std::vector<Matrix> excepted = {
            Vector(0, 0),
            Vector(0.5, 0),
            Vector(1, 0),
            Vector(1, 0.5),
            Vector(1, 1),
            Vector(0.5, 1),
            Vector(0, 1),
            Vector(0, 0.5),
    };
    auto actual = get_representation_with_mid_points(points);
    ASSERT_EQ(excepted, actual);
}

TEST(SymmetryTests, GetSymmetryLinesExample1)
{
    const std::vector<Matrix> points = {
            Vector(0, 0),
            Vector(1, 0),
            Vector(1, 1),
            Vector(0, 1)
    };

    auto actual = get_symmetry_lines(points);

    std::vector<std::pair<Matrix, Matrix>> excepted = {
            {Vector(1, 0.5), Vector(0, 0.5)},
            {Vector(1, 1), Vector(0, 0)},
            {Vector(0.5, 1), Vector(0.5, 0)},
            {Vector(0, 1), Vector(1, 0)},
    };
    ASSERT_EQ(excepted, actual);
}

TEST(SymmetryTests, GetSymmetryLinesExample2)
{
    const std::vector<Matrix> points = {
            Vector(1, 0),
            Vector(0, sqrt(3)),
            Vector(-1, 0)
    };

    auto actual = get_symmetry_lines(points);

    std::vector<std::pair<Matrix, Matrix>> excepted = {
            {(points[0]+points[1])/2.0, Vector(-1, 0)},
            {Vector(0, sqrt(3)), Vector(0, 0)},
            {(points[1]+points[2])/2.0, Vector(1, 0)},
    };
    ASSERT_EQ(excepted, actual);
}

TEST(SymmetryTests, GetSymmetryLinesExample3)
{
    const std::vector<Matrix> points = {
            Vector(0, 0),
            Vector(2, 1),
            Vector(0, 3),
            Vector(-2, 1)
    };

    auto actual = get_symmetry_lines(points);

    std::vector<std::pair<Matrix, Matrix>> excepted = {
            {Vector(0, 3), Vector(0, 0)},
    };
    ASSERT_EQ(excepted, actual);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
