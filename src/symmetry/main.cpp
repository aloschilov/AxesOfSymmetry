#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>
#include "symmetry/symmetry.h"


std::vector<Matrix> readFromFile(const std::filesystem::path& path)
{
    std::ifstream vertices_file(path);
    std::string file_text, line;

    std::vector<Matrix> res;

    while(getline(vertices_file, line))
    {
        std::stringstream ss(line);

        double value;
        std::vector<double> numbers_as_vector;

        while (ss >> value) {
            numbers_as_vector.push_back(value);
        }

        res.push_back(Vector(numbers_as_vector[0],numbers_as_vector[1]));
    }

    return res;
}

int main(int ac, char* av[]) {
    if(ac == 1) {
        std::cout << "provide the filename as argument" << std::endl;
    } else {
        const std::filesystem::path filepath(av[1]);
        const std::vector<Matrix> points = readFromFile(filepath);
        auto symmetry_lines = get_symmetry_lines(points);


        if(symmetry_lines.empty()) {
            std::cout << "non-symmetric" << std::endl;
        }
        else {
            for(auto const&symmetry_line : symmetry_lines) {
                std::cout <<
                          symmetry_line.first.at(0, 0) << " " << symmetry_line.first.at(1, 0) <<
                          " - "  <<
                          symmetry_line.second.at(0, 0) << " " << symmetry_line.second.at(1, 0) << std::endl;
            }
        }
    }

    return 0;
}
