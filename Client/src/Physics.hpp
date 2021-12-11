#pragma once

#include <vector>

#include "Constants.hpp"

struct Simulation_Result {
	size_t resolution = 0;
	std::vector<double> field;
};
struct Simulation_Parameters {
	size_t resolution = 1000;

	double pivot_distance = 0.00;
	double pivot_angle = PI / 2;
	double R_in = 0.000;
	double R_out = 0.005;
	double h = 0.065;

	double magnet_strength = 4.875;
};

extern Simulation_Result space_sim(Simulation_Parameters& state) noexcept;
