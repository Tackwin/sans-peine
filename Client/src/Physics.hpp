#pragma once

#include <array>
#include <vector>

#include "Definitions.hpp"
#include "Constants.hpp"

struct Simulation_Parameters {
	size_t resolution = 1000;

	double pivot_distance = 0.00;
	double pivot_angle = PI / 2;
	double R_in = 0.000;
	double R_out = 0.005;

	double h = 0.065;

	double magnet_strength = 4.875;

	size_t distance_resolution = 1000;
	double distance_step = 0.001;

	size_t angle_resolution = 1000;

	Reading reading;
	std::array<Beacon, N_Beacons> beacons;

	std::array<bool, N_Beacons> use_dist;
	std::array<bool, N_Beacons> use_angle;

	double sensitivity = 0.1;
};


struct Simulation_Result {
	Simulation_Parameters input_parameters;

	std::vector<std::array<double, N_Beacons>> distance_fields;
	std::vector<std::array<double, N_Beacons>> angle_fields;
};

struct State;
extern Simulation_Result _space_sim(Simulation_Parameters& state) noexcept;
extern Simulation_Result space_sim(Simulation_Parameters& state) noexcept;
extern void compute_probability_grid(State& state, const Simulation_Result& result) noexcept;