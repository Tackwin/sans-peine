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

	double magnet_strength = 10.f;

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

	size_t distance_resolution = 0;
	size_t angle_resolution = 0;
	std::vector<double> distance_fields;
	std::vector<double> angle_fields;
};

struct Input_State {
	double magnet_strength = 4.875;
	std::array<Beacon, N_Beacons> beacons;
	Reading reading;
};

struct Input_Sampling {
	static constexpr size_t N = 100000;

	Input_State input_state;
	std::array<Vector3d, N_Beacons> (*samplef)(const Input_State&, uint32_t*) = nullptr;
};

struct State;
extern Simulation_Result _space_sim(Simulation_Parameters& state) noexcept;
extern Simulation_Result space_sim(Simulation_Parameters& state) noexcept;
extern void compute_probability_grid(State& state, const Simulation_Result& result) noexcept;

extern Input_Sampling sample_input_space(const Input_State& state) noexcept;
extern void compute_probability_grid(State& state, const Input_Sampling& samplings) noexcept;