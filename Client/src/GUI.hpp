#pragma once
#include <unordered_map>
#include <SFML/Graphics.hpp>

#include "Constants.hpp"
#include "Physics.hpp"

struct Debug_Values {
	std::unordered_map<const char*, std::vector<double>> histograms;

	void reset() noexcept;
	void add_to_distribution(const char* name, double x) noexcept;
};

extern Debug_Values frame_debug_values;


struct GUI_State {
	double power = 3.0;
	double grid = 0.02;
	sf::Color grid_color = sf::Color(10, 10, 10);

	double magnet_strength = 3.15;
	double magnet_height = 0.06;

	size_t oversampling = 1;
	bool calibrating = false;
	bool display_matrix[N_Beacons * N_Beacons] = { true };

	Simulation_Parameters space_sim;

	bool display_x = true;
	bool display_y = true;
	bool display_z = true;

	bool physic_model = true;

	bool want_compute = false;
	bool want_next_reading = false;

	sf::Texture field_texture;

	size_t sample_to_display = 0;
	bool sample_live = true;

	bool use_dist_beacon[N_Beacons];
	bool use_angle_beacon[N_Beacons];

	bool debug_opened = false;
};

struct State;
extern void render(State& state, GUI_State& gui_state) noexcept;