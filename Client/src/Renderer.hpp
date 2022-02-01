#pragma once

#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>

#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <optional>

#include <Windows.h>

#include "Constants.hpp"
#include "GUI.hpp"
#include "Physics.hpp"
#include "Definitions.hpp"

struct State {
	size_t n_beacons_placed = 0;
	std::array<Beacon, N_Beacons> beacons;
	bool right_clicked = false;

	sf::ContextSettings context_settings;

	sf::RenderTarget* renderTarget = nullptr;
	sf::RenderWindow* window = nullptr;

	HANDLE mail_slot = nullptr;
	bool fullscreen = false;

	std::vector<Reading> readings;
	std::vector<Reading> new_readings;

	float zoom_level = 1.f;
	sf::Vector2f camera_pos = {};

	std::optional<Simulation_Result> space_res;

	size_t next_reading = 0;

	GUI_State gui;

	sf::Texture probability_texture;
	double* probability_grid = nullptr;
	double probability_resolution = 0.001;
	double probability_space_size = 0.5;

	std::vector<Vector2d> estimated_points;

	double last_sps_timestamp = 0;
	size_t curr_sps_counter = 0;
	size_t last_sps_counter = 0;
};

extern void render_triangulation(State& state) noexcept;
extern void render_estimation_trace(State& state) noexcept;
extern void update_probability_texture(State& state) noexcept;
