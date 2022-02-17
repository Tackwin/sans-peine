#pragma once

#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>

#include <array>
#include <cmath>
#include <vector>
#include <random>
#include <thread>
#include <optional>

#include <Windows.h>

#include "Constants.hpp"
#include "GUI.hpp"
#include "Physics.hpp"
#include "Definitions.hpp"

constexpr size_t DETAILS[] = { 51, 15, 3 };
constexpr size_t N_DETAILS = sizeof(DETAILS) / sizeof(DETAILS[0]);

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
	double probability_space_size = 0.5;
	size_t probability_resolution = DETAILS[0] * DETAILS[1] * DETAILS[2];
	size_t display_probability_resolution = 250;

	std::vector<Vector2d> estimated_points;

	double last_sps_timestamp = 0;
	size_t curr_sps_counter = 0;
	size_t last_sps_counter = 0;

	std::vector<Driver_Interface> loaded_drivers;
	std::vector<std::thread> driver_threads;
};

extern void render_triangulation(State& state) noexcept;
extern void render_estimation_trace(State& state) noexcept;
extern void update_probability_texture(State& state) noexcept;
