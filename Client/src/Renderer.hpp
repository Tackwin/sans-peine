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

constexpr const char* Mail_Name = "\\\\.\\Mailslot\\SP";
constexpr const char* App_Name = "SP Client";


struct Vector3d {
	double x, y, z = 0;
};

static Vector3d operator+(Vector3d a, Vector3d b) noexcept {
	return { a.x + b.x, a.y + b.y, a.z + b.z };
}

static Vector3d& operator+=(Vector3d& a, Vector3d b) noexcept {
	a = a + b;
	return a;
}
static Vector3d& operator/=(Vector3d& a, double n) noexcept {
	a.x /= n;
	a.y /= n;
	a.z /= n;
	return a;
}
#pragma pack(push, 1)
struct Reading {
	Vector3d beacons[N_Beacons] = { 0 };
	bool pressed;
};

struct Beacon {
	sf::Vector2f pos;

	Vector3d sum_sample = {};
	size_t calibration_sample = 0;
};
#pragma pack(pop)

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

	float zoom_level = 1.f;
	sf::Vector2f camera_pos = {};

	std::optional<Simulation_Result> space_res;

	size_t next_reading = 0;

	GUI_State gui;
};

extern void render_triangulation(State& state) noexcept;
