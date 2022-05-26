#pragma once

#include "Constants.hpp"

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
	Vector3d accel[N_Imus] = { 0 };
	Vector3d gyro[N_Imus] = { 0 };
	bool pressed;
	double timestamp = 0;
};

struct Beacon {
	Vector2d pos;

	Vector3d sum_sample = {};
	Vector3d sum2_sample = {};
	double sum_dist = 0;
	double sum2_dist = 0;

	size_t calibration_sample = 0;
	Vector3d mean;
	Vector3d std;
};

struct Pen {
	Vector3d acc_std[N_Imus];
	Vector3d acc_sum[N_Imus];
	Vector3d acc_sum2[N_Imus];
	
	size_t calibration_sample = 0;

	Vector3d current_g = { 0 };

	Vector3d velocity = { 0 };
};

#pragma pack(pop)

#include <chrono>
static double seconds() noexcept {
	auto now     = std::chrono::system_clock::now();
	auto epoch   = now.time_since_epoch();
	auto seconds = std::chrono::duration_cast<std::chrono::nanoseconds>(epoch);

	// return the number of seconds
	return seconds.count() / 1'000'000'000.0;
}

using driver_init_f = void*(*)();
using driver_shut_f = void(*)(void*);

using driver_play_f = void(*)(void*);
using driver_stop_f = void(*)(void*);
using driver_next_f = bool(*)(void*, Reading*);

struct Driver_Interface {
	void* lib = nullptr;
	void* ptr = nullptr;

	driver_init_f init = nullptr;
	driver_shut_f shut = nullptr;
	driver_play_f play = nullptr;
	driver_stop_f stop = nullptr;
	driver_next_f next = nullptr;
};