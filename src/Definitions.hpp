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
#pragma pack(pop)

#include <chrono>
static double seconds() noexcept {
	auto now     = std::chrono::system_clock::now();
	auto epoch   = now.time_since_epoch();
	auto seconds = std::chrono::duration_cast<std::chrono::nanoseconds>(epoch);

	// return the number of seconds
	return seconds.count() / 1'000'000'000.0;
}