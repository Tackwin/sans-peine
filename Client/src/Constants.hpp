#pragma once

constexpr double PI = 3.141592653589793238462643383279502884;
constexpr size_t N_Beacons = 4;


template<typename T>
struct array_view {
	T* data = nullptr;
	size_t count = 0;
};