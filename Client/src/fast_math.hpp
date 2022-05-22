#pragma once
#include <stdint.h>
#include <math.h>

inline double fast_atan(double z) {
	auto n1 = 0.97239411f;
	auto n2 = -0.19194795f;
	return (n1 + n2 * z * z) * z;
}

inline double fast_atan2(double y, double x) {
	if (x != 0.0f) {
		if (fabs(x) > fabs(y)) {
			auto z = y / x;
			if (x > 0.0) {
				// atan2(y,x) = atan(y/x) if x > 0
				return fast_atan(z);
			}
			else if (y >= 0.0) {
				// atan2(y,x) = atan(y/x) + PI if x < 0, y >= 0
				return fast_atan(z) + PI;
			}
			else {
				// atan2(y,x) = atan(y/x) - PI if x < 0, y < 0
				return fast_atan(z) - PI;
			}
		}
		else {
			// Use property atan(y/x) = PI/2 - atan(x/y) if |y/x| > 1.
			auto z = x / y;
			if (y > 0.0) {
				// atan2(y,x) = PI/2 - atan(x/y) if |y/x| > 1, y > 0
				return -fast_atan(z) + PI/2;
			}
			else {
				// atan2(y,x) = -PI/2 - atan(x/y) if |y/x| > 1, y < 0
				return -fast_atan(z) - PI/2;
			}
		}
	}
	else
	{
		if (y > 0.0f) {
			// x = 0, y > 0
			return PI/2;
		}
		else if (y < 0.0f) {
			// x = 0, y < 0
			return -PI/2;
		}
	}
	return 0.0; // x,y = 0. Could return NaN instead.
}

inline double fast_cos(double x) noexcept {
	constexpr double tp = 1./(2.*PI);
	x *= tp;
	x -= .25 + (int)(x + .25);
	x *= 16. * ((x > 0 ? x : -x) - .5);
	return x;
}


uint32_t xorshift32(uint32_t* state)
{
	/* Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs" */
	uint32_t x = *state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	return *state = x;
}

inline double fast_uniform(uint32_t* seed) noexcept {
	uint32_t rng = xorshift32(seed);
	return (double)rng / (double)UINT32_MAX;
}

inline double fast_normal(uint32_t* seed) noexcept {
	double x = 0;
	double y = 0;
	double s = 0;

	do {
		x = fast_uniform(seed) * 2 - 1;
		y = fast_uniform(seed) * 2 - 1;
		s = x + y;
	} while (s < 0 || s > 1);

	return x / sqrt(-2 * log(s) / s);
}