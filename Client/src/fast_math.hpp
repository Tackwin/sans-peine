#pragma once
inline double fast_atan(double z) {
	auto n1 = 0.97239411f;
	auto n2 = -0.19194795f;
	return (n1 + n2 * z * z) * z;
}

inline double fast_atan2(double y, double x) {
	if (x != 0.0f) {
		if (fabsf(x) > fabsf(y)) {
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