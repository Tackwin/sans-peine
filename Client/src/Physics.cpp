#include "Physics.hpp"
#include <cmath>
#include "Renderer.hpp"
#include <cmath>

extern Simulation_Result space_sim(Simulation_Parameters& state) noexcept {
	double h = state.h;
	double dist = state.pivot_distance + h / 2;

	auto C = [=](double kc, double p, double c, double s) -> double {
		if (kc == 0) return NAN;

		double errtol = 0.0001;
		double k = kc > 0 ? kc : -kc;
		double em = 1;

		if (p > 0) {
			p = std::sqrt(p);
			s = s / p;
		} else {
			double f = kc * kc;
			double q = 1 - f;
			double g = 1 - p;
			f = f - p;
			q = q * (s - c * p);
			p = std::sqrt(f / g);
			c = (c - s) / g;
			s = -q / (g * g * p) + c * p;
		}

		double f = c;
		c = c + s / p;
		double g = k / p;
		s = 2 * (s + f * g);
		p = g + p;
		g = em;
		em = k + em;
		double kk = k;

		while (std::abs(g - k) > g * errtol) {
			k = 2 * std::sqrt(kk);
			kk = k * em;
			f = c;
			c = c + s / p;
			g = kk / p;
			s = 2 * (s + f * g);
			p = g + p;
			g = em;
			em = k + em;
		}

		return PI / 2 * (s + c * em) / (em * (em + p));
	};

	double B0 = state.magnet_strength;
	double b = state.h / 2;

	struct Cylindrical_Vec {
		double p;
		double z;
	};
	auto at = [=](double a, double p, double z) -> Cylindrical_Vec {
		double zp = z + b;
		double zm = z - b;
		double ap = a / std::sqrt(zp * zp + (p + a) * (p + a));
		double am = a / std::sqrt(zm * zm + (p + a) * (p + a));
		double bp = zp / std::sqrt(zp * zp + (p + a) * (p + a));
		double bm = zm / std::sqrt(zm * zm + (p + a) * (p + a));
		double g = (a - p) / (a + p);
		double kp = std::sqrt(((zp * zp) + (a - p) * (a - p)) / ((zp * zp) + (a + p) * (a + p)));
		double km = std::sqrt(((zm * zm) + (a - p) * (a - p)) / ((zm * zm) + (a + p) * (a + p)));

		double Bp = B0 * (ap * C(kp, 1, 1, -1) - am * C(km, 1, 1, -1));
		double Bz = B0 * a / (a + p) * (bp * C(kp, g *g, 1, g) - bm * C(km, g * g, 1, g));

		return { Bp, Bz };
	};


	Simulation_Result result;

	result.field.reserve((state.resolution + 1) * (state.resolution + 1) * 6);
	result.field.resize((state.resolution + 1) * (state.resolution + 1) * 3);

	double* gridx = result.field.data();
	double* gridy = gridx + (state.resolution + 1) * (state.resolution + 1);
	double* gridz = gridy + (state.resolution + 1) * (state.resolution + 1);

	double max_ = 0;
	double min_ = 0;

	for (size_t x = 0; x <= state.resolution; x++) for (size_t y = 0; y <= state.resolution; y++) {
		// the point to probe in cartesian space
		double px = (double)x / state.resolution;
		double py = (double)y / state.resolution - 1;
		double pz = 0;
		px *= 0.26; // in meter
		py *= 0.26; // in meter
		pz *= 0.26; // in meter

		// The center of the magnet in cartesian space
		double cx = dist * std::cos(state.pivot_angle);
		double cy = 0;
		double cz = dist * std::sin(state.pivot_angle);

		// the axis of the magnet in cartesian space
		double vx = (dist - 1) * std::cos(state.pivot_angle) - cx;
		double vy = 0;
		double vz = (dist - 1) * std::sin(state.pivot_angle) - cz;

		// The axis normalized
		double l = std::sqrt(vx * vx + vy * vy + vz * vz);
		vx /= l;
		vy /= l;
		vz /= l;

		// The vector from the center of the magnet to the point to probe in cartesian space
		double dx = px - cx;
		double dy = py - cy;
		double dz = pz - cz;
		l = std::sqrt(dx * dx + dy * dy + dz * dz);

		// z is the dot product of the axis and the delta vector
		double z = dx * vx + dy * vy + dz * vz;

		// p is the magnitude of the delta vector when the z component is removed
		dx -= z * vx;
		dy -= z * vy;
		dz -= z * vz;
		double p = std::sqrt(dx * dx + dy * dy + dz * dz);

		auto [rout_p, rout_z] = at(state.R_out, p, z);
		auto [rin_p,  rin_z]  = at(state.R_in, p, z);

		double resx = vx * (rout_z - rin_z);
		double resy = vy * (rout_z - rin_z);
		double resz = vz * (rout_z - rin_z);
		
		resx += dx * (rout_p - rin_p);
		resy += dy * (rout_p - rin_p);
		resz += dz * (rout_p - rin_p);

		gridx[x + y * (state.resolution + 1)] = resx;
		gridy[x + y * (state.resolution + 1)] = resy;
		gridz[x + y * (state.resolution + 1)] = resz;

		// grid[x + y * (state.resolution + 1)] = at(state.R_out, p, z) - at(state.R_in, p, z);
	}

	result.resolution = state.resolution;

	return result;
}

void compute_probability_grid(State& state) noexcept {
	constexpr double u0 = 1.225663753e-6;
	constexpr double PI = 3.141592653589793238462643383279502884;


	if (state.readings.size() < 2) return;
	auto to_draw = state.gui.sample_live ? state.readings.size() : state.gui.sample_to_display;
	auto read = state.readings[to_draw - 1];

	auto w = (size_t)(state.probability_space_size / state.probability_resolution);
	auto h = (size_t)(state.probability_space_size / state.probability_resolution);

	auto N = [] (double x, double u, double s) -> double {
		return std::exp(-0.5 * ((x - u) / s) * ((x - u) / s)) / (s * std::sqrt(2 * PI));
	};


	auto I = [] (double x, double a) -> double {
		double res = 0.0;
		double nomi = x*x/4;
		double current_nomi = 1.0;

		double denomi = 1;

		for (size_t j = 0; j < 5; ++j) {
			res += current_nomi / (denomi * std::tgamma(a + j + 1));

			current_nomi *= nomi;
			denomi *= (j + 1);
		}

		return res * std::pow(x / 2, a);
	};

	auto QI = [I] (double x, double lambda, double k) -> double {
		return
			1/2.0 *
			std::exp(-(x + lambda) / 2) *
			std::pow(x / lambda, k / 4 - 0.5) *
			I(std::sqrt(lambda * x), k / 2 - 1);
	};

	for (size_t x = 0; x < w * h; ++x) state.probability_grid[x] = 1;

	for (size_t b_idx = 0; b_idx < N_Beacons; ++b_idx) if (state.gui.use_dist_beacon[b_idx]) {
		auto& b = state.beacons[b_idx];
		auto n = b.calibration_sample;

		auto s = std::sqrt((b.sum2_dist / n - (b.sum_dist / n) * (b.sum_dist / n)) * (n / (n - 1.0)));
		s *= 100;
		auto u = b.sum_dist / n;

		u = std::hypot(read.beacons[b_idx].x, read.beacons[b_idx].y, read.beacons[b_idx].z) - u;
		u = std::abs(u);
		
		for (size_t xi = 0; xi < w; ++xi) for (size_t yi = 0; yi < h; ++yi) {
			auto px = xi / (w - 1.0) - 0.5;
			auto py = yi / (h - 1.0) - 0.5;

			px *= state.probability_space_size * 1;
			py *= state.probability_space_size * 1;

			px -= b.pos.x;
			py -= b.pos.y;

			auto d = std::hypot(px, py, state.gui.magnet_height / 2);
			auto lambda = 0.0;
			lambda += (read.beacons[b_idx].x / b.std.x) * (read.beacons[b_idx].x / b.std.x);
			lambda += (read.beacons[b_idx].y / b.std.y) * (read.beacons[b_idx].y / b.std.y);

			auto c = state.gui.magnet_strength * u0 / (4 * PI);

			auto p = N(c / (d * d * d), u, s) * 3 * c / (d * d * d *d);
			state.probability_grid[xi + yi * w] *= p;
		}
	}

	auto Q1 = [&] (double u) -> double { return std::sqrt(PI) / 2.0 * std::erfc(u); };
	auto Q2 = [&] (double x, double y, double p) -> double {
		auto delta = (x >= 0 && y >= 0) ? 0 : 0.5;

		auto help = [&] (double x, double p) -> double {
			auto Q = [&] (double x, double p) -> double {
				auto res = 0;
				res += std::sqrt(1 - p * p) / 12 * Q1(x / std::sqrt(1 - p * p));

				res +=
					0.25 *
					std::sqrt((3 - 3 * p * p) / (3 + p * p)) *
					Q1(x * std::sqrt((3 + p * p) / (3 - 3 * p * p)));

				return res;
			};

			if (x < 0) {
				if (p < 0) return 0.5 - (Q1(-x) - Q(-x, p));
				return 0.5 - Q(-x, -p);
			} else {
				if (p < 0) return Q(x, p);
				return Q1(x) - Q(x, -p);
			}
		};

		auto signx = x >= 0 ? 1 : -1;
		auto signy = y >= 0 ? 1 : -1;
		auto px = signx * (p * x - y) / std::sqrt(x*x - 2 * p * x * y + y * y);
		auto py = signy * (p * y - x) / std::sqrt(x*x - 2 * p * x * y + y * y);


		return help(x, px) + help(y, py) - delta;
	};

	for (size_t b_idx = 0; b_idx < N_Beacons; ++b_idx) if (state.gui.use_angle_beacon[b_idx]) {
		auto& b = state.beacons[b_idx];
		auto n = b.calibration_sample;

		auto sx = b.std.x * 100;
		auto sy = b.std.x * 100;

		auto sx2 = sx*sx;
		auto sy2 = sy*sy;

		auto ux = read.beacons[b_idx].x - b.mean.x;
		auto uy = read.beacons[b_idx].y - b.mean.y;

		for (size_t xi = 0; xi < w; ++xi) for (size_t yi = 0; yi < h; ++yi) {
			auto px = xi / (w - 1.0) - 0.5;
			auto py = yi / (h - 1.0) - 0.5;

			px *= state.probability_space_size * 1;
			py *= state.probability_space_size * 1;

			px -= b.pos.x;
			py -= b.pos.y;

			auto z = PI/2 -std::atan(py / px);
			auto tz = std::tan(z);

			auto d = 2 * sx2 * sy2;
			auto a = (sy2 + sx2 * tz * tz) / d;
			auto b = 2 * (sy2 * ux + sx2 * tz * uy) / d;
			auto c = -(ux*ux * sy2 + sx2 * uy*uy) / d;


			auto f = 1.0 / (std::cos(z) * std::cos(z) * sx*sy*2*PI);
			auto p = std::exp(-c) * f;
			p *= (
				2 * std::sqrt(a) +
				std::sqrt(PI) * b * std::exp(b*b / (4*a)) * std::erf(b/(2*std::sqrt(a)))
			) / (2 * std::pow(a, 1.5));
			if (b_idx == 0) {
				frame_debug_values.add_to_distribution("a", a);
				frame_debug_values.add_to_distribution("b", b);
				frame_debug_values.add_to_distribution("p", p);
				frame_debug_values.add_to_distribution("e", 2 * std::sqrt(a) +
				std::sqrt(PI) * b * std::exp(b*b / (4*a)) * std::erf(b/(2*std::sqrt(a))));
			}
			state.probability_grid[xi + yi * w] *= p;
		}
	}
}
