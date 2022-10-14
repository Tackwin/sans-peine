#include "Physics.hpp"
#include <cmath>
#include "Renderer.hpp"
#include <cmath>
#include <omp.h>

#include "Macro.hpp"

#include "fast_math.hpp"

extern Simulation_Result _space_sim(Simulation_Parameters& state) noexcept {
	Simulation_Result result;
#if 0
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

#endif
	return result;
}


Simulation_Result space_sim(Simulation_Parameters& state) noexcept {
	// >SPEED(Tackwin): result allocate the space for the result itself but it's stupid we
	// could pass the space as a parameter to save allocation since we likely want to reuse the
	// space for new readings.
	Simulation_Result result;
	
	auto read = state.reading;

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

	result.distance_fields.resize(state.distance_resolution * N_Beacons, 1);
	result.angle_fields.resize(state.angle_resolution * N_Beacons, 1);

	for (size_t b_idx = 0; b_idx < N_Beacons; ++b_idx) if (state.use_dist[b_idx]) {
		auto& b = state.beacons[b_idx];
		auto n = b.calibration_sample;

		auto s = std::sqrt((b.sum2_dist / n - (b.sum_dist / n) * (b.sum_dist / n)) * (n / (n - 1.0)));
		auto u = b.sum_dist / n;

		u = std::hypot(read.beacons[b_idx].x, read.beacons[b_idx].y, read.beacons[b_idx].z) - u;
		u = std::abs(u);

		for (size_t i = 1; i < state.distance_resolution; ++i) {
			auto d = i * state.distance_step;

			auto lambda = 0.0;
			lambda += (read.beacons[b_idx].x / b.std.x) * (read.beacons[b_idx].x / b.std.x);
			lambda += (read.beacons[b_idx].y / b.std.y) * (read.beacons[b_idx].y / b.std.y);

			auto c = state.magnet_strength * u0 / (4 * PI);

			auto p = N(c / (d * d * d), u, s) * 3 * c / (d * d * d *d);
			result.distance_fields[i + b_idx * state.distance_resolution] = p;
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

	for (size_t b_idx = 0; b_idx < N_Beacons; ++b_idx) if (state.use_angle[b_idx]) {
		auto& b = state.beacons[b_idx];
		auto n = b.calibration_sample;

		auto ux = read.beacons[b_idx].x - b.mean.x;
		auto uy = read.beacons[b_idx].y - b.mean.y;

		auto sx = b.std.x + std::abs(ux) * state.sensitivity;
		auto sy = b.std.y + std::abs(uy) * state.sensitivity;

		auto sx2 = sx*sx;
		auto sy2 = sy*sy;

		for (size_t i = 0; i < state.angle_resolution; ++i) {

			auto z = PI * 2 * i / (state.angle_resolution - 1.0);
			auto tz = std::tan(z);

			long double d = 2 * sx2 * sy2;
			long double a = (sy2 + sx2 * tz * tz) / d;
			long double b = 2 * (sy2 * ux + sx2 * tz * uy) / d;
			long double c = (ux*ux * sy2 + sx2 * uy*uy) / d;


			long double f = 1.0 / (std::cos(z) * std::cos(z) * sx*sy*2*PI);

			#if 1
			long double p = std::exp(-c) * f;
			p *= (
				2 * std::sqrt(a) +
				std::sqrt(PI) * b * std::exp(b*b / (4*a)) * std::erf(b/(2*std::sqrt(a)))
			) / (2 * std::pow(a, 1.5));
			#elif 0

			auto bsqrta2 = b / (2 * std::sqrt(a));

			auto p = std::log(f);
			p -= c;
			p -= std::log(2 / std::sqrt(PI));
			p -= 1.5 * std::log(a);
			p += bsqrta2 * bsqrta2;
			p += std::log(b * std::erf(bsqrta2));

			#endif

			result.angle_fields[i + b_idx * state.angle_resolution] = p;
		}
	}

	result.distance_resolution = state.distance_resolution;
	result.angle_resolution = state.angle_resolution;
	return result;
}
extern Debug_Values perf_values;

void compute_probability_grid(State& state, const Simulation_Result& result) noexcept {
	auto w = state.probability_resolution;
	auto h = state.probability_resolution;
	auto distance_freq = 1.0 / result.input_parameters.distance_step;
	auto angle_freq = result.input_parameters.angle_resolution / (2 * PI);
	auto z = result.input_parameters.h * result.input_parameters.h / 4;
	auto ww = state.probability_space_size / (w - 1.0);
	auto hh = state.probability_space_size / (h - 1.0);
	auto bias = -state.probability_space_size / 2.0;

	// for (size_t i = 0; i < w * h; ++i) state.probability_grid[i] = 0;

	// auto f = [&] (size_t xi, size_t yi) {
	// 	auto px = xi * ww + bias;
	// 	auto py = yi * hh + bias;
	// 	double local_p = 1;

	// 	for (size_t b_idx = 0; b_idx < N_Beacons; ++b_idx) {
	// 		auto& b = state.beacons[b_idx];
	// 		auto x = px - b.pos.x;
	// 		auto y = py - b.pos.y;

	// 		auto d = std::sqrt(x * x + y * y + z);

	// 		int64_t t_low  = (int64_t)(d * distance_freq);
	// 		double t = d * distance_freq - t_low;

	// 		auto u = result.distance_fields[t_low + 0][b_idx];
	// 		auto v = result.distance_fields[t_low + 1][b_idx];

	// 		local_p *= u * (1 - t) + v * t;

	// 		double a = PI/2 - fast_atan2(y, x);
	// 		if (a < 0) a += 2 * PI;

	// 		t_low  = (int64_t)(a * angle_freq);
	// 		t = a * angle_freq - t_low;

	// 		u = result.angle_fields[t_low + 0][b_idx];
	// 		v = result.angle_fields[t_low + 1][b_idx];

	// 		local_p *= u * (1 - t) + v * t;
	// 	}
	// 	state.probability_grid[xi + yi * w] = std::max(0.0, std::log(local_p)) + rand() / (1.0 * RAND_MAX);

	// 	return local_p;
	// };


	// auto F = [&] (size_t lx, size_t hx, size_t ly, size_t hy, size_t d, auto F) {
	// 	size_t D = DETAILS[d];
	// 	std::array<double, DETAILS[0] * DETAILS[0]> current_level_samples;

	// 	auto at = [&] (size_t x, size_t y) -> double& { return current_level_samples[x + y * D]; };

	// 	for (size_t xi = lx; xi < hx; xi += (hx - lx) / D)
	// 	for (size_t yi = ly; yi < hy; yi += (hy - ly) / D) {
	// 		at(D * (xi - lx) / (hx - lx), D * (yi - ly) / (hy - ly)) = f(xi, yi);

	// 	}

	// 	if (d + 1 >= N_DETAILS) return;

	// 	auto deeper = [&] (size_t x, size_t y) {
	// 		F(
	// 			lx + x * (hx - lx) / D,
	// 			lx + (x + 1) * (hx - lx) / D,
	// 			ly + y * (hy - ly) / D,
	// 			ly + (y + 1) * (hy - ly) / D,
	// 			d + 1,
	// 			F
	// 		);
	// 	};

	// 	bool went_deep[DETAILS[0] * DETAILS[0]] = { false };
	// 	for (size_t x = 0; x < D; x++) for (size_t y = 0; y < D; y++)
	// 	{
	// 		auto a = at(x, y);
	// 		if (a > state.gui.epsilon) {
	// 			deeper(x, y);
	// 		}
	// 	}
	// };

	// F(0, w, 0, h, 0, F);

	// return;


	for (size_t x = 0; x < w * h; ++x) state.probability_grid[x] = 1;

	size_t idx = 0;
	#pragma omp parallel for
	for (size_t yi = 0; yi < h; ++yi)
	for (size_t xi = 0; xi < w; ++xi)
	{
		auto px = xi * ww + bias;
		auto py = yi * hh + bias;

		double local_p = 1;

		for (size_t b_idx = 0; b_idx < N_Beacons; ++b_idx) {
			auto& b = state.beacons[b_idx];
			auto x = px - b.pos.x;
			auto y = py - b.pos.y;

			auto d = std::sqrt(x * x + y * y + z);

			int64_t t_low  = (int64_t)(d * distance_freq);
			double t = d * distance_freq - t_low;

			auto u = result.distance_fields[t_low + 0 + b_idx * result.distance_resolution];
			auto v = result.distance_fields[t_low + 1 + b_idx * result.distance_resolution];

			local_p *= u * (1 - t) + v * t;

			double a = std::atan2(y, x);
			if (a < 0) a += 2 * PI;

			t_low  = (int64_t)(a * angle_freq);
			t = a * angle_freq - t_low;

			u = result.angle_fields[t_low + 0 + b_idx * result.angle_resolution];
			v = result.angle_fields[t_low + 1 + b_idx * result.angle_resolution];

			local_p *= u * (1 - t) + v * t;
		}

		state.probability_grid[xi + yi * w] *= local_p;
	}
}


Input_Sampling sample_input_space(const Input_State& input_state) noexcept {
	Input_Sampling result;

	result.input_state = input_state;
	result.samplef = [] (const Input_State& input_state, uint32_t* seed) {
		Input_Sampling::Sample sample;
		for (size_t i = 0; i < N_Beacons; ++i) {
			auto& r = input_state.reading.beacons[i];
			auto& b = input_state.beacons[i];
			Vector3d mag;

			mag.x = (r.x - b.mean.x) + b.std.x * fast_normal(seed);
			mag.y = (r.y - b.mean.y) + b.std.y * fast_normal(seed);
			mag.z = (r.z - b.mean.z) + b.std.z * fast_normal(seed);

			sample.mag[i] = mag;
		}

		const Pen& pen = input_state.pen;
		sample.acc = { 0 };
		for (size_t i = 0; i < N_Imus; ++i) {
			auto& a = input_state.reading.accel[i];
			sample.acc.x += (a.x - input_state.current_g.x) + pen.acc_std[i].x * fast_normal(seed);
			sample.acc.y += (a.y - input_state.current_g.y) + pen.acc_std[i].y * fast_normal(seed);
			sample.acc.z += (a.z - input_state.current_g.z) + pen.acc_std[i].z * fast_normal(seed);
		}
		sample.acc /= N_Imus;

		return sample;
	};

	return result;
}

void compute_probability_grid(State& state, const Input_Sampling& samplings) noexcept {
	auto w = state.probability_resolution;
	auto h = state.probability_resolution;
	auto ww = state.probability_space_size / (w - 1.0);
	auto hh = state.probability_space_size / (h - 1.0);
	auto bias = -state.probability_space_size / 2.0;

	auto m_strength = samplings.input_state.magnet_strength;

	for (size_t x = 0; x < w * h; ++x) state.probability_grid[x] = 0;

	auto solve_xy = [] (Vector3d m, Vector3d b) -> Vector2d {
		thread_local uint32_t seed[1] = { 2 };
		constexpr double u0 = 1.225663753e-6;
		constexpr double PI = 3.141592653589793238462643383279502884;

		auto c = (3 * u0) / (4 * PI);

		auto mx = m.x;
		auto my = m.y;
		auto mz = m.z;

		auto bx = b.x;
		auto by = b.y;
		auto bz = b.z;

		auto bx2 = bx * bx;
		auto by2 = by * by;
		auto bz2 = bz * bz;

		auto mx2 = mx * mx;
		auto my2 = my * my;
		auto mz2 = mz * mz;

		auto s7 = bx2 + by2 + bz2;
		auto s6 = mz * bx2 - bz * mx * bx + mz * by2 - bz * my * by;
		s6 *= s6;
		auto s5 = pow(abs(-c * (bx * mx + by * my + bz * mz) / s7), 2.0/3.0);
		auto s4 = bx2 * (my2 + mz2) - 2 * bx * mx * (by * my + bz * mz);
		s4     += by2 * (mx2 + mz2) - 2 * by * bz * my * mz + bz2 * (mx2 + my2);

		auto x = (mx * by2 - bx * my * by + mx * bz2 - bx * mz * bz) * std::sqrt(s7 * s6 * s5 * s4);
		x /= s7 * s4 * (mz * bx2 - bz * mx * bx + mz * by2 - bz * my * by);

		auto y = (my * bx2 - by * mx * bx + my * bz2 - by * mz * bz) * std::sqrt(s5 * s7 * s6 * s4);
		y /= s7 * s4 * (mz * bx2 - bz * mx * bx + mz * by2 - bz * my * by);


		// auto l = std::hypot(b.x, b.y, b.z);
		// l = m.z * 1e-6 / l;
		// l = std::pow(l, 1.0/3.0);
		// auto p = fast_uniform(seed) * 2 * 3.1415926;
		// auto x = std::cos(p) * l;
		// auto y = std::sin(p) * l;

		return Vector2d{ +x, +y };
	};

	uint32_t seeds[32];
	for (size_t i = 0; i < 32; ++i) seeds[i] = (int32_t)(i + 1);

	#pragma omp parallel for
	for (size_t n = 0; n < samplings.N; ++n) {
		uint32_t* seed = seeds + 1; // omp_get_thread_num();
		Input_Sampling::Sample sample = samplings.samplef(samplings.input_state, seed);

		if (state.gui.compute_mag) for (size_t i = 0; i < N_Beacons; ++i) {
			auto b = sample.mag[i];
			auto m = Vector3d{0.0, 0.0, m_strength};

			auto p = solve_xy(m, b);

			p.x += state.beacons[i].pos.x;
			p.y += state.beacons[i].pos.y;

			int xi = (p.x - bias) / ww;
			int yi = (p.y - bias) / ww;

			if (
				0 <= xi && xi < state.probability_resolution &&
				0 <= yi && yi < state.probability_resolution
			) {
				state.probability_grid[xi + yi * w] += 1;
			}
		}

		if (state.gui.compute_acc) {
			Vector3d v = samplings.input_state.pen.velocity;
			v.x += sample.acc.x * samplings.input_state.dt;
			v.y += sample.acc.y * samplings.input_state.dt;
			v.z += sample.acc.z * samplings.input_state.dt;

			double V = std::hypot(v.x, v.y, v.z);

			double a = 2 * 3.1415926 * n / samplings.N;
			int xi = (V * std::cos(a) - bias) / ww;
			int yi = (V * std::sin(a) - bias) / ww;

			if (
				0 <= xi && xi < state.probability_resolution &&
				0 <= yi && yi < state.probability_resolution
			) {
				state.probability_grid[xi + yi * w] += 1;
			}
		}
	}
}
