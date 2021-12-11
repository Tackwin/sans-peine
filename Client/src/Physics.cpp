#include "Physics.hpp"
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
