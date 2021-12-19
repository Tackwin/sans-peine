#include "Renderer.hpp"

std::optional<sf::Vector2f> triangulate(State& state, Reading r, size_t bi, size_t bj) noexcept;
sf::Vector2f line_line(
	float p0_x, float p0_y, float p1_x, float p1_y, float p2_x, float p2_y, float p3_x, float p3_y
) noexcept;
double mag_to_meter(State& state, double mag, Beacon b) noexcept;
std::optional<sf::Vector2f> circle_circle_intersection(
	sf::Vector2f c0, double r0, sf::Vector2f c1, double r1
) noexcept;

void render_triangulation(State& state) noexcept {
	const std::uint32_t Distinct_Colors[8] = {
		0x0000FFFF,
		0x00FF00FF,
		0xFF0000FF,
		0xFFFF00FF,
		0xFF00FFFF,
		0x00FFFFFF,
		0xFFFFFFFF,
		0xF0F0F0FF
	};

	constexpr auto N_Combi = (N_Beacons * (N_Beacons - 1)) / 2;
	static std::vector<std::array<sf::Vector2f, N_Combi>> lines;

	std::array<sf::Vector2f, N_Combi> default_set;
	for (size_t i = 0; i < N_Beacons; ++i) default_set[i] = { NAN, NAN };

	lines.clear();
	lines.resize(state.readings.size(), default_set);

	size_t bn = 0;
	auto to_draw = state.gui.sample_live ? state.readings.size() : state.gui.sample_to_display;
	for (size_t bi = 0; bi < N_Beacons; ++bi) for (size_t bj = bi + 1; bj < N_Beacons; ++bj, ++bn) {
		std::optional<sf::Vector2f> last;
		if (state.readings.size() > 0) last = triangulate(state, state.readings.front(), bi, bj);

		if (last) lines[0][bn] = *last;

		for (size_t i = 1; i < to_draw; ++i) {
			auto next = triangulate(state, state.readings[i], bi, bj);

			if (!next) continue;
			lines[i][bn] = *next;

			if (!last) { last = next; continue; }
			if (!state.gui.display_matrix[bi + bj * N_Beacons]) continue;

			auto color = sf::Color(Distinct_Colors[bn % 8]);
			sf::Vertex line[] = { sf::Vertex(*last, color), sf::Vertex(*next, color)};
			state.renderTarget->draw(line, 2, sf::Lines);

			last = next;
		}
	}

	for (size_t i = 0; i < N_Combi; ++i) {
		for (size_t x = 1; x < lines.size(); ++x) if (lines[x][i].x == NAN) {
			sf::Vector2f start = lines[x][i];

			size_t n = 0;
			for (; x + n < lines.size() && lines[x + n][i].x == NAN; ++n);

			if (x + n == lines.size()) {
				for (; x < lines.size(); ++x) lines[x][i] = start;
				continue;
			}

			sf::Vector2f end = lines[x + n][i];

			for (size_t j = 0; j < n; ++j) {
				lines[x + j][i].x = start.x + n * (end.x - start.x) / (j + 1.f);
				lines[x + j][i].y = start.y + n * (end.y - start.y) / (j + 1.f);
			}
		}
	}

	for (size_t i = 0; i + 1 < lines.size(); ++i) {
		sf::Vector2f curr = {};
		for (auto& x : lines[i]) curr += x;
		curr.x /= N_Combi;
		curr.y /= N_Combi;

		sf::Vector2f next = {};
		for (auto& x : lines[i + 1]) next += x;
		next.x /= N_Combi;
		next.y /= N_Combi;

		sf::Vertex l[] = { sf::Vertex(curr), sf::Vertex(next) };
		state.renderTarget->draw(l, 2, sf::Lines);
	}

	sf::CircleShape shape;
	shape.setRadius(0.005);

	if (state.readings.size() > 0 && to_draw > 0) {
		for (size_t i = 0; i < state.n_beacons_placed; ++i) {

			auto r = state.readings[to_draw - 1];
			auto m = std::hypot(r.beacons[i].x, r.beacons[i].y, r.beacons[i].z);
			shape.setRadius(mag_to_meter(state, m, state.beacons[i]));
			shape.setPosition(
				state.beacons[i].pos.x - shape.getRadius(),
				state.beacons[i].pos.y - shape.getRadius()
			);
			shape.setFillColor({0, 0, 0, 0});
			shape.setOutlineColor({255, 255, 255, 255});
			shape.setOutlineThickness(0.0001);

			state.renderTarget->draw(shape);
		}
	}

	if (state.readings.size() > 0 && to_draw > 0) {
		for (size_t i = 0; i < state.n_beacons_placed; ++i) {

			auto r = state.readings[to_draw - 1];

			double base_ix = state.beacons[i].sum_sample.x / state.beacons[i].calibration_sample;
			double ix = r.beacons[i].x - base_ix;
			double base_iy = state.beacons[i].sum_sample.y / state.beacons[i].calibration_sample;
			double iy = r.beacons[i].y - base_iy;

			auto alpha = (float)std::atan2(ix, -iy);


			sf::Vertex line[2] = {
				sf::Vertex(state.beacons[i].pos),
				sf::Vertex({
					state.beacons[i].pos.x + std::cosf(alpha),
					state.beacons[i].pos.y + std::sinf(alpha)
				})
			};

			state.renderTarget->draw(line, 2, sf::Lines);
		}
	}


}

std::optional<sf::Vector2f> triangulate(State& state, Reading r, size_t i, size_t j) noexcept {
	constexpr double u0 = 1.225663753e-6;
	constexpr double PI = 3.141592653589793238462643383279502884;

	if (state.n_beacons_placed < N_Beacons) return {};

	sf::Vector2f b1 = state.beacons[i].pos;
	sf::Vector2f b2 = state.beacons[j].pos;

#if 0
	double bi = std::hypot(r.beacons[i].x, r.beacons[i].y, r.beacons[i].z);
	double bj = std::hypot(r.beacons[j].x, r.beacons[j].y, r.beacons[j].z);

	double d1 = mag_to_meter(state, bi, state.beacons[i]);
	double d2 = mag_to_meter(state, bj, state.beacons[j]);

	auto res = circle_circle_intersection(b1, d1, b2, d2);
	if (!res) return {};
	return *res;
#else
	double ix =
		r.beacons[i].x - state.beacons[i].sum_sample.x / state.beacons[i].calibration_sample;
	double iy =
		r.beacons[i].y - state.beacons[i].sum_sample.y / state.beacons[i].calibration_sample;
		
	double jx =
		r.beacons[j].x - state.beacons[j].sum_sample.x / state.beacons[j].calibration_sample;
	double jy =
		r.beacons[j].y - state.beacons[j].sum_sample.y / state.beacons[j].calibration_sample;

	auto base_length = std::hypot((b1 - b2).x, (b1 - b2).y);


	// The magnetometer have the folowing axis
/*
Y
^
|         x
|         ^
|         |
|    y<---+
|
|
|         x
|         ^
|         |
|    y<---+
|
+---------------> X


*/

	// So we need to do some swapping here
	double alpha = std::atan2(ix, -iy);
	double beta  = std::atan2(jx, -jy);

	sf::Vector2f b1_unit = b1;
	sf::Vector2f b2_unit = b2;
	b1_unit.x += (float)std::cos(alpha);
	b1_unit.y += (float)std::sin(alpha);
	b2_unit.x += (float)std::cos(beta);
	b2_unit.y += (float)std::sin(beta);

	return line_line(b1.x, b1.y, b1_unit.x, b1_unit.y, b2.x, b2.y, b2_unit.x, b2_unit.y);
#endif
}

double mag_to_meter(State& state, double mag, Beacon b) noexcept {
	constexpr double u0 = 1.225663753e-6;
	constexpr double PI = 3.141592653589793238462643383279502884;
	double base_strength = std::hypot(
		b.sum_sample.x / b.calibration_sample,
		b.sum_sample.y / b.calibration_sample,
		b.sum_sample.z / b.calibration_sample
	);
	double strength = std::abs(mag - base_strength);
	
	auto res = state.gui.magnet_strength * u0 / (4 * PI * strength);
	return std::pow(res, 1 / state.gui.power);
}


std::optional<sf::Vector2f> circle_circle_intersection(
	sf::Vector2f c0, double r0,
	sf::Vector2f c1, double r1
) noexcept {
    double a, dx, dy, d, h, rx, ry;
    double x2, y2;

    /* dx and dy are the vertical and horizontal distances between
     * the circle centers.
     */
    dx = c1.x - c0.x;
    dy = c1.y - c0.y;

    /* Determine the straight-line distance between the centers. */
    //d = sqrt((dy*dy) + (dx*dx));
    d = hypot(dx,dy); // Suggested by Keith Briggs

    /* Check for solvability. */
    if (d > (r0 + r1)) {
        /* no solution. circles do not intersect. */
        return std::nullopt;
    }
    if (d < fabs(r0 - r1)) {
        /* no solution. one circle is contained in the other */
        return std::nullopt;
    }

    /* 'point 2' is the point where the line through the circle
     * intersection points crosses the line between the circle
     * centers.  
     */

    /* Determine the distance from point 0 to point 2. */
    a = ((r0*r0) - (r1*r1) + (d*d)) / (2.0 * d) ;

    /* Determine the coordinates of point 2. */
    x2 = c0.x + (dx * a/d);
    y2 = c0.y + (dy * a/d);

    /* Determine the distance from point 2 to either of the
     * intersection points.
     */
    h = sqrt((r0*r0) - (a*a));

    /* Now determine the offsets of the intersection points from
     * point 2.
     */
    rx = -dy * (h/d);
    ry = dx * (h/d);

    return sf::Vector2f{
    	(float)(x2 - rx), (float)(y2 + ry)
    };
}

sf::Vector2f line_line(
	float p0_x, float p0_y, float p1_x, float p1_y, float p2_x, float p2_y, float p3_x, float p3_y
) noexcept {
	float s1_x, s1_y, s2_x, s2_y;
	s1_x = p1_x - p0_x;
	s1_y = p1_y - p0_y;
	s2_x = p3_x - p2_x;
	s2_y = p3_y - p2_y;

	float s, t;
	s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
	t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

	return { p0_x + (t * s1_x), p0_y + (t * s1_y) };
}

void update_probability_texture(State& state) noexcept {
	static sf::Image probability_image;
	auto w = (size_t)(state.probability_space_size / state.probability_resolution);
	auto h = (size_t)(state.probability_space_size / state.probability_resolution);
	probability_image.create(w, h);

	double log_max = -INFINITY;
	double log_min = +INFINITY;

	auto scale = [] (double x) -> double {
		return std::pow(x, 1/3.0);
	};
	auto cmap = [] (double x) -> sf::Color {
		auto c = (sf::Uint8)(x * 255);
		auto [r, g, b] = Viridis_Color_Map[c];
		return {
			(sf::Uint8)(r * 255),
			(sf::Uint8)(g * 255),
			(sf::Uint8)(b * 255),
			255
		};
	};

	for (size_t x = 0; x < probability_image.getSize().x; ++x)
	for (size_t y = 0; y < probability_image.getSize().y; ++y) {
		auto it = scale(state.probability_grid[x + y * w]);
		if (log_max < it) log_max = it;
		if (log_min > it) log_min = it;
	}
	for (size_t x = 0; x < probability_image.getSize().x; ++x)
	for (size_t y = 0; y < probability_image.getSize().y; ++y) {
		auto it = scale(state.probability_grid[x + y * w]);
		auto t = (it - log_min) / (log_max - log_min);
		probability_image.setPixel(x, probability_image.getSize().y - y - 1, cmap(t));
	}

	state.probability_texture.loadFromImage(probability_image);
}
