//7 dear imgui: standalone example application for SDL2 + OpenGL
// If you are new to dear imgui, see examples/README.txt and documentation at the top of imgui.cpp.
// (SDL is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan/Metal graphics context creation, etc.)
// (GL3W is a helper library to access OpenGL functions since there is no standard header to access modern OpenGL functions easily. Alternatives are GLEW, Glad, etc.)

#include "imgui.h"
#include "imgui_ext.h"
#include "imgui-SFML.h"

#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include <SFML/Graphics/CircleShape.hpp>

#include <array>
#include <vector>
#include <cmath>
#include <random>

#include <Windows.h>
#include <strsafe.h>

#include <optional>

constexpr const char* Mail_Name = "\\\\.\\Mailslot\\SP";
constexpr const char* App_Name = "SP Client";

constexpr size_t N_Beacons = 3;

#pragma pack(1)
struct Reading {
	double beacons[N_Beacons] = { 0 };
	bool pressed;
};

struct Beacon {
	sf::Vector2f pos;

	sf::Vector3f sum_sample = {};
	size_t calibration_sample = 0;
};

struct State {
	size_t n_beacons_placed = 0;
	std::array<Beacon, N_Beacons> beacons;
	bool right_clicked = false;
	double beacon_distance = 0.1; // in meters

	size_t oversampling = 1;

	double info_distance = 0.3; // in meters

	sf::ContextSettings context_settings;

	sf::RenderTarget* renderTarget = nullptr;
	sf::RenderWindow* window = nullptr;
	double magnet_strength = 1;

	HANDLE mail_slot = nullptr;
	bool fullscreen = false;

	std::vector<Reading> readings;

	bool calibrating = false;
};


void toggle_fullscreen(State& state) noexcept;

void make_slot(State& state) noexcept;
std::optional<Reading> read_mail(State& state) noexcept;
void update(State& state) noexcept;
void render(State& state) noexcept;
void render_plot(const std::vector<Reading>& readings) noexcept;

std::optional<sf::Vector2f> triangulate(State& state, Reading r, size_t bi, size_t bj) noexcept;

#undef main
// Main code
int main(int, char**) {
	State state;
	make_slot(state);

	state.context_settings.antialiasingLevel = 4;
	sf::RenderWindow window(
		sf::VideoMode(1280, 720), App_Name, sf::Style::Default, state.context_settings
	);
	window.setFramerateLimit(60);
	ImGui::SFML::Init(window);

	state.renderTarget = &window;
	state.window = &window;

	sf::Clock deltaClock;
	while (window.isOpen()) {
		sf::Event event;
		state.right_clicked = false;
		while (window.pollEvent(event)) {
			ImGui::SFML::ProcessEvent(event);

			if (event.type == sf::Event::Closed) {
				window.close();
			}
			if (event.type == sf::Event::MouseButtonPressed) {
				if (event.mouseButton.button == sf::Mouse::Right) {
					state.right_clicked = true;
				}
			}
			if (event.type == sf::Event::KeyPressed) {
				if (event.key.code == sf::Keyboard::F11) toggle_fullscreen(state);
			}
		}

		ImGui::SFML::Update(window, deltaClock.restart());
		update(state);

		window.clear();
		render(state);
		ImGui::SFML::Render(window);
		window.display();
	}

	ImGui::SFML::Shutdown();
	
	return 0;
}

void update(State& state) noexcept {
	thread_local std::vector<Reading> avg_readings;

	if (state.right_clicked) {
		if (state.n_beacons_placed >= N_Beacons) {
			state.n_beacons_placed = 0;
		} else {
			state.beacons[state.n_beacons_placed].pos = {
				(float)sf::Mouse::getPosition(*state.window).x,
				(float)sf::Mouse::getPosition(*state.window).y
			};
			state.beacons[state.n_beacons_placed].calibration_sample = 0;
			state.beacons[state.n_beacons_placed].sum_sample = {};
			state.n_beacons_placed++;
		}
	}

	auto res = read_mail(state);
	while(res) {

		if (state.calibrating) {
			for (size_t i = 0; i < N_Beacons; ++i) {
				state.beacons[i].calibration_sample++;
				state.beacons[i].sum_sample.x += res->beacons[i];
			}
		} else {
			avg_readings.push_back(*res);
			if (avg_readings.size() >= state.oversampling) {
				Reading avg = {};
				size_t n_pressed = 0;
				for (auto& x : avg_readings) {
					for (size_t i = 0; i < N_Beacons; ++i) avg.beacons[i] += x.beacons[i];
					n_pressed += (x.pressed ? 1 : 0);
				}

				for (size_t i = 0; i < N_Beacons; ++i) avg.beacons[i] /= avg_readings.size();
				avg.pressed = n_pressed > avg_readings.size() / 2;

				avg_readings.clear();

				state.readings.push_back(avg);
			}
		}
		res = read_mail(state);
	}

	if (sf::Keyboard::isKeyPressed(sf::Keyboard::R)) state.readings.clear();
}

void render(State& state) noexcept {
	static bool open_demo = false;
	static ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_PassthruCentralNode;
	ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking;
	ImGuiViewport* viewport = ImGui::GetMainViewport();
	ImGui::SetNextWindowPos(viewport->GetWorkPos());
	ImGui::SetNextWindowSize(viewport->GetWorkSize());
	ImGui::SetNextWindowViewport(viewport->ID);
	ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
	ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
	window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse;
	window_flags |= ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
	window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;
	window_flags |= ImGuiWindowFlags_NoBackground;

	ImGui::SetNextWindowBgAlpha(0);
	ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
	ImGui::Begin("DockSpace", nullptr, window_flags);
	ImGui::PopStyleVar(3);

	ImGuiID dockspace_id = ImGui::GetID("MyDockSpace");
	ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), dockspace_flags);

	ImGui::Begin("Parameters");

	ImGui::SliderDouble("Info Distance", &state.info_distance, 0, 1);
	ImGui::SliderDouble("Beacond Distance", &state.beacon_distance, 0, 1);
	ImGui::SliderDouble("Magnet Strength", &state.magnet_strength, 0, 1);

	ImGui::SliderSize("Oversampling", &state.oversampling, 1, 100);

	if (!state.readings.empty()) {
		ImGui::Separator();
		render_plot(state.readings);
	}

	if (ImGui::Button("Clear")) {
		state.readings.clear();
		// for (auto& x : state.beacons) {
		// 	x.calibration_sample = 0;
		// 	x.sum_sample = {};
		// }
	}

	ImGui::Checkbox("Calibrating", &state.calibrating);


	ImGui::End();

	ImGui::End();

	sf::CircleShape shape;
	shape.setRadius(10);

	for (size_t i = 0; i < state.n_beacons_placed; ++i) {
		shape.setPosition(
			state.beacons[i].pos.x - shape.getRadius(), state.beacons[i].pos.y - shape.getRadius()
		);
		shape.setFillColor({255, 0, 0, 255});
		state.renderTarget->draw(shape);
	}

	const std::uint32_t Distinct_Colors[8] = {
		0x6969FF,
		0x22FF22,
		0xFF0000,
		0xFFFF00,
		0x483dFF,
		0x00FFFF,
		0x0000FF,
		0xFF00FF
	};

	constexpr auto N_Combi = (N_Beacons * (N_Beacons - 1)) / 2;
	static std::vector<std::array<sf::Vector2f, N_Combi>> lines;

	std::array<sf::Vector2f, N_Combi> default_set;
	for (size_t i = 0; i < N_Beacons; ++i) default_set[i] = { NAN, NAN };

	lines.clear();
	lines.resize(state.readings.size(), default_set);

	size_t bn = 0;
	for (size_t bi = 0; bi < N_Beacons; ++bi) for (size_t bj = bi + 1; bj < N_Beacons; ++bj, ++bn) {
		std::optional<sf::Vector2f> last;
		if (state.readings.size() > 0) last = triangulate(state, state.readings.front(), bi, bj);

		if (last) lines[0][bn] = *last;

		for (size_t i = 1; i < state.readings.size(); ++i) {
			auto next = triangulate(state, state.readings[i], bi, bj);

			if (!next) continue;
			if (!last) { last = next; continue; }

			auto color = sf::Color(Distinct_Colors[bn % 8]);
			sf::Vertex line[] = { sf::Vertex(*last, color), sf::Vertex(*next, color)};
			state.renderTarget->draw(line, 2, sf::Lines);

			last = next;
			lines[i][bn] = *next;
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


	if (state.n_beacons_placed >= N_Beacons) {
		sf::Vector2f b1 = state.beacons[0].pos;
		sf::Vector2f b2 = state.beacons[1].pos;

		double dt_beacons = (b1.x - b2.x) * (b1.x - b2.x) + (b1.y - b2.y) * (b1.y - b2.y);
		dt_beacons = std::sqrt(dt_beacons);
		shape.setRadius(
			(dt_beacons * state.info_distance / state.beacon_distance) / 2
		);
		shape.setPosition(
			(state.beacons[0].pos.x + state.beacons[1].pos.x) / 2 - shape.getRadius(),
			(state.beacons[0].pos.y + state.beacons[1].pos.y) / 2 - shape.getRadius()
		);
		shape.setFillColor({0, 0, 0, 0});
		shape.setOutlineColor({255, 255, 255, 255});
		shape.setOutlineThickness(3);

		state.renderTarget->draw(shape);

	}
}


void render_plot(const std::vector<Reading>& readings) noexcept {
	static std::vector<std::vector<float>> lines;
	static std::vector<const float*> ys;
	static size_t N_Samples = 10000;

	static std::array<char*, N_Beacons> names = { nullptr };
	for (size_t i = 0; i < N_Beacons; ++i) {
		free(names[i]);
		names[i] = (char*)malloc(sizeof("Beacons XXXXXXX"));
		sprintf(names[i], "Beacon %d", (int)i);
	}

	lines.clear();
	ys.clear();

	lines.resize(N_Beacons);

	float maxs[N_Beacons];
	size_t n = readings.size();

	for (size_t j = 0; j < N_Beacons; ++j) {
		for (size_t i = (n > N_Samples) ? (n - N_Samples) : 0; i < n; ++i) {
			auto& a = readings[i].beacons[j];
			lines[j].push_back((float)a);
			maxs[j] = (float)(a > maxs[j] ? a : maxs[j]);
		}
	}

	for (size_t i = 0; i < N_Beacons; ++i) ys.push_back(lines[i].data());

	ImGui::PlotConfig conf;
	conf.values.ys_list = ys.data();
	conf.values.ys_count = N_Beacons;
	conf.values.count = (int)std::min(readings.size(), N_Samples);

	conf.scale.min = 0;
	for (auto& x : maxs) conf.scale.max = std::max(conf.scale.max, x);
	conf.scale.max *= 1.1f;

	conf.tooltip.show = true;

	conf.grid_y.show = true;
	conf.grid_y.size = conf.scale.max / 5;
	conf.grid_y.subticks = 5;

	conf.line_thickness = 2.f;

	conf.frame_size = { ImGui::GetWindowWidth() * 0.95f, 100 };

	conf.tooltip.ys_names = names.data();
	conf.tooltip.format = "%s, %g: % 12.9f";

	ImGui::Plot("Readings", conf);

	for (size_t i = 0; i < N_Beacons; ++i) {
		conf.values.ys_list = ys.data() + i;
		conf.values.ys_count = 1;
		conf.scale.max = 1.1f * maxs[i];

		conf.grid_y.size = conf.scale.max / 5;

		ImGui::Plot(names[i], conf);
	}

	ImGui::SliderSize("#Samples", &N_Samples, 10, 1000);
}


double B(double M, sf::Vector2f r) noexcept {
	constexpr double u0 = 1.225663753e-6;
	constexpr double PI = 3.141592653589793238462643383279502884;

	double d = std::sqrt(r.x * r.x + r.y * r.y);

	return u0 * M / (4 * PI * d * d * d);
}

void make_slot(State& state) noexcept {
	state.mail_slot = CreateMailslotA(Mail_Name, 0, 1, nullptr);

	if (state.mail_slot == INVALID_HANDLE_VALUE) {
		fprintf(stderr, "Error creating mail slot\n");
		std::terminate();
	}
}

std::optional<Reading> read_mail(State& state) noexcept {
	DWORD cbRead = 0; 
	BOOL fResult; 

	Reading result;

	fResult = ReadFile(state.mail_slot, &result, sizeof(result), &cbRead, nullptr);

	if (cbRead < sizeof(Reading)) return std::nullopt;

	if (!fResult) {
		fprintf(stderr, "ReadFile failed with %d.\n", (int)GetLastError());
		return std::nullopt;
	}

	return result;
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
    rx = dy * (h/d);
    ry = dx * (h/d);

    return sf::Vector2f{
    	(float)(x2 + rx), (float)(y2 - ry)
    };
}

std::optional<sf::Vector2f> triangulate(State& state, Reading r, size_t i, size_t j) noexcept {
	constexpr double u0 = 1.225663753e-6;
	constexpr double PI = 3.141592653589793238462643383279502884;

	if (state.n_beacons_placed < 2) return {};

	sf::Vector2f b1 = state.beacons[i].pos;
	sf::Vector2f b2 = state.beacons[j].pos;

	double dt_beacons = (b1.x - b2.x) * (b1.x - b2.x) + (b1.y - b2.y) * (b1.y - b2.y);

	dt_beacons = std::sqrt(dt_beacons);

	double base_strength1 = state.beacons[i].sum_sample.x / state.beacons[i].calibration_sample;
	double strength1 = r.beacons[i] - base_strength1;

	double base_strength2 = state.beacons[j].sum_sample.x / state.beacons[j].calibration_sample;
	double strength2 = r.beacons[j] - base_strength2;

	double d1 = state.magnet_strength * u0 / (4 * PI * strength1);
	if (d1 < 0) return {};
	d1 = std::pow(d1, 1.0 / 3.0);
	d1 = dt_beacons * d1 / state.beacon_distance;

	double d2 = state.magnet_strength * u0 / (4 * PI * strength2);
	if (d2 < 0) return {};
	d2 = std::pow(d2, 1.0 / 3.0);
	d2 = dt_beacons * d2 / state.beacon_distance;

	auto res = circle_circle_intersection(b1, d1, b2, d2);
	if (!res) return {};
	return *res;
}


void toggle_fullscreen(State& state) noexcept {
	if (state.fullscreen) state.window->create(
		sf::VideoMode(1280, 720), App_Name, sf::Style::Default, state.context_settings
	);
	else state.window->create(
		sf::VideoMode::getFullscreenModes()[0],
		App_Name,
		sf::Style::Fullscreen,
		state.context_settings
	);
	state.fullscreen = !state.fullscreen;
}
