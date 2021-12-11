
#include "imgui.h"
#include "imgui_ext.h"
#include "imgui-SFML.h"

#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/Graphics/Image.hpp>
#include <SFML/Graphics/Texture.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include <SFML/Graphics/CircleShape.hpp>

#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <optional>
#include <cstdio>

#include <Windows.h>

constexpr auto PI = 3.1415926;
constexpr const char* App_Name = "SP Sim";
constexpr const char* Mail_Name = "\\\\.\\Mailslot\\SP";

constexpr size_t N_Beacons = 3;
struct Reading {
	double beacons[N_Beacons];
	bool pressed;
};

struct State {
	size_t n_beacons_placed = 0;
	std::array<sf::Vector2f, N_Beacons> beacons;

	float x = 0;
	float y = 0;

	double earth_magnetic_field = 0.000'050;
	double beacon_distance = 0.1; // in meters
	double magnet_strength = 1;
	double noise = 0.1f; // in micro tesla

	double max_range = 0.001; // in tesla
	size_t frequency = 50;

	size_t resolution = 16;

	double time_since_write = 0;

	bool recording = false;
	bool right_clicked = false;

	bool fullscreen = false;

	sf::View camera;
	sf::RenderWindow* window = nullptr;
	sf::RenderTarget* renderTarget = nullptr;

	std::vector<Reading> readings;

	HANDLE mail_slot = nullptr;

	struct Space_Simulation {
		size_t resolution = 1000;

		double pivot_distance = 0.01;
		double pivot_angle = PI / 2;
		double R_in = 0.003;
		double R_out = 0.01;
		double h = 0.1;

		double magnet_strength = 1;
	} space_sim;

	struct Simulation_Result {
		double noise = 0;
		size_t resolution = 0;
		std::vector<double> field;

		bool display_x = false;
		bool display_y = false;
		bool display_z = false;

		struct Reading {
			double x = 0;
			double y = 0;
			double z = 0;
		};
		std::vector<Reading> readings;

		struct Candidate_Point {
			size_t x = 0;
			size_t y = 0;
		};
		std::vector<Candidate_Point> candidate_points;
	};

	sf::Texture field_texture;
	std::optional<Simulation_Result> space_res;
};

// Calculate the magnetic field strength 
double B(double M, sf::Vector2f r) noexcept;

void open_slot(State& state) noexcept;
void write_mail(State& state, Reading msg) noexcept;
void update(double dt, State& state) noexcept;
void render(State& state) noexcept;
void render_plot(const std::vector<Reading>& readings) noexcept;

void toggle_fullscreen(State& state) noexcept;

State::Simulation_Result space_sim(State::Space_Simulation& state) noexcept;

#undef main
// Main code
int main(int, char**) {
	State state;
	open_slot(state);
	sf::RenderWindow window(sf::VideoMode(1280, 720), App_Name);
	ImGui::SFML::Init(window);

	state.renderTarget = &window;
	state.window = &window;

	sf::Clock deltaClock;
	while (window.isOpen()) {
		double dt = deltaClock.getElapsedTime().asSeconds();

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
				if (event.key.code == sf::Keyboard::R) state.readings.clear();
			}
		}

		ImGui::SFML::Update(window, deltaClock.restart());
		update(dt, state);

		window.clear();
		render(state);
		ImGui::SFML::Render(window);
		window.display();
	}

	ImGui::SFML::Shutdown();
	CloseHandle(state.mail_slot);
	
	return 0;
}

void update(double dt, State& state) noexcept {
	state.time_since_write += dt;

	state.x = (float)sf::Mouse::getPosition(*state.window).x;
	state.y = (float)sf::Mouse::getPosition(*state.window).y;

	if (state.right_clicked) {
		if (state.n_beacons_placed >= N_Beacons) {
			state.n_beacons_placed = 0;
		} else {
			state.beacons[state.n_beacons_placed] = {
				state.x, state.y
			};
			state.n_beacons_placed++;
		}
	}

	if (state.n_beacons_placed == N_Beacons && state.recording) {
		double max_res = (1 << (state.resolution - 1));
		double dt_beacons =
			(state.beacons[0].x - state.beacons[1].x) * (state.beacons[0].x - state.beacons[1].x) +
			(state.beacons[0].y - state.beacons[1].y) * (state.beacons[0].y - state.beacons[1].y);
		dt_beacons = std::sqrt(dt_beacons);
		static std::mt19937 gen;
		std::uniform_real_distribution<double> dist(
			-state.noise / 1'000'000, state.noise / 1'000'000
		);

		double r1 = state.earth_magnetic_field + B(
			state.magnet_strength,
			sf::Vector2f{
				(float)(state.beacon_distance * (state.x - state.beacons[0].x) / dt_beacons),
				(float)(state.beacon_distance * (state.y - state.beacons[0].y) / dt_beacons)
			}
		);
		r1 = std::min(r1, state.max_range);
		r1 += dist(gen);
		r1 = state.max_range * std::round((r1 / state.max_range) * max_res) / max_res;

		// generate a random number along a normal distribution

		double r2 = state.earth_magnetic_field + B(
			state.magnet_strength,
			sf::Vector2f{
				(float)(state.beacon_distance * (state.x - state.beacons[1].x) / dt_beacons),
				(float)(state.beacon_distance * (state.y - state.beacons[1].y) / dt_beacons)
			}
		);
		r2 = std::min(r2, state.max_range);
		r2 += dist(gen);
		r2 = state.max_range * std::round((r2 / state.max_range) * max_res) / max_res;

		Reading new_reading;
		new_reading.beacons[0] = r1;
		new_reading.beacons[1] = r2;
		new_reading.pressed = sf::Mouse::isButtonPressed(sf::Mouse::Left);

		state.readings.push_back(new_reading);

		if (state.time_since_write >= 1.0 / state.frequency) {
			state.time_since_write = 0;

			if (!state.readings.empty()) write_mail(state, state.readings.back());
		}
	}


}

void upload_field_texture(State& state) noexcept {
	double* gridx = state.space_res->field.data();
	double* gridy = gridx + state.space_res->field.size() / 3;
	double* gridz = gridy + state.space_res->field.size() / 3;

	static sf::Image field_color;
	thread_local bool field_color_init = false;
	if (!field_color_init)
		field_color.create(state.space_res->resolution + 1, state.space_res->resolution + 1);
	field_color_init = true;

	auto color = [](double t) -> auto {
		t = t > 0 ? t : -t;

		t = 1 - std::exp(-t);
		t = std::sqrt(t);

		return (uint8_t)(t * 255);
	};

	double noise = state.space_res->noise * state.space_res->noise;
	for (size_t x = 0; x < state.space_res->field.size() / 3; ++x) {
		double read_x = gridx[x];
		double read_y = gridy[x];
		double read_z = gridz[x];

		sf::Color computed_color;
		auto n_readings = state.space_res->readings.size();
		for (auto& r : state.space_res->readings) {
			if ((r.x - read_x) * (r.x - read_x) < noise) computed_color.r += 255 / n_readings;
			if ((r.y - read_y) * (r.y - read_y) < noise) computed_color.g += 255 / n_readings;
			if ((r.z - read_z) * (r.z - read_z) < noise) computed_color.b += 255 / n_readings;
		}

		computed_color.r = computed_color.r > 0 ? computed_color.r : color(read_x);
		computed_color.g = computed_color.g > 0 ? computed_color.g : color(read_y);
		computed_color.b = computed_color.b > 0 ? computed_color.b : color(read_z);

		field_color.setPixel(
			x % (state.space_res->resolution + 1),
			x / (state.space_res->resolution + 1),
			{
				(uint8_t)(state.space_res->display_x ? computed_color.r : 0),
				(uint8_t)(state.space_res->display_y ? computed_color.g : 0),
				(uint8_t)(state.space_res->display_z ? computed_color.b : 0)
			}
		);
	}

	if (state.field_texture.getSize().x != field_color.getSize().x)
		state.field_texture.create(field_color.getSize().x, field_color.getSize().y);
	state.field_texture.update(field_color);
}

void update_candidate_point(State::Simulation_Result& state) noexcept {
	state.candidate_points.clear();
	if (state.readings.empty()) return;

	double* gridx = state.field.data();
	double* gridy = gridx + state.field.size() / 3;
	double* gridz = gridy + state.field.size() / 3;

	size_t stride = state.resolution + 1;

	thread_local bool* mask = nullptr;
	thread_local size_t mask_size = 0;
	if (mask_size != state.field.size()) {
		if (mask) free(mask);
		mask = (bool*)malloc(state.field.size() / 3);
		mask_size = state.field.size() / 3;
	}

	double noise = state.noise * state.noise;
	for (size_t i = 0; i < mask_size; i++) {
		double dx = gridx[i] - state.readings[0].x;
		double dy = gridy[i] - state.readings[0].y;
		double dz = gridz[i] - state.readings[0].z;

		mask[i] = (dx * dx < noise && dy * dy < noise && dz * dz < noise);
	}

	for (size_t i = 0; i < mask_size; i++) {
		if (mask[i]) {
			size_t sx = 0;
			size_t sy = 0;
			size_t n = 0;

			auto rec = [&](size_t x, size_t y, auto rec) -> void {
				if (mask[x + y * stride]) {
					sx += x;
					sy += y;
					n += 1;
					mask[x + y * stride] = false;

					if (x > 0) rec(x - 1, y, rec);
					if (y > 0) rec(x, y - 1, rec);
					if (x + 1 < stride) rec(x + 1, y, rec);
					if (y + 1 < stride) rec(x, y + 1, rec);
				}
			};

			rec(i % stride, i / stride, rec);
			state.candidate_points.push_back({ sx / n, sy / n });
		}
	}
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
	ImGui::Checkbox("Demo", &open_demo);
	ImGui::SliderDouble("Earth Magnetic Field", &state.earth_magnetic_field, 0, 0.001);
	ImGui::SliderDouble("Beacond Distance", &state.beacon_distance, 0, 1);
	ImGui::SliderDouble("Magnet stength", &state.magnet_strength, 0, 1);
	ImGui::SliderDouble("Sensor noise", &state.noise, 0, 1);
	ImGui::SliderSize("Frequency", &state.frequency, 1, 100);
	ImGui::SliderSize("Resolution", &state.resolution, 1, 32);
	ImGui::Checkbox("Record", &state.recording);

	if (ImGui::Button("Reopen Mail")) {
		CloseHandle(state.mail_slot);
		open_slot(state);
	}

	if (!state.readings.empty()) {
		ImGui::Separator();

		render_plot(state.readings);
	}

	ImGui::Separator();

	if (ImGui::CollapsingHeader("Simulation")) {
		ImGui::PushID("Simulation");
		ImGui::SliderSize("Resolution", &state.space_sim.resolution, 1, 2000);
		ImGui::SliderAngled("Pivot angle", &state.space_sim.pivot_angle, 180, 0);
		ImGui::SliderDouble("Pivot distance", &state.space_sim.pivot_distance, 0, 0.1);
		ImGui::SliderDouble("R_in", &state.space_sim.R_in, 0, state.space_sim.R_out);
		ImGui::SliderDouble("R_out", &state.space_sim.R_out, state.space_sim.R_in, 0.05);
		ImGui::SliderDouble("h", &state.space_sim.h, 0, 0.1);
		ImGui::SliderDouble("strength", &state.space_sim.magnet_strength, 0, 0.1);
		if (ImGui::Button("Compute")) {
			state.space_res = space_sim(state.space_sim);
			state.space_res->display_x = true;
			state.space_res->display_y = true;
			state.space_res->display_z = true;

			upload_field_texture(state);
		}
		if (state.space_res) {
			ImGui::Separator();

			bool dirty = false;
			dirty |= ImGui::Checkbox("X", &state.space_res->display_x);
			dirty |= ImGui::Checkbox("Y", &state.space_res->display_y);
			dirty |= ImGui::Checkbox("Z", &state.space_res->display_z);

			ImGui::Separator();

			if (ImGui::Button("Add reading")) state.space_res->readings.push_back({});
			if (state.space_res->readings.size() > 0) {
				ImGui::SameLine();
				if (ImGui::Button("Del reading")) state.space_res->readings.pop_back();
			}

			dirty |= ImGui::SliderDouble("Noise", &state.space_res->noise, 0, 1, "%.7f", 4);

			for (auto& r : state.space_res->readings) {
				ImGui::PushID(&r);
				dirty |= ImGui::SliderDouble("X", &r.x, 0, 1, "%.7f", 4);
				dirty |= ImGui::SliderDouble("Y", &r.y, 0, 1, "%.7f", 4);
				dirty |= ImGui::SliderDouble("Z", &r.z, 0, 1, "%.7f", 4);
				ImGui::PopID();
			}

			if (dirty) {
				upload_field_texture(state);
				update_candidate_point(*state.space_res);
			}

			ImGui::Separator();

			for (auto& x : state.space_res->candidate_points) {
				ImGui::Text("%zu, %zu", x.x, x.y);
			}
		}
		ImGui::PopID();
	}

	if (state.space_res) {
		ImGui::Begin("Result");
		ImGui::Image(state.field_texture);
		ImGui::End();
	}

	ImGui::End();

	if (open_demo) ImGui::ShowDemoWindow();


	sf::CircleShape shape;
	shape.setRadius(10);
	shape.setPosition(state.x - shape.getRadius(), state.y - shape.getRadius());

	state.renderTarget->draw(shape);

	for (size_t i = 0; i < state.n_beacons_placed; ++i) {
		shape.setPosition(
			state.beacons[i].x - shape.getRadius(), state.beacons[i].y - shape.getRadius()
		);
		shape.setFillColor({255, 0, 0, 255});
		state.renderTarget->draw(shape);
	}

	ImGui::End();

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

	// We add the earth magnetic field
	return u0 * M / (4 * PI * d * d * d);
}

void open_slot(State& state) noexcept {
	state.mail_slot = CreateFileA(
		Mail_Name,
		GENERIC_WRITE,
		FILE_SHARE_READ,
		nullptr,
		OPEN_EXISTING,
		FILE_ATTRIBUTE_NORMAL,
		nullptr
	);
	if (state.mail_slot == INVALID_HANDLE_VALUE) {
		fprintf(stderr, "CreateFile failed with %d.\n", (int)GetLastError());
		return;
	}
}

void write_mail(State& state, Reading msg) noexcept {
	DWORD cbWritten;

	auto fResult = WriteFile(
		state.mail_slot,
		reinterpret_cast<const char*>(&msg),
		sizeof(msg),
		&cbWritten,
		nullptr
	);

	if (!fResult) 
	
	fprintf(stderr, "Write failed with %d.\n", (int)GetLastError());
}

void toggle_fullscreen(State& state) noexcept {
	if (state.fullscreen) state.window->create(
		sf::VideoMode(1280, 720), App_Name
	);
	else state.window->create(
		sf::VideoMode::getFullscreenModes()[0], App_Name, sf::Style::Fullscreen
	);
	state.fullscreen = !state.fullscreen;
}

State::Simulation_Result space_sim(State::Space_Simulation& state) noexcept {
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


	State::Simulation_Result result;

	result.field.resize((state.resolution + 1) * (state.resolution + 1) * 3);

	double* gridx = result.field.data();
	double* gridy = gridx + (state.resolution + 1) * (state.resolution + 1);
	double* gridz = gridy + (state.resolution + 1) * (state.resolution + 1);

	double max_ = 0;
	double min_ = 0;

	for (size_t x = 0; x <= state.resolution; x++) for (size_t y = 0; y <= state.resolution; y++) {
		// the point to probe in cartesian space
		double px = (double)x / state.resolution;
		double py = (double)y / state.resolution - 0.5;
		double pz = 0;
		px *= 0.500; // in meter
		py *= 0.500; // in meter
		pz *= 0.500; // in meter

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
