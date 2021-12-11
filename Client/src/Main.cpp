// dear imgui: standalone example application for SDL2 + OpenGL
// If you are new to dear imgui, see examples/README.txt and documentation at the top of imgui.cpp.
// (SDL is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan/Metal graphics context creation, etc.)
// (GL3W is a helper library to access OpenGL functions since there is no standard header to access modern OpenGL functions easily. Alternatives are GLEW, Glad, etc.)
#include "imgui.h"
#include "imgui-SFML.h"

#include <array>
#include <vector>
#include <cmath>
#include <random>

#include <strsafe.h>

#include "Renderer.hpp"
#include "GUI.hpp"
#include "Physics.hpp"

void toggle_fullscreen(State& state) noexcept;

void make_slot(State& state) noexcept;
std::optional<Reading> read_mail(State& state) noexcept;
void update(State& state) noexcept;
void render(State& state) noexcept;
void upload_field_texture(State& state, Vector3d readings) noexcept;
#undef main

// Main code
int main(int, char**) {
	State state;
	make_slot(state);

	state.context_settings.antialiasingLevel = 8;
	sf::RenderWindow window(
		sf::VideoMode(1280, 720), App_Name, sf::Style::Default, state.context_settings
	);
	window.setFramerateLimit(60);
	ImGui::SFML::Init(window);

	state.renderTarget = &window;
	state.window = &window;

	sf::Clock deltaClock;

	sf::View meter_view;
	meter_view.setCenter(state.camera_pos);
	state.renderTarget->setView(meter_view);

	sf::Clock delta_clock;
	std::optional<sf::Vector2i> last_mouse_pos;
	while (window.isOpen()) {
		auto dt = delta_clock.restart().asSeconds();
		sf::Event event;
		state.right_clicked = false;

		auto r = window.getSize().x / (float)window.getSize().y;
		meter_view.setCenter(state.camera_pos);
		meter_view.setSize({ state.zoom_level, -state.zoom_level / r });
		state.renderTarget->setView(meter_view);

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::D)) state.camera_pos.x += 0.02f * dt;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Q)) state.camera_pos.x -= 0.02f * dt;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Z)) state.camera_pos.y += 0.02f * dt;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::S)) state.camera_pos.y -= 0.02f * dt;

		if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Middle) && last_mouse_pos) {
			auto last = state.renderTarget->mapPixelToCoords(*last_mouse_pos);
			auto curr = state.renderTarget->mapPixelToCoords(sf::Mouse::getPosition(*state.window));
			auto dt_move = sf::Vector2f{ curr.x - last.x, curr.y - last.y };
			
			state.camera_pos.x -= dt_move.x;
			state.camera_pos.y -= dt_move.y;
		}
		last_mouse_pos = sf::Mouse::getPosition(*state.window);

		while (window.pollEvent(event)) {
			ImGui::SFML::ProcessEvent(event);
			if (ImGui::GetIO().WantCaptureMouse) continue;

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
			if (event.type == sf::Event::MouseWheelScrolled) {
				auto dt_wheel = event.mouseWheelScroll.delta;
				if (dt_wheel > 0) state.zoom_level /= dt_wheel + 1;
				if (dt_wheel < 0) state.zoom_level *= 1 - dt_wheel;
			}
		}

		ImGui::SFML::Update(window, deltaClock.restart());
		update(state);

		window.clear();
		render(state);
		render(state, state.gui);

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
			auto pos = sf::Mouse::getPosition(*state.window);
			state.beacons[state.n_beacons_placed].pos =
				state.renderTarget->mapPixelToCoords(pos);
			state.beacons[state.n_beacons_placed].calibration_sample = 0;
			state.beacons[state.n_beacons_placed].sum_sample = {};
			state.n_beacons_placed++;
		}
	}

	auto res = read_mail(state);
	while(res) {

		if (state.gui.calibrating) {
			for (size_t i = 0; i < N_Beacons; ++i) {
				state.beacons[i].calibration_sample++;
				state.beacons[i].sum_sample += res->beacons[i];
			}
		} else {
			avg_readings.push_back(*res);
			if (avg_readings.size() >= state.gui.oversampling) {
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

	if (state.gui.want_compute) {
		state.gui.want_compute = false;

		state.space_res = space_sim(state.gui.space_sim);

		upload_field_texture(state, {0.000065, -0.001500, 0.000783 - 0.000043});
		upload_field_texture(state, {0.000065, -0.000196, 0.000176 - 0.000043});
		// upload_field_texture(state, {0.000065, -0.000055, 0.000176 - 0.000043});
	}

	if (state.gui.want_next_reading) {
		state.gui.want_next_reading = false;

		upload_field_texture(state, state.readings[state.next_reading].beacons[0]);

		state.next_reading++;
		if (!state.readings.empty()) state.next_reading %= state.readings.size();
	}

	state.gui.sample_to_display = std::min(state.gui.sample_to_display, state.readings.size());
}

void render(State& state) noexcept {
	for (float x = -state.zoom_level; x <= state.zoom_level; x += state.gui.grid) {
		float xx = std::round(x / state.gui.grid) * state.gui.grid;
		
		sf::Vertex lines[2] = {
			sf::Vertex({ xx, -state.zoom_level }, state.gui.grid_color),
			sf::Vertex({ xx, +state.zoom_level }, state.gui.grid_color)
		};

		state.renderTarget->draw(lines, 2, sf::Lines);
	}
	for (float y = -state.zoom_level; y <= state.zoom_level; y += state.gui.grid) {
		float yy = std::round(y / state.gui.grid) * state.gui.grid;
		sf::Vertex lines[2] = {
			sf::Vertex({ -state.zoom_level, yy }, state.gui.grid_color),
			sf::Vertex({ +state.zoom_level, yy }, state.gui.grid_color)
		};

		state.renderTarget->draw(lines, 2, sf::Lines);
	}

	sf::CircleShape shape;
	shape.setRadius(0.005);

	for (size_t i = 0; i < state.n_beacons_placed; ++i) {
		shape.setPosition(
			state.beacons[i].pos.x - shape.getRadius(), state.beacons[i].pos.y - shape.getRadius()
		);
		shape.setFillColor({255, 0, 0, 255});
		state.renderTarget->draw(shape);
	}


	render_triangulation(state);
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

void upload_field_texture(State& state, Vector3d r) noexcept {
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

	double noise = 0.000010 * 0.000010;
	auto real = { -0.001440, -0.000560, -0.000261, -0.000131, -0.000065, -0.000035 };
	for (size_t x = 0; x < state.space_res->field.size() / 3; ++x) {
		double read_x = gridx[x];
		double read_y = gridy[x];
		double read_z = gridz[x];

		sf::Color computed_color;

		for (auto& x : real) if ((x - read_y) * (x - read_y) < noise) computed_color.g = 255;

		computed_color.r = computed_color.r > 0 ? computed_color.r : color(read_x);
		computed_color.g = computed_color.g > 0 ? computed_color.g : color(read_y);
		computed_color.b = computed_color.b > 0 ? computed_color.b : color(read_z);

		field_color.setPixel(
			x % (state.space_res->resolution + 1),
			x / (state.space_res->resolution + 1),
			{
				(uint8_t)(state.gui.display_x ? computed_color.r : 0),
				(uint8_t)(state.gui.display_y ? computed_color.g : 0),
				(uint8_t)(state.gui.display_z ? computed_color.b : 0)
			}
		);
	}

	if (state.gui.field_texture.getSize().x != field_color.getSize().x)
		state.gui.field_texture.create(field_color.getSize().x, field_color.getSize().y);
	state.gui.field_texture.update(field_color);
}
