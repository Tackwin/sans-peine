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
#include "Macro.hpp"

#include "IPC.hpp"

void toggle_fullscreen(State& state) noexcept;

void load_drivers_into(std::vector<Driver_Interface>& drivers) noexcept;

std::optional<Reading> read_mail(State& state) noexcept;
void update(State& state) noexcept;
void render(State& state) noexcept;
void upload_field_texture(State& state, Vector3d readings) noexcept;

Debug_Values perf_values;

#undef main

// Main code
int main(int, char**) {
	State state;
	state.mail_slot = IPC::open_slot(Mail_Name);
	if (!state.mail_slot) state.mail_slot = IPC::make_slot(Mail_Name);
	if (!state.mail_slot) {
		printf("couldn't open mail slot :(\n");
		return -1;
	}

	load_drivers_into(state.loaded_drivers);
	for (auto& d : state.loaded_drivers) {
		d.ptr = d.init();
		state.driver_threads.push_back(std::thread(d.play, d.ptr));
	}
	defer {
		for (auto& d : state.loaded_drivers) {
			d.shut(d.ptr);
			FreeLibrary((HINSTANCE)d.lib);
		}
	};

	state.last_sps_timestamp = seconds();

	state.probability_grid = (double*)malloc(
		state.probability_resolution *
		state.probability_resolution *
		sizeof(double)
	);
	for (size_t i = 0; i < N_Beacons; ++i) state.gui.use_dist_beacon[i] = false;
	for (size_t i = 0; i < N_Beacons; ++i) state.gui.use_angle_beacon[i] = false;

	state.context_settings.antialiasingLevel = 8;
	sf::RenderWindow window(
		sf::VideoMode(1280, 720), App_Name, sf::Style::Default, state.context_settings
	);
	window.setFramerateLimit(0);
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
		render(frame_debug_values); frame_debug_values.reset();
		render(perf_values);

		ImGui::SFML::Render(window);
		window.display();
	}

	ImGui::SFML::Shutdown();

	free(state.probability_grid);
	
	return 0;
}

void update(State& state) noexcept {
	thread_local std::vector<Reading> avg_readings;

	if (state.right_clicked) {
		if (state.n_beacons_placed >= N_Beacons) {
			state.n_beacons_placed = 0;
		} else {
			auto pos = sf::Mouse::getPosition(*state.window);
			auto& b = state.beacons[state.n_beacons_placed];
			
			b.pos.x = state.renderTarget->mapPixelToCoords(pos).x;
			b.pos.y = state.renderTarget->mapPixelToCoords(pos).y;
			b.calibration_sample = 0;
			b.sum_sample = {};
			state.n_beacons_placed++;
		}
	}


	auto handle_reading = [&] (Reading res) {
		if (state.gui.calibrating) {
			for (size_t i = 0; i < N_Beacons; ++i) {
				auto& b = state.beacons[i];

				b.calibration_sample++;
				b.sum_sample += res.beacons[i];
				b.sum_dist += std::hypot(res.beacons[i].x, res.beacons[i].y, res.beacons[i].z);
				b.sum2_sample.x += res.beacons[i].x * res.beacons[i].x;
				b.sum2_sample.y += res.beacons[i].y * res.beacons[i].y;
				b.sum2_sample.z += res.beacons[i].z * res.beacons[i].z;

				b.sum2_dist +=
					std::hypot(res.beacons[i].x, res.beacons[i].y, res.beacons[i].z) *
					std::hypot(res.beacons[i].x, res.beacons[i].y, res.beacons[i].z);

				b.mean.x = b.sum_sample.x / b.calibration_sample;
				b.mean.y = b.sum_sample.y / b.calibration_sample;
				b.mean.z = b.sum_sample.z / b.calibration_sample;

				auto n = b.calibration_sample;
				b.std.x = b.sum2_sample.x / n - (b.sum_sample.x * b.sum_sample.x) / (n * n);
				b.std.y = b.sum2_sample.y / n - (b.sum_sample.y * b.sum_sample.y) / (n * n);
				b.std.z = b.sum2_sample.z / n - (b.sum_sample.z * b.sum_sample.z) / (n * n);

				auto f = n / (n - 1.0);
				b.std.x = std::sqrt(b.std.x * f);
				b.std.y = std::sqrt(b.std.y * f);
				b.std.z = std::sqrt(b.std.z * f);
			}
			state.curr_sps_counter += 1;
		} else {
			avg_readings.push_back(res);
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

				state.new_readings.push_back(avg);

				state.curr_sps_counter += 1;
			}
		}
	};
	constexpr size_t Max_New_Per_Frame = 10;
	size_t i_new_per_frame = 0;

	auto res = read_mail(state);
	while(res) {
		handle_reading(*res);
		res = read_mail(state);
		if (++i_new_per_frame > Max_New_Per_Frame) break;
	}

	for (auto& d : state.loaded_drivers) if (d.ptr && d.next) {
		Reading read;
		while (d.next(d.ptr, &read)) {
			handle_reading(read);
			if (++i_new_per_frame > Max_New_Per_Frame) break;
		}
	}

	if (sf::Keyboard::isKeyPressed(sf::Keyboard::R)) state.readings.clear();

	if (state.gui.want_compute) {
		state.gui.want_compute = false;

		state.space_res = space_sim(state.gui.space_sim);

		// upload_field_texture(state, {0.000065, -0.001500, 0.000783 - 0.000043});
		// upload_field_texture(state, {0.000065, -0.000196, 0.000176 - 0.000043});
		// upload_field_texture(state, {0.000065, -0.000055, 0.000176 - 0.000043});
	}

	if (state.gui.want_next_reading) {
		state.gui.want_next_reading = false;

		// upload_field_texture(state, state.readings[state.next_reading].beacons[0]);

		state.next_reading++;
		if (!state.readings.empty()) state.next_reading %= state.readings.size();
	}

	if (state.gui.want_reset_driver) {
		state.gui.want_reset_driver = false;

		for (size_t i = 0; i < state.loaded_drivers.size(); ++i) {
			auto& d = state.loaded_drivers[i];
			d.shut(d.ptr);
			d.ptr = d.init();
			state.driver_threads[i].join();
			state.driver_threads[i] = std::thread(d.play, d.ptr);
		}
	}

	for (auto new_reading : state.new_readings) {
		Simulation_Parameters sim_params;
		sim_params.beacons = state.beacons;
		sim_params.use_dist = state.gui.use_dist_beacon;
		sim_params.use_angle = state.gui.use_angle_beacon;
		sim_params.magnet_strength = state.gui.magnet_strength;
		sim_params.h = state.gui.magnet_height;
		sim_params.reading = new_reading;
		sim_params.sensitivity = state.gui.sensitivity;

		auto start = seconds();
		auto t1 = seconds();
		auto sim_res = space_sim(sim_params);
		compute_probability_grid(state, sim_res);

		auto w = state.probability_resolution;
		auto h = state.probability_resolution;

		auto max_trace = [&] {
			size_t max_idx = 0;
			for (size_t i = 0; i < w * h; ++i)
				if (state.probability_grid[i] > state.probability_grid[max_idx]) max_idx = i;

			state.estimated_points.push_back({
				state.probability_space_size * ((max_idx % w) / (w - 1.0) - 0.5),
				state.probability_space_size * ((max_idx / h) / (h - 1.0) - 0.5)
			});
		};
		auto avg_trace = [&] {
			Vector2d sum = {};
			long double norm = 0;
			for (size_t i = 0; i < w * h; ++i) norm += state.probability_grid[i];

			for (size_t xi = 0; xi < w; ++xi) for (size_t yi = 0; yi < h; ++yi) {
				auto px = xi / (w - 1.0) - 0.5;
				auto py = yi / (h - 1.0) - 0.5;

				px *= state.probability_space_size * 1;
				py *= state.probability_space_size * 1;

				auto p = state.probability_grid[xi + yi * w] / norm;
				sum.x += p * px;
				sum.y += p * py;
			}

			state.estimated_points.push_back(sum);
		};
		auto avg2_trace = [&] {
			Vector2d sum = {};
			long double norm = 0;

			for (size_t xi = 0; xi < w; ++xi) for (size_t yi = 0; yi < h; ++yi) {
				auto px = xi / (w - 1.0) - 0.5;
				auto py = yi / (h - 1.0) - 0.5;

				px *= state.probability_space_size * 1;
				py *= state.probability_space_size * 1;

				auto p = state.probability_grid[xi + yi * w];
				p *= p;
				norm += p;
				sum.x += p * px;
				sum.y += p * py;
			}

			sum.x /= norm;
			sum.y /= norm;

			state.estimated_points.push_back(sum);
		};

		max_trace();
		avg_trace();
		avg2_trace();
		auto elapsed = seconds() - start;
		perf_values.add_to_distribution("All", elapsed);
		state.readings.push_back(new_reading);
	}
	state.new_readings.clear();

	if (seconds() - state.last_sps_timestamp > 1) {
		state.last_sps_counter = state.curr_sps_counter;
		state.curr_sps_counter = 0;
		state.last_sps_timestamp = seconds();
	}
}

void render(State& state) noexcept {
	thread_local double C_avg = 0;
	thread_local size_t C_n = 0;
	auto t = seconds();
	update_probability_texture(state);
	C_avg += seconds() - t;
	C_n++;
	frame_debug_values.watch("    Upload time", seconds() - t);
	frame_debug_values.watch("Avg Upload time", C_avg / C_n);

	sf::Sprite probability_sprite;
	probability_sprite.setTexture(state.probability_texture);
	probability_sprite.setOrigin(state.probability_texture.getSize().x / 2, state.probability_texture.getSize().y / 2);
	probability_sprite.setPosition(0, 0);
	probability_sprite.setScale(
		state.probability_space_size / state.probability_texture.getSize().x,
		state.probability_space_size / state.probability_texture.getSize().y
	);
	state.renderTarget->draw(probability_sprite);

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
	render_estimation_trace(state);
	// render_triangulation(state);
}

std::optional<Reading> read_mail(State& state) noexcept {
	Reading r;

	if (IPC::read(state.mail_slot, &r, sizeof(r))) return r;
	return std::nullopt;
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
#if 0
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
#endif
}

void load_drivers_into(std::vector<Driver_Interface>& drivers) noexcept {
	HINSTANCE dll = LoadLibraryA("Driver.dll");
	if (!dll) return;

	Driver_Interface driver;
	driver.lib = dll;
	driver.init = (driver_init_f)GetProcAddress(dll, "init");
	driver.shut = (driver_shut_f)GetProcAddress(dll, "shut");
	driver.play = (driver_play_f)GetProcAddress(dll, "play");
	driver.stop = (driver_play_f)GetProcAddress(dll, "stop");
	driver.next = (driver_next_f)GetProcAddress(dll, "next");

	drivers.push_back(driver);
}
