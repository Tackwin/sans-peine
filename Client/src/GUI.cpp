#include "GUI.hpp"

#include "Renderer.hpp"

#include "imgui.h"
#include "imgui_ext.h"
#include "imgui-SFML.h"

void render_plot(const std::vector<Reading>& readings) noexcept;

void render(State& state, GUI_State& gui_state) noexcept {
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
	static float temp_col[3] = { 0.1f, 0.1f, 0.1f };
	ImGui::ColorEdit3("Grid color", temp_col);
	gui_state.grid_color = {
		(uint8_t)(temp_col[0] * 255),
		(uint8_t)(temp_col[1] * 255),
		(uint8_t)(temp_col[2] * 255)
	};

	ImGui::Text("Sample per seconds: %zu", state.last_sps_counter);

	ImGui::SliderDouble("Epsilon", &gui_state.epsilon, 0, 100, "%.7f", 6);
	ImGui::SliderDouble("Power", &gui_state.power, 0, 10);
	ImGui::SliderDouble("Magnet Strength", &gui_state.magnet_strength, 0, 100, "%.7f", 6);
	ImGui::SliderDouble("Magnet Height", &gui_state.magnet_height, 0, 0.10);
	ImGui::SliderSize("Oversampling", &gui_state.oversampling, 1, 100);
	ImGui::SliderDouble("Sensitivity", &gui_state.sensitivity, 0, 0.5, "%.6f", 2);

	const char* items[] = { "Max", "Avg", "Avg2" };
	int item_current = (int)gui_state.trace_mode;
	ImGui::ListBox("Trace Mode", &item_current, items, (int)GUI_State::Trace_Mode::Count);
	gui_state.trace_mode = (GUI_State::Trace_Mode)item_current;

	ImGui::Checkbox("Live", &gui_state.sample_live);
	ImGui::SameLine();
	ImGui::SliderSize("# Sample", &gui_state.sample_to_display, 0, state.readings.size());
	ImGui::SameLine();
	if (ImGui::Button("+")) gui_state.sample_to_display += 1;
	ImGui::SameLine();
	if (ImGui::Button("-") && gui_state.sample_to_display > 0) gui_state.sample_to_display -= 1;

	if (ImGui::Button("Default placement")) {
		state.n_beacons_placed = 0;
		for (size_t n = 0; n < N_Beacons; ++n) {
			auto t = 0.5 - n / (N_Beacons - 1.0);

			state.beacons[state.n_beacons_placed] = {};
			state.beacons[state.n_beacons_placed].pos = Beacons_Pos[n];
			state.n_beacons_placed++;
		}
	}

	if (!state.readings.empty()) {
		ImGui::Separator();
		render_plot(state.readings);
	}

	if (ImGui::Button("Clear")) {
		state.readings.clear();
		state.estimated_points.clear();
		// for (auto& x : state.beacons) {
		// 	x.calibration_sample = 0;
		// 	x.sum_sample = {};
		// }
	}

	ImGui::PushID("Checkbox");
	for (size_t i = 0; i < N_Beacons; ++i) {
		for (size_t j = 0; j < N_Beacons; ++j) {
			ImGui::PushID(i + j * N_Beacons);
			ImGui::Checkbox("", gui_state.display_matrix + i + j * N_Beacons);
			ImGui::SameLine();
			ImGui::PopID();
		}
		ImGui::NewLine();
	}
	ImGui::PopID();

	thread_local char temp_buffer[256];
	ImGui::PushID("Selection measures");
	thread_local bool all_dist = false;
	thread_local bool all_angle = false;

	ImGui::PushID("Dist");
	if (ImGui::Checkbox("All Distance", &all_dist))
		for (auto& x : gui_state.use_dist_beacon) x = all_dist;
	for (size_t i = 0; i < N_Beacons; ++i) {
		ImGui::SameLine();
		sprintf(temp_buffer, "%d", (int)i);
		ImGui::Checkbox(temp_buffer, &gui_state.use_dist_beacon[i]);
	}
	ImGui::PopID();

	ImGui::PushID("Angle");
	if (ImGui::Checkbox("All Angle", &all_angle))
		for (auto& x : gui_state.use_angle_beacon) x = all_angle;
	for (size_t i = 0; i < N_Beacons; ++i) {
		ImGui::SameLine();
		sprintf(temp_buffer, "%d", (int)i);
		ImGui::Checkbox(temp_buffer, &gui_state.use_angle_beacon[i]);
	}
	ImGui::PopID();

	ImGui::Checkbox("Compute mag", &gui_state.compute_mag);
	ImGui::Checkbox("Compute acc", &gui_state.compute_acc);

	ImGui::PopID();

	for (size_t i = 0; i < N_Beacons; ++i) for (size_t j = 0; j <= i; j++) {
		gui_state.display_matrix[i + j * N_Beacons] = false;
	}

	ImGui::Separator();
	ImGui::Checkbox("Calibrating", &gui_state.calibrating);

	for (size_t i = 0; i < N_Beacons; ++i) {
		auto& b = state.beacons[i];
		ImGui::Text(
			"% 10.9lf +- % 10.9lf | % 10.9lf +- % 10.9lf | % 10.9lf +- % 10.9lf",
			b.mean.x, b.std.x,
			b.mean.y, b.std.y,
			b.mean.z, b.std.z
		);
	}
	ImGui::Separator();

	gui_state.want_compute = false;
	gui_state.want_reset_driver = false;
	gui_state.want_next_reading = false;

	gui_state.want_reset_driver = ImGui::Button("Reset driver");

	if (ImGui::CollapsingHeader("Simulation")) {
		ImGui::PushID("Simulation");
		ImGui::SliderSize("Resolution", &gui_state.space_sim.resolution, 1, 2000);
		ImGui::SliderAngled("Pivot angle", &gui_state.space_sim.pivot_angle, 180, 0);
		ImGui::SliderDouble("Pivot distance", &gui_state.space_sim.pivot_distance, 0, 0.1);
		ImGui::SliderDouble("R_in", &gui_state.space_sim.R_in, 0, gui_state.space_sim.R_out);
		ImGui::SliderDouble("R_out", &gui_state.space_sim.R_out, gui_state.space_sim.R_in, 0.05);
		ImGui::SliderDouble("h", &gui_state.space_sim.h, 0, 0.1);
		ImGui::SliderDouble("strength", &gui_state.space_sim.magnet_strength, 0, 10);
		if (ImGui::Button("Compute")) gui_state.want_compute = true;

		ImGui::Checkbox("X", &gui_state.display_x);
		ImGui::SameLine();
		ImGui::Checkbox("Y", &gui_state.display_y);
		ImGui::SameLine();
		ImGui::Checkbox("Z", &gui_state.display_z);
		ImGui::PopID();
	}

	if (state.space_res && ImGui::CollapsingHeader("Model")) {
		gui_state.want_next_reading |= ImGui::Button("Feed next readings");

		if (state.next_reading < state.readings.size()) {
			auto& next_read = state.readings[state.next_reading];

			for (size_t i = 0; i < N_Beacons; ++i) {
				auto& r = next_read.beacons[i];

				ImGui::Text("%zu: % 10.8lf, % 10.8lf, %10.8lf", i, r.x, r.y, r.z);
			}
		}


		ImGui::PlotLines("Section X", [](void* data, int idx) -> float {
			State& s = *(State*)data;

			// auto d = s.space_res->field[10 * idx + 500 * 1001];
			return (float)0;
		}, &state, 100, 0);
	}

	ImGui::End();

	if (state.space_res) {
		ImGui::Begin("Result");
		ImGui::Image(gui_state.field_texture);
		ImGui::End();
	}




	ImGui::End();
}


void render_plot(const std::vector<Reading>& readings) noexcept {
	static std::vector<std::vector<float>> lines;
	static std::vector<const float*> ys;
	static size_t N_Samples = 1000;
	constexpr size_t N_Comp = 2;

	static std::array<char*, N_Beacons * N_Comp> names = { nullptr };
	for (size_t i = 0; i < N_Beacons * N_Comp; ++i) {
		free(names[i]);
		names[i] = (char*)malloc(sizeof("Beacons XXXXXXX"));
		sprintf(names[i], "Beacon %d %d", (int)(i / N_Beacons), (int)(i % N_Beacons));
	}

	lines.clear();
	ys.clear();

	lines.resize(N_Beacons * N_Comp);

	float maxs[N_Beacons * N_Comp];
	float mins[N_Beacons * N_Comp];
	for (auto& x : maxs) x = -FLT_MAX;
	for (auto& x : mins) x = +FLT_MAX;
	size_t n = readings.size();

	for (size_t j = 0; j < N_Beacons; ++j) {
		for (size_t i = (n > N_Samples) ? (n - N_Samples) : 0; i < n; ++i) {
			auto& a = readings[i].beacons[j];
			lines[j * N_Comp + 0].push_back(a.x);
			lines[j * N_Comp + 1].push_back(a.y);
			// lines[j * 3 + 2].push_back(a.z);

			maxs[j] = (float)(lines[j].back() > maxs[j] ? lines[j].back() : maxs[j]);
			mins[j] = (float)(lines[j].back() < mins[j] ? lines[j].back() : mins[j]);
		}
	}
	for (size_t i = 0; i < lines.size(); ++i) for (auto& y : lines[i]) {
		maxs[i] = (float)(y > maxs[i] ? y : maxs[i]);
		mins[i] = (float)(y < mins[i] ? y : mins[i]);
	}

	for (size_t i = 0; i < lines.size(); ++i) ys.push_back(lines[i].data());

	ImGui::PlotConfig conf;
	conf.values.ys_list = ys.data();
	conf.values.ys_count = ys.size();
	conf.values.count = (int)std::min(readings.size(), N_Samples);

	conf.scale.min = 0;
	conf.scale.max = 0;
	for (auto& x : mins) conf.scale.min = std::min(conf.scale.min, x);
	for (auto& x : maxs) conf.scale.max = std::max(conf.scale.max, x);
	conf.scale.min *= 0.9f;
	conf.scale.max *= 1.1f;

	conf.tooltip.show = true;

	conf.grid_y.show = true;
	conf.grid_y.size = std::max(0.f, conf.scale.max) / 5 + 1e-5;
	conf.grid_y.subticks = 5;

	conf.line_thickness = 2.f;

	conf.frame_size = { ImGui::GetWindowWidth() * 0.95f, 100 };

	conf.tooltip.ys_names = names.data();
	conf.tooltip.format = "%s, %g: % 12.9f";

	ImGui::Plot("Readings", conf);

	for (size_t i = 0; i < N_Beacons * N_Comp; ++i) {
		conf.values.ys_list = ys.data() + i;
		conf.values.ys_count = 1;
		// conf.tooltip.ys_names = names.data() + i;
		// conf.scale.min = 0.9f * std::min(mins[i], 0.f);
		// conf.scale.max = 1.1f * std::max(maxs[i], 0.f);

		// conf.grid_y.size = (conf.scale.max - conf.scale.min + 1e-5) / 5;

		ImGui::Plot(names[i], conf);
	}

	ImU32 colors[3] = {
		IM_COL32(255, 0, 0, 255), IM_COL32(0, 255, 0, 255), IM_COL32(0, 0, 255, 255)
	};
	
	{
		ys.resize(2);
		lines.resize(2);
		for (auto& l : lines) l.clear();

		conf.scale.max = -FLT_MAX;
		conf.scale.min = +FLT_MAX;

		constexpr float alphag = 0.9;
		Vector3d g = { 0 };
		for (size_t j = (n > N_Samples) ? (n - N_Samples) : 0; j < n; ++j) {
			Vector3d acc = { 0 };
			for (size_t i = 0; i < N_Imus; ++i) acc += readings[j].accel[i];
			acc /= 2;

			g.x = alphag * g.x + (1 - alphag) * acc.x;
			g.y = alphag * g.y + (1 - alphag) * acc.y;
			g.z = alphag * g.z + (1 - alphag) * acc.z;

			float a = 0;
			a += (acc.x - g.x) * (acc.x - g.x);
			a += (acc.y - g.y) * (acc.y - g.y);
			a += (acc.z - g.z) * (acc.z - g.z);
			a = std::sqrtf(a);

			float G = g.x * g.x + g.y * g.y + g.z * g.z;
			G = std::sqrtf(G);

			conf.scale.max = std::max(conf.scale.max, a);
			conf.scale.max = std::max(conf.scale.max, G);
			conf.scale.min = std::min(conf.scale.min, a);
			conf.scale.min = std::min(conf.scale.min, G);

			lines[0].push_back(a);
			lines[1].push_back(G);
		}

		conf.scale.min *= 0.9f;
		conf.scale.max *= 1.1f;

		conf.grid_y.size = std::max(1e-5f, conf.scale.max - conf.scale.min) / 5;

		ys[0] = lines[0].data();
		ys[1] = lines[1].data();
		conf.values.ys_list = ys.data();
		conf.values.ys_count = ys.size();

		conf.values.colors = colors;

		ImGui::Plot("Acceleration", conf);
	}
{
		ys.resize(3);
		lines.resize(3);
		for (auto& l : lines) l.clear();

		conf.scale.max = -FLT_MAX;
		conf.scale.min = +FLT_MAX;

		for (size_t j = (n > N_Samples) ? (n - N_Samples) : 0; j < n; ++j) {
			Vector3d gyr = { 0 };
			for (size_t i = 0; i < N_Imus; ++i) gyr += readings[j].gyro[i];
			gyr /= 2;

			conf.scale.max = std::max(conf.scale.max, (float)gyr.x);
			conf.scale.max = std::max(conf.scale.max, (float)gyr.y);
			conf.scale.max = std::max(conf.scale.max, (float)gyr.z);

			conf.scale.min = std::min(conf.scale.min, (float)gyr.x);
			conf.scale.min = std::min(conf.scale.min, (float)gyr.y);
			conf.scale.min = std::min(conf.scale.min, (float)gyr.z);

			lines[0].push_back(gyr.x);
			lines[1].push_back(gyr.y);
			lines[2].push_back(gyr.z);
		}

		conf.scale.min *= 0.9f;
		conf.scale.max *= 1.1f;

		conf.grid_y.size = std::max(1e-5f, conf.scale.max - conf.scale.min) / 5;

		ys[0] = lines[0].data();
		ys[1] = lines[1].data();
		ys[2] = lines[2].data();
		conf.values.ys_list = ys.data();
		conf.values.ys_count = ys.size();

		conf.values.colors = colors;

		conf.tooltip.ys_names = names.data();

		ImGui::Plot("Gyroscope", conf);
	}

	ImGui::SliderSize("#Samples", &N_Samples, 10, 1000);
}

Debug_Values frame_debug_values;

void Debug_Values::reset() noexcept {
	histograms.clear();
}

void Debug_Values::add_to_distribution(const char* name, double x) noexcept {
	if (histograms.count(name) == 0) {
		histograms[name].reserve(1000);
	}

	histograms[name].push_back(x);
}

void Debug_Values::watch(const char* name, double x) noexcept {
	values[name] = x;
}

void render(Debug_Values& debug_values) noexcept {
	thread_local bool debug_opened = true;
	if (ImGui::Begin("Debug values", &debug_opened)) {
		for (auto& [name, values] : debug_values.histograms) if (values.size() > 1) {

			ImGui::Separator();
			ImGui::Text("%s", name);

			auto min_v = values.front();
			auto max_v = values.front();
			auto avg_v = 0.0;

			for (auto& x : values) min_v = std::min(x, min_v);
			for (auto& x : values) max_v = std::max(x, max_v);
			for (auto& x : values) avg_v += x;

			avg_v /= values.size();

			auto n_bins = (size_t)(std::log2(values.size()) + 0.5);
			auto dt = (max_v - min_v) / n_bins;

			struct User_Ptr { std::vector<double>& values; double min_v; double max_v; double dt; };

			User_Ptr user_ptr = { values, min_v, max_v, dt };

			ImGui::Text("Min Avg Max: % 10.8lf % 10.8lf % 10.8lf", min_v, avg_v, max_v);
			ImGui::PlotHistogram("", [](void* data, int idx) -> float {
				User_Ptr& user_ptr = *(User_Ptr*)data;

				double s = 0;
				for (auto& x : user_ptr.values) {
					if (
						user_ptr.min_v + idx * user_ptr.dt <= x &&
						x < user_ptr.min_v + (idx + 1) * user_ptr.dt
					) s += 1;
				}

				return s;
			}, &user_ptr, n_bins, 0, 0, FLT_MAX, FLT_MAX, {0, 200});
		} else if (values.size() == 1) {
			ImGui::Separator();
			ImGui::Text("%s %10.8lf", name, values.front());
		}

		ImGui::Separator();
		ImGui::Separator();

		for (auto& [n, x] : debug_values.values) ImGui::Text("%s %10.8lf", n, x);
	}
	ImGui::End();
}