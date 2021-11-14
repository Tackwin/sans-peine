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

#pragma pack(1)
struct Reading {
	double beacon1;
	double beacon2;
	bool pressed;
};
struct State {
	size_t n_beacons_placed = 0;
	std::array<sf::Vector2f, 2> beacons;
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
};


void toggle_fullscreen(State& state) noexcept;

void make_slot(State& state) noexcept;
std::optional<Reading> read_mail(State& state) noexcept;
void update(State& state) noexcept;
void render(State& state) noexcept;
void render_plot(const std::vector<Reading>& readings) noexcept;

sf::Vector2f triangulate(State& state, Reading r) noexcept;

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
		if (state.n_beacons_placed >= 2) {
			state.n_beacons_placed = 0;
		} else {
			state.beacons[state.n_beacons_placed] = {
				(float)sf::Mouse::getPosition(*state.window).x,
				(float)sf::Mouse::getPosition(*state.window).y
			};
			state.n_beacons_placed++;
		}
	}

	auto res = read_mail(state);

	if (res) {
		avg_readings.push_back(*res);
		if (avg_readings.size() >= state.oversampling) {
			Reading avg = {};
			size_t n_pressed = 0;
			for (auto& x : avg_readings) {
				avg.beacon1 += x.beacon1;
				avg.beacon2 += x.beacon2;
				n_pressed += (x.pressed ? 1 : 0);
			}
			avg.beacon1 /= avg_readings.size();
			avg.beacon2 /= avg_readings.size();
			avg.pressed = n_pressed > avg_readings.size() / 2;

			avg_readings.clear();

			state.readings.push_back(avg);
		}

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

	if (ImGui::Button("Clear")) state.readings.clear();

	ImGui::End();

	ImGui::End();

	sf::CircleShape shape;
	shape.setRadius(10);

	for (size_t i = 0; i < state.n_beacons_placed; ++i) {
		shape.setPosition(
			state.beacons[i].x - shape.getRadius(), state.beacons[i].y - shape.getRadius()
		);
		shape.setFillColor({255, 0, 0, 255});
		state.renderTarget->draw(shape);
	}

	for (size_t i = 0; i + 1 < state.readings.size(); ++i) {
		auto curr = state.readings[i];
		auto next = state.readings[i + 1];
		if (!next.pressed) continue;

		sf::Vertex line[] = {
			sf::Vertex(triangulate(state, curr)),
			sf::Vertex(triangulate(state, next))
		};

		state.renderTarget->draw(line, 2, sf::Lines);
	}

	if (state.n_beacons_placed >= 2) {
		sf::Vector2f b1 = state.beacons[0];
		sf::Vector2f b2 = state.beacons[1];

		double dt_beacons = (b1.x - b2.x) * (b1.x - b2.x) + (b1.y - b2.y) * (b1.y - b2.y);
		dt_beacons = std::sqrt(dt_beacons);
		shape.setRadius(
			(dt_beacons * state.info_distance / state.beacon_distance) / 2
		);
		shape.setPosition(
			(state.beacons[0].x + state.beacons[1].x) / 2 - shape.getRadius(),
			(state.beacons[0].y + state.beacons[1].y) / 2 - shape.getRadius()
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

	static const char* names[2] = {
		"Beacon 1",
		"Beacon 2"
	};

	lines.clear();
	ys.clear();

	lines.resize(2);

	float max1 = 0;
	float max2 = 0;
	size_t n = readings.size();
	for (size_t i = (n > N_Samples) ? (n - N_Samples) : 0; i < n; ++i) {
		auto& [a, _, __] = readings[i];
		lines[0].push_back((float)a);
		max1 = (float)(a > max1 ? a : max1);
	}
	for (size_t i = (n > N_Samples) ? (n - N_Samples) : 0; i < n; ++i) {
		auto& [_, b, __] = readings[i];
		lines[1].push_back((float)b);
		max2 = (float)(b > max2 ? b : max2);
	}

	ys.push_back(lines[0].data());
	ys.push_back(lines[1].data());

	ImGui::PlotConfig conf;
	conf.values.ys_list = ys.data();
	conf.values.ys_count = 2;
	conf.values.count = (int)std::min(readings.size(), N_Samples);

	conf.scale.min = 0;
	conf.scale.max = 1.1f * (max1 > max2 ? max1 : max2);

	conf.tooltip.show = true;

	conf.grid_y.show = true;
	conf.grid_y.size = conf.scale.max / 5;
	conf.grid_y.subticks = 5;

	conf.line_thickness = 2.f;

	conf.frame_size = { ImGui::GetWindowWidth() * 0.95f, 100 };

	conf.tooltip.ys_names = names;
	conf.tooltip.format = "%s, %g: % 12.9f";

	ImGui::Plot("Readings", conf);

	conf.values.ys_list = ys.data();
	conf.values.ys_count = 1;
	conf.scale.max = 1.1f * max1;

	conf.grid_y.size = conf.scale.max / 5;

	ImGui::Plot("Reading 1", conf);

	conf.values.ys_list = ys.data() + 1;
	conf.values.ys_count = 1;
	conf.scale.max = 1.1f * max2;

	conf.grid_y.size = conf.scale.max / 5;

	ImGui::Plot("Reading 2", conf);

	ImGui::SliderSize("#Samples", &N_Samples, 10, 1000);
}


double B(double M, sf::Vector2f r) noexcept {
	constexpr double u0 = 1.225663753e-6;
	constexpr double PI = 3.141592653589793238462643383279502884;

	double d = std::sqrt(r.x * r.x + r.y * r.y);

	return u0 * M / (4 * PI * d * d * d);
}

void make_slot(State& state) noexcept {
	state.mail_slot = CreateMailslotA(Mail_Name, 0, MAILSLOT_WAIT_FOREVER, nullptr);

	if (state.mail_slot == INVALID_HANDLE_VALUE) {
		fprintf(stderr, "Error creating mail slot\n");
		std::terminate();
	}
}

std::optional<Reading> read_mail(State& state) noexcept {
	DWORD cbMessage = 0;
	DWORD cMessage = 0;
	DWORD cbRead = 0; 
    BOOL fResult; 
    LPTSTR lpszBuffer; 
    TCHAR achID[80]; 
    DWORD cAllMessages; 
    OVERLAPPED ov;

    HANDLE hEvent = CreateEvent(NULL, FALSE, FALSE, TEXT("SP_Event"));
    if (!hEvent) return std::nullopt;

    ov.Offset = 0;
    ov.OffsetHigh = 0;
    ov.hEvent = hEvent;

    fResult = GetMailslotInfo(
    	state.mail_slot, // mailslot handle 
    	nullptr,
        &cbMessage,                   // size of next message 
        &cMessage,                    // number of messages 
        nullptr
    );              // no read time-out 
 
    if (!fResult) {
        fprintf(stderr, "GetMailslotInfo failed with %d.\n", (int)GetLastError()); 
        return std::nullopt;
    }

    if (cbMessage == MAILSLOT_NO_MESSAGE) {
        return std::nullopt;
    } 

    cAllMessages = cMessage;

    // Allocate memory for the message.

    lpszBuffer = (LPTSTR)GlobalAlloc(GPTR, lstrlen((LPTSTR) achID) * sizeof(TCHAR) + cbMessage);

    if(!lpszBuffer) return std::nullopt;

    lpszBuffer[0] = '\0'; 

    fResult = ReadFile(state.mail_slot, lpszBuffer, cbMessage, &cbRead, &ov);

    if (!fResult) {
        fprintf(stderr, "ReadFile failed with %d.\n", (int)GetLastError());
        GlobalFree((HGLOBAL) lpszBuffer);
        return std::nullopt;
    }

    // Concatenate the message and the message-number string. 


    GlobalFree((HGLOBAL) lpszBuffer); 

    fResult = GetMailslotInfo(
    	state.mail_slot,              // mailslot handle 
        (LPDWORD) NULL,               // no maximum message size 
        &cbMessage,                   // size of next message 
        &cMessage,                    // number of messages 
        (LPDWORD) NULL                // no read time-out
    );

    if (!fResult) {
        fprintf(stderr, "GetMailslotInfo failed (%d)\n", (int)GetLastError());
        return std::nullopt;
    }

    CloseHandle(hEvent);

    return *reinterpret_cast<Reading*>(lpszBuffer);
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

sf::Vector2f triangulate(State& state, Reading r) noexcept {
	constexpr double u0 = 1.225663753e-6;
	constexpr double PI = 3.141592653589793238462643383279502884;

	if (state.n_beacons_placed < 2) return {};

	sf::Vector2f b1 = state.beacons[0];
	sf::Vector2f b2 = state.beacons[1];

	double dt_beacons = (b1.x - b2.x) * (b1.x - b2.x) + (b1.y - b2.y) * (b1.y - b2.y);

	dt_beacons = std::sqrt(dt_beacons);

	double d1 = state.magnet_strength * u0 / (4 * PI * r.beacon1);
	d1 = std::pow(d1, 1.0 / 3.0);
	d1 = dt_beacons * d1 / state.beacon_distance;

	double d2 = state.magnet_strength * u0 / (4 * PI * r.beacon2);
	d2 = std::pow(d2, 1.0 / 3.0);
	d2 = dt_beacons * d2 / state.beacon_distance;

	auto res = *circle_circle_intersection(b1, d1, b2, d2);
	return res;
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
