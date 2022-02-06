#define NAME usb
#include "Ease.hpp"

EASE_WATCH_ME;

Build usb(Flags flags) noexcept {
	auto b = Build::get_default(flags);

	b.target = Ease::Build::Target::Shared;

	b.name = "Driver";

	b.add_source("Driver/USB/src/main.cpp");
	b.add_source_recursively("src/");
	b.add_header("src/");

	b.add_define("_CRT_SECURE_NO_WARNINGS");
	b.add_library("opengl32");
	b.add_library("gdi32");
	b.add_library("Advapi32");
	b.add_library("Winmm");

	return b;
}