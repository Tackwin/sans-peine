#include "Ease.hpp"

EASE_WATCH_ME;

Build build(Flags flags) noexcept {
	auto b = Build::get_default(flags);

	b.name = "Client";

	b.add_define("SFML_STATIC");
	b.add_source_recursively("./src/");
	b.add_header("./src/");
	b.add_header("./include/");

	b.add_library_path("./lib/");
	b.add_library("sfml-graphics-s");
	b.add_library("sfml-system-s");
	b.add_library("sfml-window-s");
	b.add_library("flac");
	b.add_library("openal32");
	b.add_library("vorbis");
	b.add_library("opengl32");
	b.add_library("gdi32");
	b.add_library("Advapi32");
	b.add_library("Winmm");

	return b;
}
