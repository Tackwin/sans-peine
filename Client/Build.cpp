#define NAME client
#include "Ease.hpp"

Build client(Flags flags) noexcept {
	auto b = Build::get_default(flags);

	b.name = "Client";
	b.flags.fast_math = true;

	b.add_define("SFML_STATIC");
	b.add_source_recursively("Client/src/");
	b.add_header("Client/src/");
	b.add_header("Client/include/");
	b.add_source_recursively("src/");
	b.add_header("src/");

	b.add_library_path("Client/lib/");
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

