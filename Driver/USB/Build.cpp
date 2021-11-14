#include "Ease.hpp"


EASE_WATCH_ME;

Build build(Flags flags) noexcept {
	auto b = Build::get_default(flags);

	b.name = "Driver";

	b.add_source("src/main.cpp");

	return b;
}