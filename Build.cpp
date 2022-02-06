#define NAME sans_peine
#include "Ease.hpp"

EASE_WATCH_ME;

#undef EASE_WATCH_ME
#define EASE_WATCH_ME int i##__COUNTER__

#include "Driver/USB/Build.cpp"
#include "Client/Build.cpp"

Build sans_peine(Flags flags) noexcept {
	bool run_last = flags.run_after_compilation;
	flags.run_after_compilation = false;
	auto usb_build = usb(flags);
	usb_build.target = Build::Target::Shared;

	flags.run_after_compilation = run_last;
	auto client_build = client(flags);

	return Build::sequentials({usb_build, client_build});
}
