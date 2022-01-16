#pragma once

#include "Windows.h"

namespace IPC {
	extern HANDLE make_slot(const char* mail_name) noexcept;
	extern HANDLE open_slot(const char* mail_name) noexcept;

	extern bool read(HANDLE mail_slot, void* dst, size_t sz) noexcept;
	extern bool write(HANDLE mail_slot, void* src, size_t sz) noexcept;
}
