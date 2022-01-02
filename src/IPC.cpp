#include "IPC.hpp"

#include <stdio.h>

HANDLE open_slot(const char* mail_name) noexcept {
	auto mail_slot = CreateFileA(
		mail_name,
		GENERIC_WRITE,
		FILE_SHARE_READ,
		nullptr,
		OPEN_EXISTING,
		FILE_ATTRIBUTE_NORMAL,
		nullptr
	);

	if (mail_slot == INVALID_HANDLE_VALUE) {
		fprintf(stderr, "Error while opening mail slot\n");
		fprintf(stderr, "CreateFile failed with %d.\n", GetLastError());
		return 0;
	}

	return mail_slot;
}

HANDLE make_slot(const char* mail_name) noexcept {
	auto mail_slot = CreateMailslotA(mail_name, 0, 1, nullptr);

	if (mail_slot == INVALID_HANDLE_VALUE) {
		fprintf(stderr, "Error creating mail slot\n");
		return 0;
	}

	return mail_slot;
}