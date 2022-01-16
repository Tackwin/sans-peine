#include "IPC.hpp"

#include <stdio.h>

HANDLE IPC::open_slot(const char* mail_name) noexcept {
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

HANDLE IPC::make_slot(const char* mail_name) noexcept {
	auto mail_slot = CreateMailslotA(mail_name, 0, 0, nullptr);

	if (mail_slot == INVALID_HANDLE_VALUE) {
		fprintf(stderr, "Error creating mail slot\n");
		return 0;
	}

	return mail_slot;
}

bool IPC::read(HANDLE mail_slot, void* dst, size_t sz) noexcept {
	DWORD cbRead = 0;
	BOOL fResult;

	fResult = ReadFile(mail_slot, dst, sz, &cbRead, nullptr);

	if (!fResult) return false;

	return cbRead == sz;
}

bool IPC::write(HANDLE mail_slot, void* src, size_t sz) noexcept {
	DWORD cbWritten;

	auto fResult = WriteFile(mail_slot, src, sz, &cbWritten, nullptr);

	if (!fResult) return false;
	return cbWritten == sz;
}
