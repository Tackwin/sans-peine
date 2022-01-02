#pragma once

#include "Windows.h"

extern HANDLE make_slot(const char* mail_name) noexcept;
extern HANDLE open_slot(const char* mail_name) noexcept;
extern HANDLE Mail_Slot;