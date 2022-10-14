#pragma once

#include "stdint.h"

enum Register {
	POLL   = 0x00,
	CMM    = 0x01,
	CCX    = 0x04,
	CCY    = 0x06,
	CCZ    = 0x08,
	TMRC   = 0x0B,
	MX     = 0x24,
	MY     = 0x27,
	MZ     = 0x2A,
	BIST   = 0x33,
	STATUS = 0x34,
	HSHAKE = 0x35,
	REVID  = 0x36,
	READ   = 0x80
};

enum Data_Rate {
	HZ_600   = 0x92,
	HZ_300   = 0x93,
	HZ_150   = 0x94,
	HZ_75    = 0x95,
	HZ_37    = 0x96,
	HZ_18    = 0x97,
	HZ_9     = 0x98,
	HZ_4_5   = 0x99,
	HZ_2_3   = 0x9A,
	HZ_1_2   = 0x9B,
	HZ_0_6   = 0x9C,
	HZ_0_3   = 0x9D,
	HZ_0_15  = 0x9E,
	HZ_0_075 = 0x9F
};

struct RM3100 {
	bool is_connected() noexcept;
	bool set_cycle_count(uint16_t ccx = 0, uint16_t ccy = 0, uint16_t ccz = 0) noexcept;
	bool get_cycle_count(uint16_t* ccx = nullptr, uint16_t* ccy = nullptr, uint16_t* ccz = nullptr) noexcept;
	bool toggle_continuous_mode(bool cmx = true, bool cmy = true, bool cmz = true) noexcept;
	bool read_continuous_mode() noexcept;
	bool set_update_rate(Data_Rate data_rate) noexcept;
	bool is_data_ready() noexcept;

	bool revid(uint8_t& out) noexcept;

	bool read_ut(float* x = nullptr, float* y = nullptr, float* z = nullptr) noexcept;

	static bool self_test(uint8_t port) noexcept;

	uint8_t address = 0b0100000;
	bool in_continuous_mode = false;
	double LSB_ut_x = 75;
	double LSB_ut_y = 75;
	double LSB_ut_z = 75;
};