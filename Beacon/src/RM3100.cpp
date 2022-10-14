#include "RM3100.hpp"
#include <Wire.h>
#include <SPI.h>
#include <Arduino.h>

// #define USE_SPI

extern void serial_printf(const char *fmt, ...);
bool RM3100::is_connected() noexcept {
	#ifdef USE_SPI
		uint8_t out;
		if (!revid(out)) return false;
		return true;
	#else
		Wire.beginTransmission(address);
		int end = Wire.endTransmission();
		return end == 0;
	#endif
}

bool RM3100::self_test(uint8_t port) noexcept {
	auto write = [port] (uint8_t addr, uint8_t data) -> bool {
		digitalWrite(port, LOW);
		SPI.transfer(addr);
		SPI.transfer(data);
		digitalWrite(port, HIGH);
	};

	auto read = [port] (uint8_t addr) -> uint8_t {
		digitalWrite(port, LOW);
		uint8_t ret = SPI.transfer(addr | Register::READ);
		digitalWrite(port, HIGH);
		return ret;
	};
	write(Register::CMM, 112);
	serial_printf("-");
	write(Register::BIST, 0x08);
	serial_printf("A");
	write(Register::POLL, 0x70);
	serial_printf("B");

	while (true) {
		uint8_t status = read(Register::STATUS);
		serial_printf("C");

		if (!(status & (1 << 7))) continue;
		status = read(Register::BIST);
		serial_printf("D");
		if (!(status & (1 << 7))) continue;

		serial_printf("\nTest result of %d: %d\n", (int)port, (int)status);
		write(Register::BIST, 0);
		break;
	}
	return true;
}

bool RM3100::set_cycle_count(uint16_t ccx, uint16_t ccy, uint16_t ccz) noexcept {
	bool success = true;

	#ifdef USE_SPI
		uint8_t data[6];
		size_t count = 0;
		if (ccx > 0) {
			data[count++] = ((uint8_t)((ccx & 0xFF00) >> 8));
			data[count++] = ((uint8_t)((ccx & 0x00FF) >> 0));
			LSB_ut_x = 0.3671 * ccx + 1.5;
		}
		if (ccy > 0) {
			data[count++] = ((uint8_t)((ccy & 0xFF00) >> 8));
			data[count++] = ((uint8_t)((ccy & 0x00FF) >> 0));
			LSB_ut_y = 0.3671 * ccy + 1.5;
		}
		if (ccz > 0) {
			data[count++] = ((uint8_t)((ccz & 0xFF00) >> 8));
			data[count++] = ((uint8_t)((ccz & 0x00FF) >> 0));
			LSB_ut_z = 0.3671 * ccz + 1.5;
		}
		SPI.transfer(Register::CCX);
		SPI.transfer(data[0]);
		SPI.transfer(data[1]);
		SPI.transfer(data[2]);
		SPI.transfer(data[3]);
		SPI.transfer(data[4]);
		SPI.transfer(data[5]);
	#else
		if (ccx > 0) {
			
			Wire.beginTransmission(address);
			Wire.write((uint8_t)Register::CCX);
			Wire.write((uint8_t)((ccx & 0xFF00) >> 8));
			Wire.write((uint8_t)((ccx & 0x00FF) >> 0));
			// success &= Wire.endTransmission() == 0;
			LSB_ut_x = 0.3671 * ccx + 1.5;
		}
		if (ccy > 0) {
			// Wire.beginTransmission(address);
			// Wire.write((uint8_t)Register::CCY);
			Wire.write((uint8_t)((ccy & 0xFF00) >> 8));
			Wire.write((uint8_t)((ccy & 0x00FF) >> 0));
			// success &= Wire.endTransmission() == 0;
			LSB_ut_y = 0.3671 * ccy + 1.5;
		}
		if (ccz > 0) {
			// Wire.beginTransmission(address);
			// Wire.write((uint8_t)Register::CCZ);
			Wire.write((uint8_t)((ccz & 0xFF00) >> 8));
			Wire.write((uint8_t)((ccz & 0x00FF) >> 0));
			success &= Wire.endTransmission() == 0;
			LSB_ut_z = 0.3671 * ccz + 1.5;
		}
	#endif

	return success;
}

bool RM3100::get_cycle_count(uint16_t* ccx = nullptr, uint16_t* ccy = nullptr, uint16_t* ccz = nullptr) noexcept {
	bool success = true;
	
	#ifdef USE_SPI
		uint8_t data[6];
		size_t count = 0;
		uint8_t x = SPI.transfer(Register::CCX | Register::READ);
		serial_printf("A %d\n", (int)x);
		SPI.transfer(data, count);
		serial_printf("azeaze");
		if (ccx) {
			*ccx = 0;
			*ccx |= data[0];
			*ccx |= data[1];
		}
		if (ccy) {
			*ccy = 0;
			*ccy |= data[2];
			*ccy |= data[3];
		}
		if (ccz) {
			*ccz = 0;
			*ccz |= data[4];
			*ccz |= data[5];
		}
	#else
		if (ccx) {
			Wire.beginTransmission(address);
			Wire.write((uint8_t)Register::CCX);
			success &= Wire.endTransmission() == 0;
			Wire.requestFrom(address, 2);
			*ccx = 0;
			*ccx |= (Wire.read() << 8);
			*ccx |= (Wire.read() << 0);
		}
		if (ccz) {
			Wire.beginTransmission(address);
			Wire.write((uint8_t)Register::CCY);
			success &= Wire.endTransmission() == 0;
			Wire.requestFrom(address, 2);
			*ccy = 0;
			*ccy |= (Wire.read() << 8);
			*ccy |= (Wire.read() << 0);
		}
		if (ccz) {
			Wire.beginTransmission(address);
			Wire.write((uint8_t)Register::CCZ);
			success &= Wire.endTransmission() == 0;
			Wire.requestFrom(address, 2);
			*ccz = 0;
			*ccz |= (Wire.read() << 8);
			*ccz |= (Wire.read() << 0);
		}
	#endif
	return success;
}


bool RM3100::toggle_continuous_mode(bool cmx, bool cmy, bool cmz) noexcept {
	bool start = !in_continuous_mode;

	uint8_t byte = 0;
	byte |= (start ? 1 : 0) << 0; // Continuous On/Off
	byte |= (0)             << 1; // Reserved
	byte |= (0)             << 2; // DRDY 0 => wait for full axis 1 => wait for any
	byte |= (1)             << 3; // Reserved
	byte |= (cmx ? 1 : 0)   << 4; // CMX
	byte |= (cmy ? 1 : 0)   << 5; // CMY
	byte |= (cmz ? 1 : 0)   << 6; // CMZ
	byte |= (0)             << 7; // Reserved

	#ifdef USE_SPI
		SPI.transfer(Register::CMM);
		SPI.transfer(byte);
		return true;
	#else
		Wire.beginTransmission(address);
		Wire.write((uint8_t)Register::CMM);
		Wire.write(byte);
		if (Wire.endTransmission() == 0) {
			in_continuous_mode ^= true;
		}
	#endif

	#ifdef USE_SPI
		// SPI.transfer(Register::CMM | Register::READ);
		// int read = SPI.transfer(0);
	#else
		Wire.beginTransmission(address);
		Wire.write((uint8_t)(Register::CMM));
		Wire.requestFrom(address, (uint8_t)1);
		while(!Wire.available()) {};
		int read = Wire.read();
	#endif
	return Wire.endTransmission() == 0;
}

bool RM3100::read_continuous_mode() noexcept {

	for (uint8_t a = 0; a < 0xfe; ++a) {
		Wire.beginTransmission(address);
		Wire.write(a);
		Wire.endTransmission();
		delay(100);
		int read = Wire.requestFrom(address, (uint8_t)1);
		read = Wire.read();
		serial_printf("XX 0x%x Read 0x%x\n", (int)a, (int)read);
		Wire.endTransmission() == 0;
	}	

	return Wire.endTransmission() == 0;
}

bool RM3100::set_update_rate(Data_Rate data_rate) noexcept {
	uint8_t TMRC = (uint8_t)data_rate;
	uint8_t byte = 0b1001 << 4;
	byte |= TMRC & 0b1111;

	#ifdef USE_SPI
		SPI.transfer(Register::TMRC);
		SPI.transfer(byte);
	#else
		Wire.beginTransmission(address);
		Wire.write((uint8_t)Register::TMRC);
		Wire.write(byte);
		return Wire.endTransmission() == 0;
	#endif
}

bool RM3100::is_data_ready() noexcept {
	#ifdef USE_SPI
		uint8_t x = SPI.transfer(Register::STATUS | Register::READ);
		uint8_t value = SPI.transfer(0);
	#else
		Wire.beginTransmission(address);
		Wire.write((uint8_t)Register::STATUS);
		if (Wire.endTransmission() != 0) return false;
		Wire.requestFrom(address, (uint8_t)1);
		uint8_t value = Wire.read();
	#endif
	return value & 0x80;
}

bool RM3100::read_ut(float* x, float* y, float* z) noexcept {
	#ifdef USE_SPI
		uint8_t bus = Register::MX | Register::READ;
		SPI.transfer(bus++);
		uint8_t x2 = SPI.transfer(bus++);
		uint8_t x1 = SPI.transfer(bus++);
		uint8_t x0 = SPI.transfer(bus++);
		uint8_t y2 = SPI.transfer(bus++);
		uint8_t y1 = SPI.transfer(bus++);
		uint8_t y0 = SPI.transfer(bus++);
		uint8_t z2 = SPI.transfer(bus++);
		uint8_t z1 = SPI.transfer(bus++);
		uint8_t z0 = SPI.transfer(bus++);
		bool success = true;
	#else
		Wire.beginTransmission(address);
		Wire.write((uint8_t)(Register::MX));
		bool success = Wire.endTransmission() == 0;

		// Request 9 bytes from the measurement results registers
		Wire.requestFrom(address, (uint8_t)9);
		while (Wire.available() != 9);
		uint8_t x2 = Wire.read();
		uint8_t x1 = Wire.read();
		uint8_t x0 = Wire.read();
		uint8_t y2 = Wire.read();
		uint8_t y1 = Wire.read();
		uint8_t y0 = Wire.read();
		uint8_t z2 = Wire.read();
		uint8_t z1 = Wire.read();
		uint8_t z0 = Wire.read();
	#endif
	int32_t ix = 0;
	int32_t iy = 0;
	int32_t iz = 0;

	//special bit manipulation since there is not a 24 bit signed int data type
	if (x2 & 0x80) ix = 0xFF;
	if (y2 & 0x80) iy = 0xFF;
	if (z2 & 0x80) iz = 0xFF;

	//format results into single 32 bit signed value
	ix = (ix * 256 * 256 * 256) | (int32_t)(x2) * 256 * 256 | (uint16_t)(x1) * 256 | x0;
	iy = (iy * 256 * 256 * 256) | (int32_t)(y2) * 256 * 256 | (uint16_t)(y1) * 256 | y0;
	iz = (iz * 256 * 256 * 256) | (int32_t)(z2) * 256 * 256 | (uint16_t)(z1) * 256 | z0;


	if (x) *x = (float)(ix / LSB_ut_x);
	if (y) *y = (float)(iy / LSB_ut_y);
	if (z) *z = (float)(iz / LSB_ut_z);

	return success;
}

bool RM3100::revid(uint8_t& out) noexcept {
	#ifdef USE_SPI
		uint8_t read = SPI.transfer(Register::REVID | Register::READ);
		out = SPI.transfer(0);
	#else
		Wire.beginTransmission(address);
		Wire.write((uint8_t)(Register::REVID));
		bool succ = Wire.endTransmission() == 0;
		if (!succ) return false;
		Wire.requestFrom(address, (uint8_t)1);
		while(!Wire.available()) {};
		out = Wire.read();
	#endif
	return true;
}
 