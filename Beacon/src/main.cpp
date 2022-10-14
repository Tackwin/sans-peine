#include "Arduino.h"
#include "Wire.h"
#include "HMC5883L.h"
#include "TCA9548.h"
#include "MLX90393.h"
#include "RM3100.hpp"
#include "GY521.h"

#include <math.h>

constexpr size_t N_Beacons = 2;
constexpr size_t N_Imus = 0;
constexpr size_t N_Sync_Seq = 16;

// HMC5883L
// 2 mG RMS
// 75 hz
// High availability, bad perf

// MMC5983MA
// 0.4/0.6/0.8/1.2 mG RMS
// 50/100/225/580 hz

// MMC34160PJ
// 1.5 mG RMS
// 125 hZ
// +- 16 G

void serial_printf(const char *fmt, ...) {
	va_list va;
	va_start(va, fmt);
	char buf[vsnprintf(NULL, 0, fmt, va) + 1];
	vsprintf(buf, fmt, va);
	Serial.print(buf);
	Serial.flush();
	va_end(va);
}

void use_bus(uint8_t bus) noexcept {
	Wire.beginTransmission(0x70);  // TCA9548A address is 0x70
	Wire.write(1 << bus);          // send byte to select bus
	Wire.endTransmission();
}


RM3100 beacons[N_Beacons];
Adafruit_MLX90393 mlx_beacons[N_Beacons];
bool beacon_healthy[N_Beacons] = { };

TCA9548 multiplexer(0x70);
size_t BEACON_BUS_MAP[] = {2, 7, 7, 6, 6, 7, 6, 7};

GY521 imus[N_Imus];
size_t IMU_BUS_MAP[] = {2, 3, 3, 3, 3, 3, 3, 3};
bool imu_healthy[N_Imus] = { };

constexpr size_t ROLLING = 1;
int last_vectors_idx[N_Beacons];
Vectori16 last_vectors[N_Beacons * ROLLING];

void select(uint8_t port) noexcept {
	multiplexer.selectChannel(port);
	// for (size_t i = 0; i < N_Beacons; i += 1) {
	// 	digitalWrite(BEACON_BUS_MAP[i], BEACON_BUS_MAP[i] != port);
	// }
	// delay(1);
}
void unselect() noexcept{
	// for (size_t i = 0; i < N_Beacons; i += 1) digitalWrite(BEACON_BUS_MAP[i], HIGH);
	// delay(1);
}

void setup() {
	SPI.begin(); // Initiate the SPI library
	SPI.beginTransaction(SPISettings(1000000, MSBFIRST, SPI_MODE0));  
	Serial.begin(128000);
	Wire.setWireTimeout(0);
	Wire.begin();
	bool res = multiplexer.begin();
	if (!res) {
		serial_printf("Multiplexer error\n");
		// return;
	}

	for (size_t i = 0; i < N_Beacons; i += 1) pinMode(BEACON_BUS_MAP[i], OUTPUT);
	// for (size_t i = 0; i < N_Beacons; ++i) {
	// 	multiplexer.selectChannel(BEACON_BUS_MAP[i]);
	// 	auto& b = mlx_beacons[i];
	// 	bool h = b.begin_I2C();
	// 	if (h) {
	// 		b.setFilter(mlx90393_filter::MLX90393_FILTER_2);
	// 		b.setOversampling(mlx90393_oversampling::MLX90393_OSR_0);
	// 		b.setResolution(mlx90393_axis::MLX90393_X, mlx90393_resolution::MLX90393_RES_16);
	// 		b.setResolution(mlx90393_axis::MLX90393_Y, mlx90393_resolution::MLX90393_RES_16);
	// 		b.setGain(mlx90393_gain::MLX90393_GAIN_1X);
	// 	} else {
	// 		serial_printf("Beacon %d error\n", (int)i);
	// 	}
	// 	beacon_healthy[i] = h;
	// }
	for (size_t i = 0; i < N_Imus; ++i) {
		select(IMU_BUS_MAP[i]);
		auto& b = imus[i];
		b = GY521(0x68);
		bool h = b.begin();
		if (h) {
			b.setAccelSensitivity(0);
			b.setGyroSensitivity(0);
			b.setThrottle(true);

			b.axe = 0;
			b.aye = 0;
			b.aze = 0;
			b.gxe = 0;
			b.gye = 0;
			b.gze = 0;
		} else {
			serial_printf("IMU %d error\n", (int)i);
		}
		imu_healthy[i] = h;
	}


	for (size_t i = 0; i < N_Beacons; ++i) {
		select(BEACON_BUS_MAP[i]);
		uint8_t revid = 0;
		bool res = beacons[i].revid(revid);
		unselect();

		// RM3100::self_test(BEACON_BUS_MAP[i]);
		beacon_healthy[i] = res;
		if (!res || revid != 0x22) {
			serial_printf("\nError init %d\n", (int)i);
			continue;
		}
		
		serial_printf("Revid of %d = %d\n", (int)i, (int)revid);
		select(BEACON_BUS_MAP[i]);
		beacons[i].set_cycle_count(200, 200, 200);
		unselect();
		select(BEACON_BUS_MAP[i]);
		beacons[i].set_update_rate(HZ_37);
		unselect();
		select(BEACON_BUS_MAP[i]);
		beacons[i].toggle_continuous_mode();
		unselect();
	}

	for (auto& x : last_vectors_idx) x = 0;
}

#define TYPE_MAGNETOMETER 0
#define TYPE_ACCELEROMETER 1
#define TYPE_GYROSCOPE 2

struct __attribute__((packed)) Output {
	uint8_t sync_start[N_Sync_Seq];
	float x;
	float y;
	float z;
	uint8_t id;
	uint8_t type;
	uint32_t time;
	uint8_t sync_end[N_Sync_Seq];
};
#undef round

#define FOR_HUMAN 0
uint32_t tick = 0;
void send_mag(uint8_t id, float x, float y, float z) noexcept {
#if FOR_HUMAN
	serial_printf(
		"Beacon [%d] mag %d (%d %d %d) ut %d %d\n",
		(int)id,
		(int)(sqrtf(x * x + y * y + z * z)),
		(int)(x),
		(int)(y),
		(int)(z),
		(int)tick,
		sizeof(Output)
	);
	// serial_printf(
	// 	"%d\n",
	// 	(int)(100 * sqrtf(x * x + y * y + z * z))
	// );
#else
	Output out;

	for (uint8_t i = 0; i < N_Sync_Seq; ++i) out.sync_start[i] = i;
	for (uint8_t i = 0; i < N_Sync_Seq; ++i) out.sync_end[N_Sync_Seq - i - 1] = i;
	out.time = tick;
	out.type = TYPE_MAGNETOMETER;
	out.id = id;
	out.x = x;
	out.y = y;
	out.z = z;

	Serial.write(reinterpret_cast<uint8_t*>(&out), sizeof(Output));
#endif
	tick += 1;
}
void send_acc(uint8_t id, float x, float y, float z) noexcept {
	Output out;

	for (uint8_t i = 0; i < N_Sync_Seq; ++i) out.sync_start[i] = i;
	for (uint8_t i = 0; i < N_Sync_Seq; ++i) out.sync_end[N_Sync_Seq - i - 1] = i;
	out.time = tick;
	out.type = TYPE_ACCELEROMETER;
	out.id = id;
	out.x = x;
	out.y = y;
	out.z = z;

	Serial.write(reinterpret_cast<uint8_t*>(&out), sizeof(Output));
	tick += 1;
}
void send_gyr(uint8_t id, float x, float y, float z) noexcept {
	Output out;

	for (uint8_t i = 0; i < N_Sync_Seq; ++i) out.sync_start[i] = i;
	for (uint8_t i = 0; i < N_Sync_Seq; ++i) out.sync_end[N_Sync_Seq - i - 1] = i;
	out.time = tick;
	out.type = TYPE_GYROSCOPE;
	out.id = id;
	out.x = x;
	out.y = y;
	out.z = z;

	Serial.write(reinterpret_cast<uint8_t*>(&out), sizeof(Output));
	tick += 1;
}

void loop() {
	for (size_t i = 0; i < N_Beacons; ++i) if (beacon_healthy[i]) {
		select(BEACON_BUS_MAP[i]);
		bool ready = beacons[i].is_data_ready();
		unselect();
		if (!ready) continue;

		Vector read;
		select(BEACON_BUS_MAP[i]);
		beacons[i].read_ut(&read.XAxis, &read.YAxis, &read.ZAxis);
		unselect();
		send_mag((uint8_t)i, read.XAxis * 10, read.YAxis * 10, read.ZAxis * 10);
	}
	for (size_t i = 0; i < N_Imus; ++i) if (imu_healthy[i]) {
		select(IMU_BUS_MAP[i]);

		Vector acc;
		Vector gyr;

		int err = imus[i].read();
		if (err != GY521_OK) {
			serial_printf("Error %d, %d\n", (int)i, (int)err);
			continue;
		}
		acc.XAxis = imus[i].getAccelX();
		acc.YAxis = imus[i].getAccelY();
		acc.ZAxis = imus[i].getAccelZ();
		gyr.XAxis = imus[i].getGyroX();
		gyr.YAxis = imus[i].getGyroY();
		gyr.ZAxis = imus[i].getGyroZ();
		send_acc((uint8_t)i, acc.XAxis, acc.YAxis, acc.ZAxis);
		send_gyr((uint8_t)i, gyr.XAxis, gyr.YAxis, gyr.ZAxis);
		// serial_printf("X %d Y %d Z %d\n", (int)(1000000 * gyr.XAxis), (int)(1000000 * gyr.YAxis), (int)(1000000 * gyr.ZAxis));
	}
}
