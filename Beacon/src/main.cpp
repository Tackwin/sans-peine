#include "Arduino.h"
#include "Wire.h"
#include "HMC5883L.h"
#include "TCA9548.h"
#include "MLX90393.h"
#include "GY521.h"

#include <math.h>

constexpr size_t N_Beacons = 4;
constexpr size_t N_Imus = 2;
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

HMC5883L beacons[N_Beacons];
Adafruit_MLX90393 mlx_beacons[N_Beacons];
bool beacon_healthy[N_Beacons] = { false };

TCA9548 multiplexer(0x70);
size_t BEACON_BUS_MAP[] = {0, 1, 7, 6, 6, 7, 6, 7};

GY521 imus[N_Imus];
size_t IMU_BUS_MAP[] = {2, 3, 3, 3, 3, 3, 3, 3};
bool imu_healthy[N_Imus] = { false };

constexpr size_t ROLLING = 1;
int last_vectors_idx[N_Beacons];
Vectori16 last_vectors[N_Beacons * ROLLING];

void setup() {
	// Serial.begin(128000);
	Serial.begin(115200);
	Wire.setWireTimeout(0);
	Wire.begin();
	bool res = multiplexer.begin();
	if (!res) {
		serial_printf("Multiplexer error\n");
		// return;
	}
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
		multiplexer.selectChannel(IMU_BUS_MAP[i]);
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
		multiplexer.selectChannel(BEACON_BUS_MAP[i]);
		bool res = beacons[i].begin();
		beacon_healthy[i] = res;
		if (!res) {
			serial_printf("\nError init %d\n", (int)i);
			continue;
		}
		
		beacons[i].setSamples(HMC5883L_SAMPLES_1);
		beacons[i].setRange(HMC5883L_RANGE_1_3GA);
		beacons[i].setDataRate(HMC5883L_DATARATE_30HZ);
		// beacons[i].setMeasurementMode(HMC5883L_SINGLE);
	}

	for (auto& x : last_vectors_idx) x = 0;
}

#define TYPE_MAGNETOMETER 0
#define TYPE_ACCELEROMETER 1
#define TYPE_GYROSCOPE 2

struct __attribute__((packed)) Output {
	uint8_t sync[N_Sync_Seq];
	float x;
	float y;
	float z;
	uint8_t id;
	uint8_t type;
};
#undef round
void send_mag(uint8_t id, float x, float y, float z) noexcept {
	Output out;

	for (uint8_t i = 0; i < N_Sync_Seq; ++i) out.sync[i] = i;
	out.type = TYPE_MAGNETOMETER;
	out.id = id;
	out.x = x;
	out.y = y;
	out.z = z;

	Serial.write(reinterpret_cast<uint8_t*>(&out), sizeof(Output));
}
void send_acc(uint8_t id, float x, float y, float z) noexcept {
	Output out;

	for (uint8_t i = 0; i < N_Sync_Seq; ++i) out.sync[i] = i;
	out.type = TYPE_ACCELEROMETER;
	out.id = id;
	out.x = x;
	out.y = y;
	out.z = z;

	Serial.write(reinterpret_cast<uint8_t*>(&out), sizeof(Output));
}
void send_gyr(uint8_t id, float x, float y, float z) noexcept {
	Output out;

	for (uint8_t i = 0; i < N_Sync_Seq; ++i) out.sync[i] = i;
	out.type = TYPE_GYROSCOPE;
	out.id = id;
	out.x = x;
	out.y = y;
	out.z = z;

	Serial.write(reinterpret_cast<uint8_t*>(&out), sizeof(Output));
}

void loop() {
	for (size_t i = 0; i < N_Beacons; ++i) if (beacon_healthy[i]) {
		continue;
		multiplexer.selectChannel(BEACON_BUS_MAP[i]);

		auto read = beacons[i].readRawi16();
		send_mag(
			(uint8_t)i,
			read.XAxis * beacons[i].mgPerDigit,
			read.YAxis * beacons[i].mgPerDigit,
			read.ZAxis * beacons[i].mgPerDigit
		);
		
	}
	for (size_t i = 0; i < N_Imus; ++i) if (imu_healthy[i]) {
		multiplexer.selectChannel(IMU_BUS_MAP[i]);

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
		// send_acc((uint8_t)i, acc.XAxis, acc.YAxis, acc.ZAxis);
		// send_gyr((uint8_t)i, gyr.XAxis, gyr.YAxis, gyr.ZAxis);
		serial_printf("X %d Y %d Z %d\n", (int)(1000000 * gyr.XAxis), (int)(1000000 * gyr.YAxis), (int)(1000000 * gyr.ZAxis));
	}
}
