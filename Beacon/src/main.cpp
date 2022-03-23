#include "Arduino.h"
#include "Wire.h"
#include "HMC5883L.h"
#include "TCA9548.h"

#include <math.h>

constexpr size_t N_Beacons = 2;
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
bool healthy[N_Beacons] = { false };

TCA9548 multiplexer(0x70);
size_t BUS_MAP[] = {7, 6, 4, 5, 6, 7, 6, 7};


void setup() {
	Wire.setWireTimeout(1000);
	Wire.begin();
	Serial.begin(128000);

	bool res = multiplexer.begin();
	if (!res) {
		serial_printf("Multiplexer error\n");
		return;
	}
	for (size_t i = 0; i < N_Beacons; ++i) {
		multiplexer.selectChannel(BUS_MAP[i]);
		serial_printf("Trying to init %d", (int)i);
		bool res = beacons[i].begin();
		healthy[i] = res;
		if (!res) {
			serial_printf("\nError init %d\n", (int)i);
			continue;
		}
		serial_printf(".");
		
		beacons[i].setSamples(HMC5883L_SAMPLES_8);
		serial_printf(".");
		beacons[i].setRange(HMC5883L_RANGE_4GA);
		serial_printf(".");
		beacons[i].setDataRate(HMC5883L_DATARATE_75HZ);
		serial_printf(".\n");
		serial_printf("Initialized %d\n", (int)i);

	}
}

struct __attribute__((packed)) Output {
	uint8_t sync[N_Sync_Seq];
	float x;
	float y;
	float z;
	uint8_t id;
};
#undef round
void send_mag(uint8_t id, float x, float y, float z) noexcept {
	Output out;

	for (uint8_t i = 0; i < N_Sync_Seq; ++i) out.sync[i] = i;
	out.id = id;
	//out.x = (int16_t)round(cbrt(100000 * (double)x));
	//out.y = (int16_t)round(cbrt(100000 * (double)y));
	//out.z = (int16_t)round(cbrt(100000 * (double)z));
	out.x = x;
	out.y = y;
	out.z = z;

	Serial.write(reinterpret_cast<uint8_t*>(&out), sizeof(Output));
}

Vector last_vectors[N_Beacons];

void loop() {
	for (size_t i = 0; i < N_Beacons; ++i) if (healthy[i]) {
		multiplexer.selectChannel(BUS_MAP[i]);
		auto read = beacons[i].readNormalize();
		// if (memcmp(&last_vectors[i], &read, sizeof(read)) == 0) continue;

		// last_vectors[i] = read;
		send_mag((uint8_t)i, read.XAxis, read.YAxis, read.ZAxis);
	}
}
