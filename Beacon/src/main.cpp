#include "Arduino.h"
#include "Wire.h"
#include "HMC5883L.h"
#include "TCA9548.h"

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


constexpr auto N_Sync_Seq = 32;

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

constexpr size_t N_Beacons = 6;
HMC5883L beacons[N_Beacons];
bool healthy[N_Beacons] = { false };

TCA9548 multiplexer(0x70);
size_t BUS_MAP[] = {2, 3, 4, 5, 6, 7, 6, 7};


void setup() {
	Wire.setWireTimeout(1000);
	Wire.begin();
	Serial.begin(115200);

	HMC5883L beacon_test;
	if (!beacon_test.begin()) {
		serial_printf("Beacon test FAIL !");
	} else {
		serial_printf("Beacon test OK !");
	}

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
		beacons[i].setSamples(HMC5883L_SAMPLES_1);
		serial_printf(".");
		beacons[i].setRange(HMC5883L_RANGE_8_1GA);
		serial_printf(".");
		beacons[i].setDataRate(HMC5883L_DATARATE_75HZ);
		serial_printf(".\n");
		serial_printf("Initialized %d\n", (int)i);

	}
}

struct __attribute__((packed)) Output {

	uint8_t sync[N_Sync_Seq];
	uint8_t id;
	float x;
	float y;
	float z;
};

void send_mag(uint8_t id, float x, float y, float z) noexcept {
	Output out;

	for (uint8_t i = 0; i < N_Sync_Seq; ++i) out.sync[i] = i;
	out.id = id;
	out.x = x;
	out.y = y;
	out.z = z;

	Serial.write(reinterpret_cast<uint8_t*>(&out), sizeof(Output));
}

void loop() {
	for (size_t i = 0; i < N_Beacons; ++i) if (healthy[i]) {
		multiplexer.selectChannel(BUS_MAP[i]);
		auto read = beacons[i].readNormalize();
		send_mag((uint8_t)i, read.XAxis, read.YAxis, read.ZAxis);
	}
}
