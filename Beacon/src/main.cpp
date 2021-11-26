#include "Arduino.h"
#include "Wire.h"
#include "HMC5883L.h"

constexpr auto N_Sync_Seq = 32;

void serial_printf(const char *fmt, ...) {
	va_list va;
	va_start(va, fmt);
	char buf[vsnprintf(NULL, 0, fmt, va) + 1];
	vsprintf(buf, fmt, va);
	Serial.print(buf);
	va_end(va);
}

void use_bus(uint8_t bus) noexcept {
	Wire.beginTransmission(0x70);  // TCA9548A address is 0x70
	Wire.write(1 << bus);          // send byte to select bus
	Wire.endTransmission();
}

constexpr size_t N_Beacons = 2;
HMC5883L beacons[N_Beacons];
bool healthy = false;
size_t BUS_MAP[] = {2, 3};

void setup() {
	Wire.begin();
	Serial.begin(115200);

	for (size_t i = 0; i < N_Beacons; ++i) {
		use_bus(BUS_MAP[i]);
		bool res = beacons[i].begin();
		beacons[i].setSamples(HMC5883L_SAMPLES_8);
		beacons[i].setRange(HMC5883L_RANGE_0_88GA);
		beacons[i].setDataRate(HMC5883L_DATARATE_75HZ);
		if (!res) {
			serial_printf("Error init %d\n", (int)i);
			return;
		}
	}

	healthy = true;
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
	out.x = x

	Serial.write(reinterpret_cast<uint8_t*>(&out), sizeof(Output));
}

void loop() {
	for (size_t i = 0; i < N_Beacons; ++i) {
		use_bus(BUS_MAP[i]);
		auto read = beacons[i].readNormalize();
		send_mag((uint8_t)i, read.XAxis, read.YAxis, read.ZAxis);
	}
}
