#include "Arduino.h"
#include "SparkFun_MAG3110.h"

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

MAG3110 beacons[2];

void setup() {
	Wire.begin();
	Serial.begin(115200);

	for (size_t i = 0; i < 2; ++i) {
		use_bus(i);
		beacons[i].initialize();
		beacons[i].setDR_OS(MAG3110_DR_OS_80_16);
		beacons[i].start();
	}
}

struct __attribute__((packed)) Output {
	uint8_t sync[N_Sync_Seq];
	uint8_t id;
	int32_t x;
	int32_t y;
	int32_t z;
	int32_t sum_error;
};

void send_mag(uint8_t id, int x, int y, int z) noexcept {
	Output out;

	for (uint8_t i = 0; i < N_Sync_Seq; ++i) out.sync[i] = i;
	out.id = id;
	out.x = x;
	out.y = y;
	out.z = z;
	out.sum_error = x + y + z;

	Serial.write(reinterpret_cast<uint8_t*>(&out), sizeof(Output));
}


void loop() {
	for (size_t i = 0; i < 2; ++i) {
		use_bus(i);
		if (!beacons[i].isCalibrated()) {
			if (!beacons[i].isCalibrating()) beacons[i].enterCalMode();
			else                             beacons[i].calibrate();
			continue;
		}
		if (!beacons[i].dataReady()) continue;

		int x = 0, y = 0, z = 0;
		beacons[i].readMag(&x, &y, &z);
		send_mag((uint8_t)i, x, y, z);
	}
}
