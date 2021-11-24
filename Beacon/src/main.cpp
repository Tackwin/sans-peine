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
// size_t BUS_MAP[] = {7, 6};

void setup() {
	Wire.begin();
	Serial.begin(115200);

	for (size_t i = 0; i < 1; ++i) {
		// use_bus(BUS_MAP[i]);
		bool res = beacons[i].initialize();
		beacons[i].setDR_OS(MAG3110_DR_OS_80_16);
		beacons[i].start();
	}
}

struct __attribute__((packed)) Output {
	uint8_t sync[N_Sync_Seq];
	uint8_t id;
	int32_t x0;
	int32_t x1;
	int32_t x2;
	int32_t y0;
	int32_t y1;
	int32_t y2;
	int32_t z0;
	int32_t z1;
	int32_t z2;
	int32_t sum_error;
};

void send_mag(uint8_t id, int x, int y, int z) noexcept {
	Output out;

	for (uint8_t i = 0; i < N_Sync_Seq; ++i) out.sync[i] = i;
	out.id = id;
	out.x0 = x;
	out.x1 = x;
	out.x2 = x;
	out.y0 = y;
	out.y1 = y;
	out.y2 = y;
	out.z0 = z;
	out.z1 = z;
	out.z2 = z;
	out.sum_error = x + y + z;

	Serial.write(reinterpret_cast<uint8_t*>(&out), sizeof(Output));
}

size_t cal_samples[2] = {0, 0};

void loop() {
	for (size_t i = 0; i < 1; ++i) {
		// use_bus(BUS_MAP[i]);
		if (!beacons[i].isCalibrated()) {
			if (!beacons[i].isCalibrating()) {
				beacons[i].enterCalMode();
			}
			else {
				beacons[i].calibrate();
				cal_samples[i]++;
			}


			if (cal_samples[i] > 100000) {
				beacons[i].exitCalMode();
			}
			continue;
		}
		if (!beacons[i].dataReady()) {
			// serial_printf("%d NOT READY\n", (int)i);
			continue;
		}
		int x = 0, y = 0, z = 0;
		beacons[i].readMag(&x, &y, &z);
		send_mag((uint8_t)i, x, y, z);
	}
}
