#include <vector>
#include <string>
#include <cstdint>
#include <cmath>
#include "Windows.h"

struct Serial_Port_Window {
	HANDLE handle{ INVALID_HANDLE_VALUE };
	bool connected{ false };
	COMSTAT status;
	DWORD errors;
	size_t baud_rate;
};

struct Serial_Port {
	std::string name;

	Serial_Port(std::string name) noexcept;
	~Serial_Port() noexcept;

	Serial_Port(Serial_Port&) = delete;
	Serial_Port& operator=(Serial_Port&) = delete;
	Serial_Port(Serial_Port&&) = default;
	Serial_Port& operator=(Serial_Port&&) = default;

	std::vector<uint8_t> read(size_t n) noexcept;
	std::vector<uint8_t> wait_read(size_t n) noexcept;
	size_t write(const std::vector<uint8_t>& data) noexcept;

private:
	Serial_Port_Window platform;
};

Serial_Port::Serial_Port(std::string name) noexcept : name(std::move(name)) {
	platform.connected = false;

	platform.handle = CreateFileA(
		this->name.c_str(),
		GENERIC_READ | GENERIC_WRITE,
		0,
		NULL,
		OPEN_EXISTING,
		FILE_ATTRIBUTE_NORMAL,
		NULL
	);

	if (platform.handle == INVALID_HANDLE_VALUE) {
		printf(
			"ERROR: Handle was not attached. Reason: %s not available\n",
			this->name.c_str()
		);
		return;
	}

	DCB dcb = { 0 };

	if (!GetCommState(platform.handle, &dcb)) {
		printf(
			"ERROR: GetCommState failed on: %s\n",
			this->name.c_str()
		);
		return;
	}

	platform.baud_rate = CBR_115200;
	dcb.BaudRate = platform.baud_rate;
	dcb.ByteSize = 8;
	dcb.StopBits = ONESTOPBIT;
	dcb.Parity = NOPARITY;
	dcb.fDtrControl = DTR_CONTROL_ENABLE;

	if (!SetCommState(platform.handle, &dcb)) {
		printf(
			"ERROR: SetCommState failed on: %s\n",
			this->name.c_str()
		);
		return;
	}
	
	platform.connected = true;

	PurgeComm(platform.handle, PURGE_RXCLEAR | PURGE_TXCLEAR);
}

Serial_Port::~Serial_Port() noexcept {
	if (platform.connected) {
		platform.connected = false;
		CloseHandle(platform.handle);
	}
}

std::vector<uint8_t> Serial_Port::read(size_t n) noexcept {
	DWORD bytes_read{ 0 };

	ClearCommError(platform.handle, &platform.errors, &platform.status);

	if ((size_t)platform.status.cbInQue < n)
		n = (size_t)platform.status.cbInQue;

	std::vector<uint8_t> result(n);
	ReadFile(platform.handle, result.data(), n, &bytes_read, NULL);
	result.resize(bytes_read);

	return result;
}

std::vector<uint8_t> Serial_Port::wait_read(size_t n) noexcept {
	size_t bytes_read{ 0 };

	std::vector<uint8_t> result(n);
	while (bytes_read < n) {
		ClearCommError(platform.handle, &platform.errors, &platform.status);

		DWORD byte_read_this_time{ 0 };


		auto to_read = 
			((size_t)platform.status.cbInQue < n - bytes_read)
				? (size_t)platform.status.cbInQue
				: (n - bytes_read);
		ReadFile(
			platform.handle, result.data() + bytes_read, to_read, &byte_read_this_time, nullptr
		);
		bytes_read += byte_read_this_time;

		Sleep(100);
	}
	result.resize(bytes_read);
	return result;
}

size_t Serial_Port::write(const std::vector<uint8_t>& data) noexcept {
	DWORD send{ 0 };

	if (!WriteFile(platform.handle, (void*)data.data(), data.size(), &send, 0))
		ClearCommError(platform.handle, &platform.errors, &platform.status);

	return send;
}

constexpr auto N_Sync_Seq = 32;

#pragma pack(1)
struct Inputs_DRV425 {
	std::uint8_t sync[N_Sync_Seq];
	std::uint32_t digital_1;
	std::uint32_t digital_2;
};

#pragma pack(1)
struct Inputs_MAG3110 {
	std::uint8_t sync[N_Sync_Seq];
	std::uint8_t id;
	float x;
	float y;
	float z;
};

struct Opts {
	size_t read_at_a_time = sizeof(Inputs_MAG3110) * 32;
	std::string com_name = "\\\\.\\COM3";

	static Opts parse(int argc, char** argv) noexcept {
		Opts res;

		auto eq = [](auto a, auto b) {
			return strcmp(a, b) == 0;
		};

		for (int i = 0; i < argc; ++i) {
			auto it = argv[i];

			if (eq(it, "-r")) res.read_at_a_time = std::stoi(argv[i + 1]);
			if (eq(it, "-c")) res.com_name = argv[i + 1];
		}

		return res;
	}
};

double lsb_to_volt(std::uint32_t lsb) noexcept {
	return lsb * 5 / 1023.0;
}

double volt_to_tesla(double volt) noexcept {
	return volt / (4 * 12.2 * 1000);
}

double mmag_to_tesla(float mag) noexcept {
	return (1 / 1'000'000.0) * mag / 1'000.0;
}

int32_t read_int32(const std::vector<std::uint8_t>& buf, size_t i) noexcept {
	return *reinterpret_cast<const int32_t*>(buf.data() + i);
}

#pragma pack(1)
struct Reading {
	double beacon1;
	double beacon2;
	bool pressed;
};
constexpr const char* Mail_Name = "\\\\.\\Mailslot\\SP";

HANDLE open_slot() noexcept {
	auto mail_slot = CreateFileA(
		Mail_Name,
		GENERIC_WRITE,
		FILE_SHARE_READ,
		nullptr,
		OPEN_EXISTING,
		FILE_ATTRIBUTE_NORMAL,
		nullptr
	);
	if (mail_slot == INVALID_HANDLE_VALUE) {
		fprintf(stderr, "Error while opening mail slot\n");
		fprintf(stderr, "CreateFile failed with %d.\n", GetLastError());
	}

	return mail_slot;
}

void send_over_mail(
	HANDLE mail_slot, Reading readings
) noexcept {
	DWORD cbWritten;

	auto fResult = WriteFile(
		mail_slot,
		reinterpret_cast<const char*>(&readings),
		sizeof(readings),
		&cbWritten,
		nullptr
	);

	if (!fResult) {
		// fprintf(stderr, "Write failed with %d.\n", GetLastError());
	}
}
constexpr auto Reading_Byte_Size =
#if 0 // DRV425
	sizeof(Inputs_DRV425);
#else // MAG3110
	sizeof(Inputs_MAG3110);
#endif

int main(int argc, char** argv) {
	Opts opts = Opts::parse(argc - 1, argv + 1);

	Serial_Port serial(opts.com_name);
	HANDLE mail_slot = open_slot();
	std::vector<Reading> readings;

	uint8_t Sync_Seq[N_Sync_Seq];
	for (uint8_t i = 0; i < N_Sync_Seq; ++i) Sync_Seq[i] = i;

	while (true) {
		auto vec = serial.wait_read(opts.read_at_a_time);

		readings.clear();
		size_t offset{ 0 };

		for (; offset + Reading_Byte_Size < vec.size(); ++offset)
			if (memcmp(vec.data() + offset, Sync_Seq, N_Sync_Seq) == 0) break;

		if (offset > 0) {
			printf("Ignored %.*s\n", offset, vec.data());
		}

		Reading r;
		bool reading1_ready = false;
		bool reading2_ready = false;
		for (size_t i = offset; i < vec.size(); i += sizeof(Inputs_MAG3110)) {
#if 0 // DRV425
			Inputs_DRV425 in = *reinterpret_cast<Inputs*>(vec.data() + i);
			printf(
				"Read % 5u LSB % 5u LSB % 10.9llf V % 5llf V % 5llf uT % 10.9llf uT\n",
				in.digital_1,
				in.digital_2,
				lsb_to_volt(in.digital_1),
				lsb_to_volt(in.digital_2),
				1'000'000 * volt_to_tesla(lsb_to_volt(in.digital_1)),
				1'000'000 * volt_to_tesla(lsb_to_volt(in.digital_2))
			);
#else // MAG3110
			Inputs_MAG3110 in =
				*reinterpret_cast<Inputs_MAG3110*>(vec.data() + i);

			double tx = mmag_to_tesla(in.x);
			double ty = mmag_to_tesla(in.y);
			double tz = mmag_to_tesla(in.z);
			double t = std::hypot(tx, ty, tz);

			if (in.id == 0) {
				r.beacon1 = t;
				reading1_ready = true;
			}
			if (in.id == 1) {
				r.beacon2 = t;
				reading2_ready = true;
			}
			if (reading1_ready && reading2_ready) {
				r.pressed = true;
				readings.push_back(r);
				reading1_ready = false;
				reading2_ready = false;
			}
			printf(
				"Read[%d] % 10.4f MAG(x) % 10.4f MAG(y) % 10.4lf MAG(z)\n",
				(int)in.id,
				in.x,
				in.y,
				in.z
			);

#endif
		}
		for (auto& x : readings) send_over_mail(mail_slot, x);
	}

	return 0;
}