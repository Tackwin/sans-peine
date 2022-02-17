


#include <vector>
#include <string>
#include <cstdint>
#include <cmath>

#define NOMINMAX
#include "Windows.h"

#include "Definitions.hpp"
#include "IPC.hpp"
#include "Macro.hpp"

template<typename T>
T min(T a, T b) noexcept { return a > b ? a : b; }

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
	void wait_read(std::vector<uint8_t>& buffer) noexcept;
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

	platform.baud_rate = CBR_128000;
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
	result.resize(min((size_t)bytes_read, n));

	return result;
}

std::vector<uint8_t> Serial_Port::wait_read(size_t n) noexcept {
	std::vector<uint8_t> result(n);
	wait_read(result);
	return result;
}
void Serial_Port::wait_read(std::vector<uint8_t>& buffer) noexcept {
	size_t bytes_read{ 0 };

	while (bytes_read < buffer.size()) {
		ClearCommError(platform.handle, &platform.errors, &platform.status);

		DWORD byte_read_this_time{ 0 };

		auto to_read = 
			((size_t)platform.status.cbInQue < buffer.size() - bytes_read)
				? (size_t)platform.status.cbInQue
				: (buffer.size() - bytes_read);
		ReadFile(
			platform.handle, buffer.data() + bytes_read, to_read, &byte_read_this_time, nullptr
		);
		bytes_read += byte_read_this_time;

		Sleep(10);
	}
	buffer.resize(min(bytes_read, buffer.size()));
}

size_t Serial_Port::write(const std::vector<uint8_t>& data) noexcept {
	DWORD send{ 0 };

	if (!WriteFile(platform.handle, (void*)data.data(), data.size(), &send, 0))
		ClearCommError(platform.handle, &platform.errors, &platform.status);

	return send;
}


#pragma pack(1)
struct Inputs_DRV425 {
	std::uint8_t sync[N_Sync_Seq];
	std::uint32_t digital_1;
	std::uint32_t digital_2;
};

#pragma pack(1)
struct Inputs_MAG3110 {
	std::uint8_t sync[N_Sync_Seq];
	float x;
	float y;
	float z;
	std::uint8_t id;
};

struct Opts {
	std::string com_name = "\\\\.\\COM3";

	const char* record = nullptr;
	const char* replay = nullptr;

	static Opts parse(int argc, char** argv) noexcept {
		Opts res;

		auto eq = [](auto a, auto b) {
			return strcmp(a, b) == 0;
		};

		for (int i = 0; i < argc; ++i) {
			auto it = argv[i];

			if (eq(it, "-c")) res.com_name = argv[++i];
			if (eq(it, "--record")) res.record = argv[++i];
			if (eq(it, "--replay")) res.replay = argv[++i];
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
	return (1 / 1'000'000.0) * mag / 10.0;
}

int32_t read_int32(const std::vector<std::uint8_t>& buf, size_t i) noexcept {
	return *reinterpret_cast<const int32_t*>(buf.data() + i);
}

void append_to(const char* filename, void* data, size_t n_data) {
	auto f = fopen(filename, "ab");
	if (!f) return;
	fseek(f, 0, SEEK_END);
	fwrite(data, 1, n_data, f);
	fclose(f);
}

std::vector<Reading> replay(const char* filename) {
	auto f = fopen(filename, "rb");
	if (!f) return {};

	fseek(f, 0, SEEK_END);
	auto size = ftell(f);
	fseek(f, 0, SEEK_SET);

	std::vector<Reading> result;
	result.resize(size / sizeof(Reading));

	fread(result.data(), 1, size, f);

	fclose(f);
	return result;
}

constexpr auto Reading_Byte_Size =
#if 0 // DRV425
	sizeof(Inputs_DRV425);
#else // MAG3110
	sizeof(Inputs_MAG3110);
#endif

#include <thread>
#include <atomic>

struct State {
	Serial_Port port = Serial_Port("\\\\.\\COM3");

	size_t ring_buffer_size = 0;
	std::atomic<size_t> first = 0;
	std::atomic<size_t> final = 0;
	Reading* ring_buffer = nullptr;

	std::atomic<bool> critical_section = false;

	std::atomic<bool> request_stop = false;
	std::atomic<bool> is_stop = false;
};
bool False = false;

__declspec(dllexport) extern "C" void stop(void* ptr) {
	State* state = (State*)ptr;
	state->request_stop = true;
	while(!state->is_stop);
}

__declspec(dllexport) extern "C" void play(void* ptr) {
	thread_local std::vector<uint8_t> buffer;
	size_t S = sizeof(Inputs_MAG3110);
	uint8_t Sync_Seq[N_Sync_Seq];
	
	State* state = (State*)ptr;

	for (uint8_t i = 0; i < N_Sync_Seq; ++i) Sync_Seq[i] = i;
	buffer.resize(32 * sizeof(Inputs_MAG3110));

	while (!state->request_stop) {
		state->port.wait_read(buffer);

		size_t offset{ 0 };

		for (; offset + Reading_Byte_Size < buffer.size(); ++offset)
			if (memcmp(buffer.data() + offset, Sync_Seq, N_Sync_Seq) == 0) break;

		bool reading_ready[N_Beacons] = { false };
		Reading r;
		r.pressed = true;

		for (size_t i = offset; i + S < buffer.size(); i += S) {
			Inputs_MAG3110 in = *reinterpret_cast<Inputs_MAG3110*>(buffer.data() + i);

			if (in.id >= N_Beacons) continue;

			double tx = mmag_to_tesla(in.x);
			double ty = mmag_to_tesla(in.y);
			double tz = mmag_to_tesla(in.z);
			double t = std::hypot(tx, ty, tz);

			reading_ready[in.id] = true;
			r.beacons[in.id].x = tx;
			r.beacons[in.id].y = ty;
			r.beacons[in.id].z = tz;

			if (memchr(reading_ready, 0, N_Beacons) != NULL) continue;
			r.timestamp = seconds();

			while (!std::atomic_compare_exchange_strong(&state->critical_section, &False, true));
			defer { state->critical_section = false; };

			size_t idx = state->final++;
			size_t size = state->ring_buffer_size;
			
			state->ring_buffer[idx % size] = r;

			if ((state->final % size) <= (state->first % size)) state->first = state->final + 1;

			memset(reading_ready, 0, N_Beacons);
		}
	}
	state->is_stop = true;
}

__declspec(dllexport) extern "C" bool next(void* ptr, Reading* read) {
	State* state = (State*)ptr;

	while (!std::atomic_compare_exchange_strong(&state->critical_section, &False, true));
	defer { state->critical_section = false; };

	if (state->first == state->final) return false;

	size_t idx = state->first % state->ring_buffer_size;
	*read = state->ring_buffer[idx];
	state->first++;

	return true;
}

__declspec(dllexport) extern "C" void* init() {
	State* state = new State;
	state->ring_buffer_size = 1024 * 1024;
	state->ring_buffer = new Reading[state->ring_buffer_size];
	return state;
}

__declspec(dllexport) extern "C" void shut(void* ptr) {
	stop(ptr);
	State* state = (State*)ptr;
	delete[] state->ring_buffer;
	delete state;
}

extern "C" void read_n_samples(void* ptr, Reading* out_readings, u32 n) {
	thread_local std::vector<uint8_t> buffer;
	size_t S = sizeof(Inputs_MAG3110);
	uint8_t Sync_Seq[N_Sync_Seq];
	
	Serial_Port& port = *(Serial_Port*)ptr;

	for (uint8_t i = 0; i < N_Sync_Seq; ++i) Sync_Seq[i] = i;
	buffer.resize(256 * sizeof(Inputs_MAG3110));

	while (n > 0) {
		port.wait_read(buffer);

		size_t offset{ 0 };

		for (; offset + Reading_Byte_Size < buffer.size(); ++offset)
			if (memcmp(buffer.data() + offset, Sync_Seq, N_Sync_Seq) == 0) break;

		bool reading_ready[N_Beacons] = { false };
		Reading r;

		for (size_t i = offset; i + S < buffer.size(); i += S) {
			Inputs_MAG3110 in = *reinterpret_cast<Inputs_MAG3110*>(buffer.data() + i);

			if (in.id >= N_Beacons) continue;

			double tx = mmag_to_tesla(in.x);
			double ty = mmag_to_tesla(in.y);
			double tz = mmag_to_tesla(in.z);
			double t = std::hypot(tx, ty, tz);

			reading_ready[in.id] = true;
			r.beacons[in.id].x = tx;
			r.beacons[in.id].y = ty;
			r.beacons[in.id].z = tz;

			if (memchr(reading_ready, 0, N_Beacons) == NULL) {
				r.pressed = true;
				r.timestamp = seconds();
				*out_readings = r;
				++out_readings;
				--n;
				memset(reading_ready, 0, N_Beacons);
			}
		}
	}
}

int main(int argc, char** argv) {
	Opts opts = Opts::parse(argc - 1, argv + 1);
	auto mail_slot = IPC::open_slot(Mail_Name);
	if (!mail_slot) mail_slot = IPC::make_slot(Mail_Name);
	if (!mail_slot) {
		printf("Couldn't open mail slot :'(\n");
		return -1;
	}

	if (opts.replay) {
		auto read = replay(opts.replay);
		for (auto& x : read) IPC::write(mail_slot, &x, sizeof(x));

		return 0;
	}

	Serial_Port serial(opts.com_name);
	std::vector<Reading> readings;

	uint8_t Sync_Seq[N_Sync_Seq];
	for (uint8_t i = 0; i < N_Sync_Seq; ++i) Sync_Seq[i] = i;

	while (true) {
		readings.clear();
		readings.resize(10);
		read_n_samples(&serial, readings.data(), readings.size());
		
		for (auto& x : readings) {
			for (size_t i = 0; N_Beacons; ++i) {
				double t = std::hypot(x.beacons[i].x, x.beacons[i].y, x.beacons[i].z);
				auto format = "Read[%zu] % 10.5f IN(x) % 10.5f IN(y) % 10.5f IN(z) % 10.9lf\n";
				printf(format, i, x.beacons[i].x, x.beacons[i].y, x.beacons[i].z, t);
			}
			IPC::write(mail_slot, &x, sizeof(x));
		}

		if (opts.record)
			append_to(opts.record, readings.data(), readings.size() * sizeof(Reading));

		auto vec = serial.wait_read(256 * sizeof(Inputs_MAG3110));

		readings.clear();
		size_t offset{ 0 };

		for (; offset + Reading_Byte_Size < vec.size(); ++offset)
			if (memcmp(vec.data() + offset, Sync_Seq, N_Sync_Seq) == 0) break;

		if (offset > 0) {
			printf("Ignored %.*s\n", (int)offset, vec.data());
		}

		Reading r;
		bool reading_ready[N_Beacons] = { false };
		for (
			size_t i = offset;
			i + sizeof(Inputs_MAG3110) < vec.size();
			i += sizeof(Inputs_MAG3110)
		) {
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

			if (in.id >= N_Beacons) {
				printf("Input corrupted. Skip.\n");
				continue;
			}

			double tx = mmag_to_tesla(in.x);
			double ty = mmag_to_tesla(in.y);
			double tz = mmag_to_tesla(in.z);
			double t = std::hypot(tx, ty, tz);

			reading_ready[in.id] = true;
			r.beacons[in.id].x = tx;
			r.beacons[in.id].y = ty;
			r.beacons[in.id].z = tz;

			if (memchr(reading_ready, 0, N_Beacons) == NULL) {
				r.pressed = true;
				r.timestamp = seconds();
				readings.push_back(r);
				memset(reading_ready, 0, N_Beacons);
			}

			printf(
				"Read[%d] % 10.5f IN(x) % 10.5f IN(y) % 10.5f IN(z) % 10.9lf\n",
				(int)in.id,
				in.x,
				in.y,
				in.z,
				t
			);

#endif
		}
		for (auto& x : readings) IPC::write(mail_slot, &x, sizeof(x));

		if (opts.record)
			append_to(opts.record, readings.data(), readings.size() * sizeof(Reading));
	}

	return 0;
}
