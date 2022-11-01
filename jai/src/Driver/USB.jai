#scope_file

#import "Basic";

#if OS == .WINDOWS {
	Windows :: #import "Windows";
	#import "Windows_Utf8";
	BYTE :: u8;

	DCB :: struct {
		DCBlength: Windows.DWORD;
		BaudRate:  Windows.DWORD;
		Bitfield:  u32;
		reserved:  Windows.WORD;
		XonLim:    Windows.WORD;
		XoffLim:   Windows.WORD;
		ByteSize:  BYTE;
		Parity:    BYTE;
		StopBits:  BYTE;
		XonChar:   u8;
		XoffChar:  u8;
		ErrorChar: u8;
		EofChar:   u8;
		EvtChar:   u8;
		reserved2: Windows.WORD;
	};

	COMSTAT :: struct {
		Bitfields: u32;
		cbInQue:   Windows.DWORD;
		cbOutQue:  Windows.DWORD;
	};
	kernel32    :: #foreign_system_library "kernel32";

	GetCommState :: (handle: Windows.HANDLE, dcb: *DCB) -> Windows.BOOL #foreign kernel32;
	SetCommState :: (handle: Windows.HANDLE, dcb: *DCB) -> Windows.BOOL #foreign kernel32;
	PurgeComm :: (handle: Windows.HANDLE, flags: Windows.DWORD) -> Windows.BOOL #foreign kernel32;
	ClearCommError :: (
		handle: Windows.HANDLE, errors: *Windows.DWORD, comstat: *COMSTAT
	) -> Windows.BOOL #foreign kernel32;

	Native_Handle :: Windows.HANDLE;

	make_windows :: (serial: *Serial) {
		namew, success, namew_len := utf8_to_wide(serial.name);
		serial.handle = Windows.CreateFileW(
			namew,
			Windows.GENERIC_READ | Windows.GENERIC_WRITE,
			0,
			null,
			Windows.OPEN_EXISTING,
			Windows.FILE_ATTRIBUTE_NORMAL,
			null
		);

		dcb: DCB;
		if !GetCommState(serial.handle, *dcb) {
			log_error("Error: Handle was not attached. Reason: % not available\n", serial.name);
			serial.connected = false;
			return;
		}

		ONESTOPBIT :: 0;
		NOPARITY   :: 0;
		DTR_CONTROL_ENABLE :: 0x10;
		dcb.BaudRate = 128_000;
		dcb.ByteSize = 8;
		dcb.StopBits = ONESTOPBIT;
		dcb.Parity = NOPARITY;
		dcb.Bitfield = DTR_CONTROL_ENABLE;

		if !SetCommState(serial.handle, *dcb) {
			log_error("Error: SetCommState failed on %\n", serial.name);
			serial.connected = false;
			return;
		}

		PURGE_RXCLEAR :: 0x8;
		PURGE_TXCLEAR :: 0x4;
		PurgeComm(serial.handle, PURGE_RXCLEAR | PURGE_TXCLEAR);
		serial.connected = true;
	}

	close_windows :: (serial: *Serial) {
		if !serial.connected return;
		serial.connected = false;
		Windows.CloseHandle(handle);
	}

	instant_read_windows :: (serial: Serial, buffer: [] u8) -> int {
		read := 0;
		errors: Windows.DWORD = 0;
		stat: COMSTAT;
		ClearCommError(serial.handle, *errors, *stat);

		Windows.ReadFile(serial.handle, buffer.data, xx buffer.count, xx *read, null);
		return read;
	};

	write_windows :: (serial: Serial, buffer: [] u8) -> bool {
		sent : u32 = 0;
		if !Windows.WriteFile(serial.handle, buffer.data, xx buffer.count, *sent, xx 0)
		{
			ClearCommError(serial.handle, null, null);
			return false;
		}
		return sent == buffer.count;
	}
}

#scope_export

Serial :: struct {
	name: string;
	connected: bool;

	handle: Native_Handle;

	make :: (name: string) -> Serial {
		serial: Serial = ---;
		serial.name = name;
		#if OS == .WINDOWS make_windows(*serial);
		return serial;
	}
}

close :: (serial: *Serial) {
	#if OS == .WINDOWS close_windows(serial);
}

instant_read :: (serial: Serial, buffer: [] u8) -> int {
	#if OS == .WINDOWS return instant_read_windows(serial, buffer);
}

read_to :: (serial: Serial, buffer: [] u8) {
	loader := buffer;

	while loader.count > 0 {
		read := instant_read(serial, loader);
		loader.data += read;
		loader.count -= read;

		sleep_milliseconds(10);
	}
}

write :: (serial: Serial, buffer: [] u8) -> bool {
	#if OS == .WINDOWS return write_windows(serial, buffer);
	return false;
}