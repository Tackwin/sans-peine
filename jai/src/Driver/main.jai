#scope_file
#import "Basic";

USB :: #import "USB";
#import "Reading";
Interface :: #import "Interface";
#import "Thread";
#import "Physics";

#if OS == .WINDOWS {
	Windows :: #import "Windows";
	#import "Windows_Utf8";
	SetThreadDescription :: (thread: Windows.HANDLE, name: Windows.LPCWSTR) #foreign kernel32;

	kernel32    :: #foreign_system_library "kernel32";
}

INTERFACE_STRING_NO_RETURN :: #string DONE
	#program_export "%1"
	_%1 :: (memory: *u8) {
		state: *State = cast(*State) memory;
		%1(state);
	}
DONE;
INTERFACE_STRING :: #string DONE
	#program_export "%1"
	_%1 :: (memory: *u8) -> %2 {
		state: *State = cast(*State) memory;
		return %1(state);
	}
DONE;
INTERFACE_STRING2 :: #string DONE
	#program_export "%1"
	_%1 :: (memory: *u8) -> %2 {
		state: *State = cast(*State) memory;
		a, b := %1(state);
		return a, b;
	}
DONE;

make_interface :: () -> string {
	builder: String_Builder;
	defer free_buffers(*builder);

	print_to_builder(*builder, INTERFACE_STRING, "size", "int");
	print_to_builder(*builder, INTERFACE_STRING, "init", "Interface.Init_Result");
	print_to_builder(*builder, INTERFACE_STRING2, "next", "Interface.Input, bool");
	for string.["shut", "play", "stop"]
		print_to_builder(*builder, INTERFACE_STRING_NO_RETURN, it);

	return builder_to_string(*builder);
}

#insert #run make_interface();

#scope_export

State :: struct {
	serial: USB.Serial;
	thread: *Thread = null;
	mutex: Mutex;

	thread_launched := false;
	running := false;

	INPUT_RING_BUFFER_SIZE : u64 : 256;
	inputs : [INPUT_RING_BUFFER_SIZE] Packet;
	beg_index : u64 = 0;
	end_index : u64 = 0;
}

size :: (state: *State) -> int {
	return size_of(State);
}

init :: (state: *State) -> Interface.Init_Result {
	<<state = .{};
	result: Interface.Init_Result;
	state.serial = USB.Serial.make("\\\\.\\COM3");
	result.success = state.serial.connected;
	return result;
}

shut :: (state: *State) {
	state.running = false;
	thread_destroy(state.thread);
	destroy(*state.mutex);
}

play :: (state: *State) {
	if !state.thread {
		state.thread = thread_create(background);
		init(*state.mutex, "Mutex Driver Sans Peine");
		state.thread.data = cast(*void) state;
		thread_start(state.thread);
	}
	state.running = true;
}

stop :: (state: *State) {
	state.running = false;
}

next :: (state: *State) -> Interface.Input, bool {
	if state.beg_index >= state.end_index return .{}, false;

	packet : Packet = ---;
	{
		lock(*state.mutex);
		defer unlock(*state.mutex);

		if state.beg_index >= state.end_index return .{}, false;
		packet = state.inputs[state.beg_index % State.INPUT_RING_BUFFER_SIZE];
		state.beg_index += 1;
	}

	out: Interface.Input;
	out.id = packet.id;
	out.kind = xx packet.kind;
	out.time = packet.time;
	if packet.kind == {
		case .Magnetometer;  out.mag = packet.mag;
		case .Accelerometer; out.acc = packet.acc;
		case .Gyroscope;     out.gyr = packet.gyr;
	}
	return out, true;
}

background :: (thread: *Thread) -> s64 {
	#if OS == .WINDOWS {
		SetThreadDescription(
			Windows.GetCurrentThread(),
			utf8_to_wide("Background usb driver Sans peine")
		);
	}

	state: *State = cast(*State) thread.data;
	while state.running {
		packet, has_packet := read_next(state.serial, 8);
		if !has_packet continue;

		lock(*state.mutex);
		defer unlock(*state.mutex);
		
		for 0..packet.count-1
			state.inputs[(state.end_index + xx it) % State.INPUT_RING_BUFFER_SIZE] = packet[it];
		state.end_index += xx packet.count;
	}
	return 0;
}

usb_buffer : [..] u8;
read_next :: (serial: USB.Serial, $n: int, allocator: Allocator = .{}) -> [] Packet, bool {
	a := ifx allocator.proc then allocator else context.allocator;

	buffer : [n * size_of(Packet)] u8;
	USB.read_to(serial, buffer);

	array_add(*usb_buffer, ..buffer);

	for 0..usb_buffer.count-size_of(Packet) {
		read : *Packet = cast(*Packet) (usb_buffer.data + it);
		if memcmp(read.sync_start.data, SYNC_START_SEQ.data, SYNC_START_SEQ.count) != 0 continue;

		packets := NewArray(n - (it / size_of(Packet)), Packet, false, a);
		for * packets { <<it = read[it_index]; }

		beg := it + size_of(Packet) * n;
		for 0..usb_buffer.count-beg-1 usb_buffer[it] = usb_buffer[beg + it];
		usb_buffer.count -= beg;

		return packets, true;
	}

	empty : [n] Packet = ---;
	return empty, false;
}

main :: () {
	state: State;
	init(*state);
	play(*state);

	while true {
		packet, has_packet := next(*state);
		if has_packet && packet.id == 0 print("%\n", length(packet.mag.B));
	}
}