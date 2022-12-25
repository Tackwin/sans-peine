#scope_file
#import "Basic";

USB :: #import "USB";
Interface :: #import "Interface";
#import "Thread";
#import "Physics";

#if OS == .WINDOWS {
	Windows :: #import "Windows";
	#import "Windows_Utf8";
	SetThreadDescription :: (thread: Windows.HANDLE, name: Windows.LPCWSTR) #foreign kernel32;

	kernel32    :: #system_library "kernel32";
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
#program_export "send"
_send :: (memory: *u8, c: Interface.Command) -> bool {
	state: *State = cast(*State) memory;
	return send(state, c);
}

#scope_export

State :: struct {
	serial: USB.Serial;
	thread: *Thread = null;
	mutex: Mutex;

	thread_launched := false;
	running := false;

	INPUT_RING_BUFFER_SIZE : u64 : 256;
	inputs : [INPUT_RING_BUFFER_SIZE] Interface.Packet;
	beg_index : u64 = 0;
	end_index : u64 = 0;

	usb_buffer : [..] u8;
}

size :: (state: *State) -> int {
	return size_of(State);
}

init :: (state: *State) -> Interface.Init_Result {
	<<state = .{};
	result: Interface.Init_Result;
	state.serial = USB.Serial.make("\\\\.\\COM5");
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
		array_reserve(*state.usb_buffer, 1024*1024*32);
	}
	state.running = true;
}

stop :: (state: *State) {
	state.running = false;
}

next :: (state: *State) -> Interface.Input, bool {
	if state.beg_index >= state.end_index return .{}, false;

	packet : Interface.Packet = ---;
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
			out.mag.B /= 1_000_000;
		case .Accelerometer; out.acc = packet.acc;
		case .Gyroscope;     out.gyr = packet.gyr;
	}
	return out, true;
}

send :: (state: *State, c: Interface.Command) -> bool {
	Packet :: struct {
		size: u8;
		kind: u8;
		payload: union {
			data_rate: Interface.Command.DataRate;
		}
		using payload;
	}

	packet: Packet = ---;
	packet.kind = 0;
	memcpy(*packet.payload, *c.data_rate, size_of(type_of(packet.data_rate)));

	buffer: [] u8;
	buffer.data = cast(*u8) *packet;
	buffer.count = size_of(Packet);

	buffer[0] = size_of(Packet);
	{
		lock(*state.mutex);
		defer unlock(*state.mutex);
		USB.write(state.serial, buffer);
	}
	return true;
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
		packet, has_packet := read_next(state, 8);
		if !has_packet continue;

		lock(*state.mutex);
		defer unlock(*state.mutex);
		
		for 0..packet.count-1
			state.inputs[(state.end_index + xx it) % State.INPUT_RING_BUFFER_SIZE] = packet[it];
		state.end_index += xx packet.count;
	}
	return 0;
}

read_next :: (
	state: *State, $n: int, allocator: Allocator = .{}
) -> [] Interface.Packet, bool {
	a := ifx allocator.proc then allocator else context.allocator;

	serial := *state.serial;

	buffer : [n * size_of(Interface.Packet)] u8;
	USB.read_to(serial, buffer);

	// if (state.usb_buffer.count + buffer.count < state.usb_buffer.allocated)
	array_add(*state.usb_buffer, ..buffer);

	packets := NewArray(n, Interface.Packet, false, a);
	packets_added := 0;

	beg := 0;
	for 0..state.usb_buffer.count-size_of(Interface.Packet)-1 {
		beg += 1;
		read : *Interface.Packet = cast(*Interface.Packet) (state.usb_buffer.data + it);
		cmp := memcmp(
			read.sync_start.data, Interface.SYNC_START_SEQ.data, Interface.SYNC_START_SEQ.count
		);
		if cmp != 0 continue;

		packets[packets_added] = <<read;
		packets_added += 1;
		if packets_added >= packets.count break;

		// state.usb_buffer.count = max(state.usb_buffer.count, 0);

	}
	for 0..state.usb_buffer.count-beg-1 state.usb_buffer[it] = state.usb_buffer[beg + it];
	state.usb_buffer.count -= beg;

	return packets, true;
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