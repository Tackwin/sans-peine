#scope_file

#import "Physics";
#import "Data";

#scope_export

SYNC_START_SEQ :: u8.[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 ];
SYNC_END_SEQ   :: u8.[ 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 ];

Packet :: struct {
	sync_start: [SYNC_START_SEQ.count] u8 #align 1;
	union {
		mag: Magnetometer #align 1;
		acc: Accelerometer #align 1;
		gyr: Gyroscope #align 1;
	} #no_padding;
	id: u8 #align 1;
	kind: enum u8 {
		Magnetometer  :: 0;
		Accelerometer :: 1;
		Gyroscope     :: 2;
	} #align 1;
	time: u32 #align 1;
	sync_end: [SYNC_START_SEQ.count] u8 #align 1;
} #no_padding;

Input :: struct {
	union {
		mag: Magnetometer;
		acc: Accelerometer;
		gyr: Gyroscope;
	}
	kind: enum u8 {
		Magnetometer  :: 0;
		Accelerometer :: 1;
		Gyroscope     :: 2;
	};
	id: u8;
	time: u32;
};




Init_Result :: struct {
	success := false;
}

Driver :: struct {
	memory : *u8 = null;

	size : (mem: *u8) -> int = null;
	init : (mem: *u8) -> Init_Result = null;
	shut : (mem: *u8) = null;
	play : (mem: *u8) = null;
	stop : (mem: *u8) = null;
	next : (mem: *u8) -> Input, bool = null;
}

size :: (driver: Driver) -> int{
	if !driver.size return 0;
	return driver.size(null);
}
init :: (driver: Driver) -> Init_Result {
	if !driver.init return .{};
	return driver.init(driver.memory);
}
shut :: (driver: Driver) {
	if driver.shut driver.shut(driver.memory);
}
play :: (driver: Driver) {
	if driver.play driver.play(driver.memory);
}
stop :: (driver: Driver) {
	if driver.stop driver.stop(driver.memory);
}
next :: (driver: Driver) -> Input, bool {
	if !driver.next return .{}, false;
	a, b := driver.next(driver.memory);
	return a, b;
}
