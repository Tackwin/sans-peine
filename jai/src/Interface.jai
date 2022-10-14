#scope_file

#import "Reading";

#scope_export

Init_Result :: struct {
	success := false;
}

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
	return driver.size(null);
}
init :: (driver: Driver) -> Init_Result {
	return driver.init(driver.memory);
}
shut :: (driver: Driver) {
	driver.shut(driver.memory);
}
play :: (driver: Driver) {
	driver.play(driver.memory);
}
stop :: (driver: Driver) {
	driver.stop(driver.memory);
}
next :: (driver: Driver) -> Input, bool {
	a, b := driver.next(driver.memory);
	return a, b;
}
