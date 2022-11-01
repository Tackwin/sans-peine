
#scope_file

#import "Physics";
M :: #import "Math";

#scope_export

Magnetometer :: struct {
	B: V3f;
}

Accelerometer :: struct {
	g: V3f;
}

Gyroscope :: struct {
	a: V3f;
}


Beacon :: struct {
	pos: V2f;

	mean_prev: V3d;
	mean: V3d;
	s_prev: V3d;
	s: V3d;
	var: V3d;
	std: V3d;
	samples_count := 0;
}

Path :: struct {
	points : [..] V2f;
}

Reading :: struct(T: Type) {
	payload: T;
	valid := false;
	measured := false;
	interpolated_beg := M.S64_MIN;
	interpolated_end := M.S64_MAX;
	beg: T;
	end: T;
}

System_Reading :: struct {
	magnetometer_readings: [2] Reading(Magnetometer);
	time := 0;
	valid := false;

	estimated_pos: V2f;
}

system_is_complete :: (using this: System_Reading) -> bool {
	for magnetometer_readings if !it.valid return false;
	return true;
}
