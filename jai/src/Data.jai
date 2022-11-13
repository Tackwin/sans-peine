
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

Tracker :: struct {
	offset: V2f;
	speed: V2f;

	acc: V3f;
	gyr: V3f;

	acc_mean_prev: V3d;
	acc_mean: V3d;
	acc_s_prev: V3d;
	acc_s: V3d;
	acc_var: V3d;
	acc_std: V3d;
	acc_samples_count := 0;
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
	accelerometer_readings: [2] Reading(Accelerometer);
	gyroscope_readings: [2] Reading(Gyroscope);

	time := 0;
	valid := false;

	estimated_pos: V2f;

}

system_is_complete :: (using this: System_Reading) -> bool {
	for magnetometer_readings if !it.valid return false;
	for accelerometer_readings if !it.valid return false;
	return true;
}

interpolate :: (a: Magnetometer, b: Magnetometer, t: float) -> Magnetometer {
	c: Magnetometer = ---;
	c.B = (b.B - a.B) * t + a.B;
	return c;
}