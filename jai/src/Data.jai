
#scope_file

#import "Physics";

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