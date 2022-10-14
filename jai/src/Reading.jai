
#scope_file

#scope_export

V2 :: struct(T: Type) {
	x: T;
	y: T;
} #no_padding;

V3 :: struct(T: Type) {
	x: T;
	y: T;
	z: T;
	#place x;
	xy: V2(T);
	#place y;
	yz: V2(T);
} #no_padding;

V2f :: V2(float32);
V2d :: V2(float64);
V3f :: V3(float32);
V3d :: V3(float64);

Reading :: struct {}

Magnetometer :: struct {
	using reading: Reading;
	B: V3f;
}

Accelerometer :: struct {
	using reading: Reading;
	g: V3f;
}

Gyroscope :: struct {
	using reading: Reading;
	a: V3f;
}

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


