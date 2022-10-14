#scope_file

#import "Reading";
#import "Math";

#scope_export

Matrix :: struct(T: Type) {
	elements: [] T;
	rows: int;
	cols: int;
}

Matrixd :: Matrix(float64);
Matrixf :: Matrix(float32);

length :: (v: V2($T)) -> T {
	return sqrt(v.x * v.x + v.y * v.y);
}

length :: (v: V3($T)) -> T {
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

