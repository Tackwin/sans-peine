#scope_file

#import "Math";
#import "Random";
#import "Basic";
#import "Sort";
Simulation :: #import "Simulation";

#scope_export

V2 :: struct(T: Type) {
	x: T;
	y: T;

	make :: (x: T, y: T) -> V2(T) {
		v: V2(T) = ---;
		v.x = x;
		v.y = y;
		return v;
	}
} #no_padding;

V3 :: struct(T: Type) {
	x: T;
	y: T;
	z: T;
	#place x;
	xy: V2(T);
	#place y;
	yz: V2(T);
	make :: (x: T, y: T, z: T) -> V3(T) {
		v: V3(T) = ---;
		v.x = x;
		v.y = y;
		v.z = z;
		return v;
	}
} #no_padding;

V2f :: V2(float32);
V2d :: V2(float64);
V3f :: V3(float32);
V3d :: V3(float64);

Cast :: (x: V3($T), $U: Type) -> V3(U) {
	return V3(U).make(cast(U) x.x, cast(U) x.y, cast(U) x.z);
}

Matrix :: struct(T: Type) {
	elements: [] T;
	rows: int;
	cols: int;
}

Matrixd :: Matrix(float64);
Matrixf :: Matrix(float32);

copy :: (mat: Matrix($T), allocator : Allocator = .{}) -> Matrix(T) {
	new : Matrix(T);
	new.elements = NewArray(mat.elements.count, T, false, allocator);
	new.rows = mat.rows;
	new.cols = mat.cols;
	memcpy(new.elements.data, mat.elements.data, mat.elements.count * size_of(T));
	return new;
}

print :: (m: Matrix($T), format: FormatFloat = .{}) {
	f := format;
	for yi: 0..m.rows-1 for xi: 0..m.cols-1 {
		f.value = m.elements[xi + yi * m.cols];
		print("% ", f);
	}
}

at :: (m: Matrix($T), x: int, y: int) -> *T {
	return m.elements.data + x + y * m.cols;
}

inv_interpolate :: (
	samples: [] $T,
	sample_pos: [] V2(T),
	region: *Simulation.Region,
	allocator: Allocator = .{}
) -> Matrix(T) {
	mat: Matrix(T);
	mat.rows = region.matrix.rows;
	mat.cols = region.matrix.cols;
	mat.elements = NewArray(mat.rows * mat.cols, T, false, allocator);
	for y: 0..mat.rows-1 for x: 0..mat.cols-1 {
		p := V2(T).make(
			region.space_origin.x + ((cast(T) x) / mat.cols - 0.5) * region.space_size.x,
			region.space_origin.y + ((cast(T) y) / mat.rows - 0.5) * region.space_size.y
		);

		min_dist := FLOAT64_MAX;
		min_it : T = 0;

		for samples {
			d := length(p - sample_pos[it_index]);
			if d < min_dist {
				min_dist = d;
				min_it = it;
				break;
			}
		}

		mat.elements[x + y * mat.cols] = min_it;
	}
	return mat;
}

length :: (v: V2($T)) -> T {
	return sqrt(v.x * v.x + v.y * v.y);
}

length :: (v: V3($T)) -> T {
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

dot :: (a: V3($T), b: V3(T)) -> T {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

operator * :: (a: V3($T), b: T) -> V3(T) #symmetric {
	return V3(T).make(a.x * b, a.y * b, a.z * b);
}
operator / :: (a: V3($T), b: T) -> V3(T) {
	return V3(T).make(a.x / b, a.y / b, a.z / b);
}
operator - :: (a: V3($T), b: V3(T)) -> V3(T) {
	return V3(T).make(a.x - b.x, a.y - b.y, a.z - b.z);
}
operator + :: (a: V3($T), b: V3(T)) -> V3(T) {
	return V3(T).make(a.x + b.x, a.y + b.y, a.z + b.z);
}


operator - :: (a: V2($T), b: V2(T)) -> V2(T) {
	return V2(T).make(a.x - b.x, a.y - b.y);
}
operator + :: (a: V2($T), b: V2(T)) -> V2(T) {
	return V2(T).make(a.x + b.x, a.y + b.y);
}
operator * :: (a: V2($T), b: T) -> V2(T) #symmetric {
	return V2(T).make(a.x * b, a.y * b);
}
operator / :: (a: V2($T), b: T) -> V2(T) {
	return V2(T).make(a.x / b, a.y / b);
}

random :: ($T: Type) -> T {
	#if T == float32 {
		return cast(float32) random_get_zero_to_one();
	} else #if T == float64 {
		return cast(float64) random_get_zero_to_one();
	} else #if T == V3f {
		return V3f.make(random(float32), random(float32), random(float32));
	} else #if T == V3d {
		return V3d.make(random(float64), random(float64), random(float64));
	} else return T.{};
}

sample_from :: (mat: Matrix($T), n: int, allocator : Allocator = .{}) -> [] V2(T) {
	sum_mat : T = 0;
	for mat.elements sum_mat += it;

	samples := NewArray(n, V2(T), false, allocator);

	indices := NewArray(n, int, false, temp);
	for * indices { <<it = it_index; }

	Pair :: struct {
		rng: T;
		idx: int;
	}
	mc_samples := NewArray(n, Pair, false, temp);

	for * mc_samples { it.rng = random(T); it.idx = it_index; }

	quick_sort(mc_samples, x => x.rng);

	cursor := 0;
	sum : T = 0;
	for 0..mat.rows*mat.cols-1 {
		x := it % mat.cols;
		y := it / mat.cols;
		
		sum += mat.elements[it] / sum_mat;

		for cursor..mc_samples.count-1 {
			if mc_samples[it].rng > sum break;

			samples[mc_samples[it].idx].x = cast(T) x;
			samples[mc_samples[it].idx].y = cast(T) y;
			cursor += 1;
		}

		if cursor >= mc_samples.count break;
	}

	return samples;
}

