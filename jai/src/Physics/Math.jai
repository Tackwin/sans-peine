#scope_file

#import "Math";
#import "Random";
#import "Basic";
#import "Sort";
#import "xDebug";
// MKL :: #import "mkl";
Simulation :: #import "Simulation";

#scope_export

V2 :: struct(T: Type) {
	x: T;
	y: T;

	make :: (x: T, y: T) -> V2(T) #c_call {
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
	make :: (x: T, y: T, z: T) -> V3(T) #c_call {
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

V3f_vec :: struct {
	xs: [8] float32;
	ys: [8] float32;
	zs: [8] float32;

	fill :: (x: float32, y: float32, z: float32) -> V3f_vec {
		vec: V3f_vec = ---;
		for 0..7 {
			vec.xs[it] = x;
			vec.ys[it] = y;
			vec.zs[it] = z;
		}
		return vec;
	}
} #no_padding;

Cast :: (x: V3($T), $U: Type) -> V3(U) {
	return V3(U).make(cast(U) x.x, cast(U) x.y, cast(U) x.z);
}

Matrix :: struct(T: Type) {
	elements: [] T;
	rows: int;
	cols: int;

	make :: (r: int, c: int, $to_init: bool, a: Allocator = .{}) -> Matrix(T) {
		m: Matrix(T) = ---;
		m.rows = r;
		m.cols = r;
		m.elements = NewArray(r * c, T, to_init, a);
		return m;
	}
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

solve_square :: (A: Matrix($T), b: [] T, x: [] T) -> [] T {
	return x;
	// assert(A.rows == A.cols);

	// status: s64 = 0;

	// lu := NewArray(A.rows * A.cols, T, false, temp);
	// ipiv := NewArray(max(1, min(A.rows, A.cols)), s64, false, temp);

	// memcpy(lu.data, A.elements.data, A.elements.count * size_of(T));
	// status = MKL.LAPACKE_sgetrf(
	// 	MKL.LAPACK_ROW_MAJOR,
	// 	A.rows,
	// 	A.cols,
	// 	lu.data,
	// 	A.cols,
	// 	ipiv.data
	// );

	// memcpy(x.data, b.data, b.count * size_of(T));
	// status = MKL.LAPACKE_sgetrs(
	// 	MKL.LAPACK_ROW_MAJOR,
	// 	#char "N",
	// 	b.count,
	// 	1,
	// 	lu.data,
	// 	A.rows,
	// 	ipiv.data,
	// 	x.data,
	// 	1
	// );

	// return x;
}

RBF :: (
	samples_x: [] V2($T),
	samples_y: [] T,
	phi: (T) -> T,
	region: *Simulation.Region,
	allocator: Allocator = .{}
) -> Matrix(T) {
	mat := Matrix(T).make(region.matrix.rows, region.matrix.cols, false, allocator);

	A := Matrix(T).make(samples_x.count, samples_x.count, false, temp);
	f := NewArray(samples_y.count, T, false, temp);

	for a, i: samples_x for b, j: samples_x A.elements[i * A.cols + j] = phi(length(b - a));
	for samples_y f[it_index] = it;

	w := NewArray(samples_y.count, T, false, temp);
	{
		scope_timer("solve");
		solve_square(A, f, w);
	}

	for y: 0..mat.rows-1 for x: 0..mat.cols-1 {
		p := Simulation.transform(region, x, y);
		
		s : float32 = 0;
		for w s += it * phi(length(p - samples_x[it_index]));
		mat.elements[x + y * mat.cols] = s;
	}

	mi := +FLOAT32_MAX;
	ma := -FLOAT32_MAX;
	for mat.elements ma = max(ma, it);
	for mat.elements mi = min(mi, it);
	for * mat.elements { <<it = (<<it - mi) / (ma - mi); }
	return mat;
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
operator * :: (a: *V3($T), b: T) -> V3(T) #symmetric {
	return V3(T).make(a.x * b, a.y * b, a.z * b);
}
operator / :: (a: V3($T), b: T) -> V3(T) {
	return V3(T).make(a.x / b, a.y / b, a.z / b);
}
operator / :: (a: *V3($T), b: T) -> V3(T) {
	return V3(T).make(a.x / b, a.y / b, a.z / b);
}
operator - :: (a: V3($T), b: V3(T)) -> V3(T) {
	return V3(T).make(a.x - b.x, a.y - b.y, a.z - b.z);
}
operator - :: (a: *V3($T), b: V3(T)) -> V3(T) {
	return V3(T).make(a.x - b.x, a.y - b.y, a.z - b.z);
}
operator + :: (a: *V3($T), b: V3(T)) -> V3(T) {
	return V3(T).make(a.x + b.x, a.y + b.y, a.z + b.z);
}
operator + :: (a: V3($T), b: V3(T)) -> V3(T) {
	return V3(T).make(a.x + b.x, a.y + b.y, a.z + b.z);
}


operator - :: (a: V2($T), b: V2(T)) -> V2(T) {
	return V2(T).make(a.x - b.x, a.y - b.y);
}
operator + :: (a: *V2($T), b: V2(T)) -> V2(T) {
	return V2(T).make(a.x + b.x, a.y + b.y);
}
operator + :: (a: V2($T), b: V2(T)) -> V2(T) {
	return V2(T).make(a.x + b.x, a.y + b.y);
}
operator * :: (a: V2($T), b: T) -> V2(T) #symmetric {
	return V2(T).make(a.x * b, a.y * b);
}
operator * :: (a: *V2($T), b: T) -> V2(T) #symmetric {
	return V2(T).make(a.x * b, a.y * b);
}
operator / :: (a: V2($T), b: T) -> V2(T) {
	return V2(T).make(a.x / b, a.y / b);
}
operator / :: (a: *V2($T), b: T) -> V2(T) {
	return V2(T).make(a.x / b, a.y / b);
}
operator == :: (a: V2($T), b: T) -> bool {
	return a.x == b.x && a.y == b.y;
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

	samples := NewArray(n, V2(T), true, allocator);

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

		added := false;
		for cursor..mc_samples.count-1 {
			if mc_samples[it].rng > sum break;

			samples[mc_samples[it].idx].x = cast(T) x;
			samples[mc_samples[it].idx].y = cast(T) y;
			cursor += 1;
		}

		if cursor >= mc_samples.count break;
	}

	samples.count = cursor;
	return samples;
}

