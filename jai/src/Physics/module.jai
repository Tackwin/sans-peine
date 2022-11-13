#scope_file

#import "Math";
Tracy :: #import "Tracy";
#import "Basic";
#import "Random";

#scope_export

Math :: #import "Physics/Math";
using Math;

u0 :: 1.25663706e-6;
PI :: 3.14159265359;
e  :: 2.71828182846;

magnetic_field :: (m: V3($T), r: V3(T)) -> V3(T) {
	C :: cast(T) (3 * u0 / (4 * PI));

	lr := length(r);
	rr := r / lr;
	b := dot(m, rr);
	B := b * rr;
	B -= m;
	B /= lr * lr * lr;
	B *= C;
	return B;
}

indices :: u32.[
	0, 0, 0, 0, 0, 0, 0, 0,
	1, 1, 1, 1, 1, 1, 1, 1,
	2, 2, 2, 2, 2, 2, 2, 2,
	3, 3, 3, 3, 3, 3, 3, 3,
	4, 4, 4, 4, 4, 4, 4, 4,
	5, 5, 5, 5, 5, 5, 5, 5,
	6, 6, 6, 6, 6, 6, 6, 6,
	7, 7, 7, 7, 7, 7, 7, 7
] #align 32;

identity :: float32.[
	1, 0, 0, 0, 0, 0, 0, 0,
	0, 1, 0, 0, 0, 0, 0, 0,
	0, 0, 1, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 0, 0, 0, 0,
	0, 0, 0, 0, 1, 0, 0, 0,
	0, 0, 0, 0, 0, 1, 0, 0,
	0, 0, 0, 0, 0, 0, 1, 0,
	0, 0, 0, 0, 0, 0, 0, 1
] #align 32;

magnetic_fields_avx2 :: (pr: *float32, pm: *float, pB: *float) {
	C :: cast(float32) (3 * u0 / (4 * PI));

	pidx := cast(u64) indices.data;
	pidentity := cast(u64) identity.data;

	c := C;
	pc := *c;

	#asm AVX, AVX2, FMA {
		movups.y rx:, [pr + 0x00];
		movups.y ry:, [pr + 0x20];
		movups.y rz:, [pr + 0x40];

		movups.y mx:, [pm + 0x00];
		movups.y my:, [pm + 0x20];
		movups.y mz:, [pm + 0x40];
		pxor.y Bx:, Bx, Bx;
		pxor.y By:, By, By;
		pxor.y Bz:, Bz, Bz;

		pxor.y f:, f, f;

		mulps.y lengths:, rx, rx;
		fmadd231ps.y lengths, ry, ry;
		fmadd231ps.y lengths, rz, rz;
		sqrtps.y lengths, lengths;
		rcpps.y lengths, lengths;
	}

	it := 0;
	while it < 64 {
		#asm AVX, AVX2, FMA {
			movdqu.y rr:, [pidx + it * 4];
			permps.y lr:, rr, lengths;

			mulps.y b:, mx, rx;
			fmadd231ps.y b, my, ry;
			fmadd231ps.y b, mz, rz;
			mulps.y b, b, lr;

			pxor.y rr, rr, rr;
			mulps.y rr, rx, lr;
			mulps.y rr, b, rr;
			subps.y rr, rr, mx;
			fmadd231ps.y Bx, rr, [pidentity + it * 4];

			mulps.y rr, ry, lr;
			mulps.y rr, b, rr;
			subps.y rr, rr, my;
			fmadd231ps.y By, rr, [pidentity + it * 4];

			mulps.y rr, rz, lr;
			mulps.y rr, b, rr;
			subps.y rr, rr, mz;
			fmadd231ps.y Bz, rr, [pidentity + it * 4];

			mulps.y rr, lr, lr;
			mulps.y rr, lr, rr;
			fmadd231ps.y f, rr, [pidentity + it * 4];

			add it, 8;
		}
	}

	#asm AVX, AVX2, FMA {
		broadcastss lengths, [pc];
		mulps.y f, f, lengths;
		mulps.y Bx, Bx, f;
		mulps.y By, By, f;
		mulps.y Bz, Bz, f;
		movups.y [pB + 0x00], Bx;
		movups.y [pB + 0x20], By;
		movups.y [pB + 0x40], Bz;
		zeroupper;
	}
}

magnetic_field_strength :: (m: $T, r: T) -> T {
	C :: cast(T) (3 * u0 / (4 * PI));

	return C * m / (r * r * r);
}

Normal :: struct(T: Type) {
	mean: T;
	std:  T;

	make :: (mean: T, std: T) -> Normal(T) {
		n: Normal(T) = ---;
		n.mean = mean;
		n.std = std;
		return n;
	}
}

Normalf :: Normal(float32);
Normald :: Normal(float64);
Normalv3f :: Normal(V3f);
Normalv3d :: Normal(V3d);

fast_exp :: (x: $T) -> T {
	S :: union {
		f: float64;
		i: s64;
	}
	s : S = ---;
	s.i = xx ((cast(float64) x) * 6497320848556798 + 4606794787981043020);
	return xx s.f;
}

pdf :: (law: Normal($T), x: T) -> T {
	C :: #run sqrt(2 * PI);
	p : T = 1;
	p  = -0.5 * (x - law.mean) * (x - law.mean) / (law.std * law.std);
	p  = exp(p);
	p /= cast(T) C;
	return p;
}

pdfv3 :: (law: Normal($T), x: T) -> T.T {
	C :: cast(T.T) #run 1.0 / (sqrt(2 * PI) * 2 * PI);

	a := (x.x - law.mean.x) / law.std.x;
	y := (x.y - law.mean.y) / law.std.y;
	z := (x.z - law.mean.z) / law.std.z;
	a *= a;
	y *= y;
	z *= z;

	e := -0.5 * (a + y + z);

	// We don't really need to multiply by C here since it's a constant and in the end we will
	// normalize the matrix, but we do it anyway.
	p := ifx e < -10 then 0 else C * exp(e);
	return p;
}

sample :: (law: Normal($T)) -> T {
	u := random_get_zero_to_one();
	v := random_get_zero_to_one();

	n := sqrt(-2 * log(u) / log(e)) * cos(2 * PI * v);

	return cast(T) (law.mean + n * law.std);
}