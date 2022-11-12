#scope_file

#import "Math";
Tracy :: #import "Tracy";
#import "Basic";

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

magnetic_fields_avx2 :: (ms: V3f_vec, r: V3f_vec, Bs: *V3f_vec) {
	C :: cast(float32) (3 * u0 / (4 * PI));

	pr := r.xs.data;
	pm := ms.xs.data;
	pB := Bs.xs.data;
	pidx := cast(u64) indices.data;
	pidentity := cast(u64) identity.data;

	// print("% % %\n", formatInt(cast(int) pr, 16), formatInt(cast(int) pm, 16), formatInt(cast(int) pB, 16));

	c := C;
	pc := *c;

	#asm AVX, AVX2, FMA {
		broadcastss Cs:, [pc];

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
			movups.y id:, [pidentity + it * 4];

			movdqu.y rr:, [pidx + it * 4];
			permps.y lr:, rr, lengths;

			pxor.y b:, b, b;
			fmadd231ps.y b, mx, rx;
			fmadd231ps.y b, my, ry;
			fmadd231ps.y b, mz, rz;
			mulps.y b, b, lr;

			pxor.y rr, rr, rr;
			mulps.y rr, rx, lr;
			mulps.y rr, b, rr;
			subps.y rr, rr, mx;
			fmadd231ps.y Bx, rr, id;

			mulps.y rr, ry, lr;
			mulps.y rr, b, rr;
			subps.y rr, rr, my;
			fmadd231ps.y By, rr, id;

			mulps.y rr, rz, lr;
			mulps.y rr, b, rr;
			subps.y rr, rr, mz;
			fmadd231ps.y Bz, rr, id;

			mulps.y rr, lr, Cs;
			mulps.y rr, lr, rr;
			mulps.y rr, lr, rr;
			fmadd231ps.y f, rr, id;

			add it, 8;
		}
	}

	#asm AVX, AVX2, FMA {
		mulps.y Bx, Bx, f;
		mulps.y By, By, f;
		mulps.y Bz, Bz, f;

		movups.y [pB + 00], Bx;
		movups.y [pB + 32], By;
		movups.y [pB + 64], Bz;
		zeroupper;
	}


	// for 0..7 {

	// 	// J'essaie de remplacer cette ligne par une passe en ASM
	// 	// lr := r.xs[it] * r.xs[it] + r.ys[it] * r.ys[it] + r.zs[it] * r.zs[it];
	// 	// lr = 1.0 / sqrt(lr);

	// 	// ça peut être remplacé par ça grâce au block ASM juste avant, mais j'ai un slowdown de 2x
	// 	lr := l[it];


	// 	rr := inline V3f.make(r.xs[it] * lr, r.ys[it] * lr, r.zs[it] * lr);
	// 	b := ms.xs[it] * r.xs[it] + ms.ys[it] * r.ys[it] + ms.zs[it] * r.zs[it];
	// 	b *= lr;
	// 	Bs.xs[it] = b * rr.x;
	// 	Bs.ys[it] = b * rr.y;
	// 	Bs.zs[it] = b * rr.z;
	// 	Bs.xs[it] -= ms.xs[it];
	// 	Bs.ys[it] -= ms.ys[it];
	// 	Bs.zs[it] -= ms.zs[it];
	// 	Bs.xs[it] *= lr * lr * lr;
	// 	Bs.ys[it] *= lr * lr * lr;
	// 	Bs.zs[it] *= lr * lr * lr;
	// 	Bs.xs[it] *= C;
	// 	Bs.ys[it] *= C;
	// 	Bs.zs[it] *= C;
	// }


	// #asm AVX, AVX2, FMA {
	// 	movups.y rx:, [pr];
	// 	movups.y ry:, [pr + 32];
	// 	movups.y rz:, [pr + 64];
		
	// 	mulps.y lengths:, rx, rx;
	// 	fmadd231ps.y lengths, ry, ry;
	// 	fmadd231ps.y lengths, rz, rz;
	// 	sqrtps.y lengths, lengths;
	// 	rcpps.y lengths, lengths;

	// 	movups.y [plr], lengths;
	// 	mulps.y lengths3:, lengths, lengths;
	// 	mulps.y lengths3, lengths3, lengths;
	// 	movups.y [pl3], lengths3;
	// }

	// for 0..7 {
	// 	lr := l[it];
	// 	rr := inline V3f.make(r.xs[it] * lr, r.ys[it] * lr, r.zs[it] * lr);
	// 	b := ms.xs[it] * r.xs[it] + ms.ys[it] * r.ys[it] + ms.zs[it] * r.zs[it];
	// 	b *= lr;
	// 	Bs.xs[it] = b * rr.x;
	// 	Bs.ys[it] = b * rr.y;
	// 	Bs.zs[it] = b * rr.z;
	// 	factor[it] = l3[it] * C;
	// }

	// #asm AVX, AVX2, FMA {
	// 	movups.y rx, [pB + 0x00];
	// 	movups.y ry, [pB + 0x20];
	// 	movups.y rz, [pB + 0x40];

	// 	movups.y mx:, [pm + 0x00];
	// 	movups.y my:, [pm + 0x20];
	// 	movups.y mz:, [pm + 0x40];

	// 	movups.y factors:, [pfactor];

	// 	subps.y rx, rx, mx;
	// 	subps.y ry, ry, my;
	// 	subps.y rz, rz, mz;

	// 	mulps.y rx, rx, factors;
	// 	mulps.y ry, ry, factors;
	// 	mulps.y rz, rz, factors;

	// 	movups.y [pB + 00], rx;
	// 	movups.y [pB + 32], ry;
	// 	movups.y [pB + 64], rz;

	// 	zeroupper;
	// }

	// #asm AVX, AVX2, FMA {
	// 	mulps Bx, Bx, f;
	// 	mulps By, By, f;
	// 	mulps Bz, Bz, f;

	// 	movups.y [pB + 00], Bx;
	// 	movups.y [pB + 32], By;
	// 	movups.y [pB + 64], Bz;
	// }
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