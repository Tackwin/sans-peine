#scope_file

#import "Math";

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

magnetic_fields_avx2 :: (ms: V3f_vec, r: V3f_vec, Bs: V3f_vec) {
	C :: cast(T) (3 * u0 / (4 * PI));

	for 0..7 {
		lr := r.xs[it] * r.xs[it] + r.ys[it] * r.ys[it]+ r.zs[it] * r.zs[it];
		lr = sqrt(lr);
		rr : V3f = .{ r.xs[it] / lr, r.xs[it] / lr, r.xs[it] / lr };
		b := ms.xs[it] * rr.x + ms.ys[it] * rr.y + ms.zs[it] * rr.z;
		Bs.xs[i] = b * rr.x;
		Bs.ys[i] = b * rr.y;
		Bs.zs[i] = b * rr.z;
		Bs.xs[i] -= ms.xs[it];
		Bs.ys[i] -= ms.ys[it];
		Bs.zs[i] -= ms.zs[it];
		Bs.xs[i] /= lr * lr * lr;
		Bs.ys[i] /= lr * lr * lr;
		Bs.zs[i] /= lr * lr * lr;
		Bs.xs[i] *= C;
		Bs.ys[i] *= C;
		Bs.zs[i] *= C;
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

pdf :: (law: Normal($T), x: T) -> T {
	p : T = 1;
	p  = -0.5 * (x - law.mean) * (x - law.mean) / (law.std * law.std);
	p  = pow(cast(T) e, p);
	p /= sqrt(2 * PI);
	return p;
}

pdfv3 :: (law: Normal($T), x: T) -> T.T {
	p : T.T = 1;
	p *= pdf(Normal(T.T).make(law.mean.x, law.std.x), x.x);
	p *= pdf(Normal(T.T).make(law.mean.y, law.std.y), x.y);
	p *= pdf(Normal(T.T).make(law.mean.z, law.std.z), x.z);
	return p;
}