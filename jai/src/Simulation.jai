#scope_file

#import "Basic";
#import "Math";


Tracy :: #import "Tracy";

#import "xDebug";
Physics :: #import "Physics";
Data :: #import "Data";

using Physics;
using Data;

#scope_export

Region :: struct {
	matrix: Matrixf;

	space_origin: V2f;
	space_size:   V2f;
}

transform :: (using region: Region, xi: int, yi: int) -> V2f {
	return V2f.make(
		(xi / (matrix.cols - 1.0) - 0.5) * space_size.x + space_origin.x,
		(yi / (matrix.rows - 1.0) - 0.5) * space_size.y + space_origin.y
	);
}

estimate_to :: (magnet_strength: float64, beacon: Beacon, mag: Magnetometer, region: *Region) {
	Tracy.ZONE("Estimate");
	m := V3d.make(0, 0, 1) * magnet_strength;

	USE_AVX2 :: true;
	RINGS :: 3;
	STARS :: 8;
	probes :: #run -> [RINGS * STARS] V3f {
		radius : float32 : 0.005;
		height : float32 : 0.06;
		h_off  : float32 : 0.015;
		ps : [RINGS * STARS] V3f;

		for r: 0..RINGS-1 for s: 0..STARS-1 {
			p: V3f;
			p.x = cast(float32) (radius * cos(2 * PI * s / cast(float32) STARS));
			p.y = cast(float32) (radius * sin(2 * PI * s / cast(float32) STARS));
			p.z = cast(float32) (h_off + r * height / (RINGS - 1.0));
			ps[s + r * STARS] = p;
		}

		return ps;
	};

	ms : V3f_vec;
	rs : V3f_vec;
	Bs : V3f_vec;
	for 0..7 {
		ms.xs[it] = 0;
		ms.ys[it] = 0;
		ms.zs[it] = 1.0 / probes.count;
	}

	lawv3 := Normalv3f.make(mag.B - Cast(beacon.mean, float32), Cast(beacon.std, float32));
	for yi: 0..region.matrix.rows-1 for xi: 0..region.matrix.cols-1 {
		r2d := transform(region, xi, yi);
		r := V3f.make(r2d.x, r2d.y, 0);
		r.x = r.x - beacon.pos.x;
		r.y = r.y - beacon.pos.y;
	
		b := V3f.{};

		count := 0;
		#if USE_AVX2 for 0..probes.count/8-1 {
			local_p := 1;

			for off: 0..7 {
				it := probes[count + off];

				rs.xs[off] = r.x + it.x;
				rs.ys[off] = r.y + it.y;
				rs.zs[off] = r.z + it.z;
			}
			inline magnetic_fields_avx2(ms, rs, *Bs);

			for 0..7 {
				b.x += Bs.xs[it];
				b.y += Bs.ys[it];
				b.z += Bs.zs[it];
			}

			count += 8;
		} else for count..probes.count-1 {
			mm := m / probes.count;
			rr := r + probes[it];
			rrr := Cast(rr, float64);
			bb := Cast(magnetic_field(mm, rrr), float32);
			b += bb;
		}

		p := inline pdfv3(lawv3, b);
		<<at(region.matrix, xi, yi) *= cast(float32) p;
	}
}

estimate_to :: (
	magnet_strength: float64,
	beacon: Beacon,
	mag: Magnetometer,
	region: *Region,
	pos_samples: [] V2f
) {
	m := V3d.make(0, 0, 1) * magnet_strength;

	RINGS :: 8;
	STARS :: 3;
	probes :: #run -> [RINGS * STARS] V3d {
		radius :: 0.005;
		height :: 0.06;
		h_off  :: 0.01;
		ps : [RINGS * STARS] V3d;

		for r: 0..RINGS-1 for s: 0..STARS-1 {
			p: V3d;
			p.x = radius * cos(2 * PI * s / cast(float64) STARS);
			p.y = radius * sin(2 * PI * s / cast(float64) STARS);
			p.z = h_off + r * height / (RINGS - 1.0);
			ps[s + r * STARS] = p;
		}

		return ps;
	};

	lawv3 := Normalv3d.make(Cast(mag.B, float64) - beacon.mean, beacon.std);

	samples := NewArray(pos_samples.count, float32, false, temp);

	for pos: pos_samples {
		r := V3d.make(pos.x, pos.y, 0);
		r.x -= beacon.pos.x;
		r.y -= beacon.pos.y;

		b := V3d.{};
		for probes b += magnetic_field(m / probes.count, r + it);

		p := pdfv3(lawv3, b);

		samples[it_index] = xx p;
	}

	phiexp :: (r: float32) -> float32 {
		eps := context.state.ui.interactive.eps;
		if r >= 1 / eps return 0;

		return 1.5 - exp(1.0 -1.0 / (1+ eps*eps*r*r)) / 2;
	}
	relu :: (r: float32) -> float32 {
		if r >= context.state.ui.interactive.eps
			return 0;
		return r;
	}
	// regular := RBF(pos_samples, samples, phiexp, region, temp);
	// regular := RBF(pos_samples, samples, relu, region, temp);
	{
		scope_timer("RBF");
		regular := RBF(pos_samples, samples, r => 1.0 / (r * r + context.state.ui.interactive.eps), region, temp);
		for regular.elements region.matrix.elements[it_index] *= it;
	}
}

Averaging :: enum {
	Avg :: 0;
	Avg2 :: 1;
	Avg3 :: 1;
}

estimate_position_avg :: (region: Region) -> V2f {
	sumP: V2f = .{0, 0};
	sumW : float32 = 0;

	for yi: 0..region.matrix.rows-1 for xi: 0..region.matrix.cols-1 {
		p := transform(region, xi, yi);
		w := <<at(region.matrix, xi, yi);

		sumP += p * w;
		sumW += w;
	}

	return sumP / sumW;
}

estimate_position_avg2 :: (region: Region) -> V2f {
	sumP: V2f = .{0, 0};
	sumW : float32 = 0;

	for yi: 0..region.matrix.rows-1 for xi: 0..region.matrix.cols-1 {
		p := transform(region, xi, yi);
		w := <<at(region.matrix, xi, yi);
		w *= w;

		sumP += p * w;
		sumW += w;
	}

	return sumP / sumW;
}


estimate_position_avg3 :: (region: Region) -> V2f {
	sumP: V2f = .{0, 0};
	sumW : float32 = 0;

	for yi: 0..region.matrix.rows-1 for xi: 0..region.matrix.cols-1 {
		p := transform(region, xi, yi);
		w := <<at(region.matrix, xi, yi);
		w *= w * w;

		sumP += w * p;
		sumW += w;
	}

	return sumP / sumW;
}


