#scope_file

#import "Basic";
#import "Math";
#import "Random";

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

estimate_to :: (dt: float64, tracker: Tracker, acc: Accelerometer, swarm: *[] V2f) {
	Tracy.ZONE("Estimate acc");

	lawv3 := Normalv3f.make(acc.g - Cast(tracker.acc_mean, float32), Cast(tracker.acc_std, float32));

	for * <<swarm {
		sx := sample(Normalf.make(lawv3.mean.x, lawv3.std.x));
		sy := sample(Normalf.make(lawv3.mean.y, lawv3.std.y));

		sx = lawv3.mean.x;
		sy = lawv3.mean.y;

		it.x += cast(float32) (dt * sx);
		it.y += cast(float32) (dt * sy);
	}
}

estimate_to :: (dt: float64, velocity_swarm: [] V2f, position_swarm: *[] V2f) {
	Tracy.ZONE("Estimate pos");

	vn := velocity_swarm.count;
	for * <<position_swarm {
		i := cast(int) (random_get_zero_to_one_open() * vn);
		it.x += cast(float32) (dt * velocity_swarm[i].x);
		it.y += cast(float32) (dt * velocity_swarm[i].y);
	}
}

estimate_to_avx2 :: (magnet_strength: float64, beacon: Beacon, mag: Magnetometer, region: *Region) {
	Tracy.ZONE("Estimate AVX2 mag");

	RINGS :: 4;
	STARS :: 8;
	probes :: #run -> [RINGS * STARS / 8] V3f_vec {
		radius : float32 : 0.005;
		height : float32 : 0.03;
		h_off  : float32 : 0.035;
		probe_aos : [RINGS * STARS] V3f;

		for r: 0..RINGS-1 for s: 0..STARS-1 {
			p: V3f;
			p.x = cast(float32) (radius * cos(2 * PI * s / cast(float32) STARS));
			p.y = cast(float32) (radius * sin(2 * PI * s / cast(float32) STARS));
			p.z = cast(float32) (h_off + r * height / (RINGS - 1.0));
			probe_aos[s + r * STARS] = p;
		}

		probe_soa : [RINGS * STARS / 8] V3f_vec;
		for probe_aos {
			probe_soa[it_index / 8].xs[it_index % 8] = it.x;
			probe_soa[it_index / 8].ys[it_index % 8] = it.y;
			probe_soa[it_index / 8].zs[it_index % 8] = it.z;
		}

		return probe_soa;
	};

	ms : V3f_vec;
	rs : V3f_vec;
	Bs : V3f_vec;
	for 0..7 {
		ms.xs[it] = 0;
		ms.ys[it] = 0;
		ms.zs[it] = cast(float32) (magnet_strength / (probes.count * 8));
	}

	pm := ms.xs.data;
	pB := Bs.xs.data;

	law := Normalf.make(length(mag.B - Cast(beacon.mean, float32)), xx length(beacon.std));
	lawv3 := Normalv3f.make(mag.B - Cast(beacon.mean, float32), Cast(beacon.std, float32));

	beg_y := -0.5 * region.space_size.y + region.space_origin.y - beacon.pos.y;
	end_y := +0.5 * region.space_size.y + region.space_origin.y - beacon.pos.y;
	ste_y := (end_y - beg_y) / (region.matrix.rows - 1.0);

	beg_x := -0.5 * region.space_size.x + region.space_origin.x - beacon.pos.x;
	end_x := +0.5 * region.space_size.x + region.space_origin.x - beacon.pos.x;
	ste_x := (end_x - beg_x) / (region.matrix.rows - 1.0);

	origin := V3f_vec.fill(beg_x, beg_y, 0);
	for yi: 0..region.matrix.rows-1 {
		for 0..7 origin.xs[it] = beg_x;

		for xi: 0..region.matrix.cols-1 {
			b := V3f.make(0, 0, 0);

			for probe_idx: 0..probes.count-1 {
				porgins := origin.xs.data;
				pprobes := (cast(s64) probes.data) + probe_idx * 3 * 8 * size_of(float32);

				r: V3f_vec = ---;
				pr := r.xs.data;
				for 0..7 {
					r.xs[it] = origin.xs[it] + probes[probe_idx].xs[it];
					r.ys[it] = origin.ys[it] + probes[probe_idx].ys[it];
					r.zs[it] = origin.zs[it] + probes[probe_idx].zs[it];
				}

				inline magnetic_fields_avx2(pr, pm, pB);

				for 0..7 {
					b.x += Bs.xs[it];
					b.y += Bs.ys[it];
					b.z += Bs.zs[it];
				}
			}

			// p := pdfv3(lawv3, b);
			p := pdf(law, length(b));
			<<at(region.matrix, xi, yi) *= cast(float32) p;

			for 0..7 origin.xs[it] += ste_x;
		}

		for 0..7 origin.ys[it] += ste_y;
	}
}

estimate_to :: (magnet_strength: float64, beacon: Beacon, mag: Magnetometer, region: *Region) {
	Tracy.ZONE("Estimate");
	m := V3d.make(0, 0, 1) * magnet_strength;

	RINGS :: 3;
	STARS :: 8;
	probes :: #run -> [RINGS * STARS] V3f {
		radius : float32 : 0.005;
		height : float32 : 0.03;
		h_off  : float32 : 0.035;
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

	law := Normalf.make(length(mag.B - Cast(beacon.mean, float32)), xx length(beacon.std));
	lawv3 := Normalv3f.make(mag.B - Cast(beacon.mean, float32), Cast(beacon.std, float32));
	for yi: 0..region.matrix.rows-1 for xi: 0..region.matrix.cols-1 {
		r2d := transform(region, xi, yi);
		r := V3f.make(r2d.x, r2d.y, 0);
		r.x = r.x - beacon.pos.x;
		r.y = r.y - beacon.pos.y;
		r.z = -0.035;
	
		// b := V3f.{};
		b := cylindrycal_field(m, Cast(r, float64), 0.025, 0.01);

		// count := 0;
		// for count..probes.count-1 {
		// 	mm := m / probes.count;
		// 	rr := r + probes[it];
		// 	rrr := Cast(rr, float64);
		// 	bb := Cast(magnetic_field(mm, rrr), float32);
		// 	b += bb;
		// }

		bl := length(b);
		// p := inline pdfv3(lawv3, Cast(b, float32));
		// print("Mesured % simulated %\n", mag.B - Cast(beacon.mean, float32), Cast(b, float32));
		p := inline pdf(law, xx bl);
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

Spaces :: enum {
	Magnetic :: 0;
	Velocity :: 1;
	Position :: 2;
	Fused :: 3;
}

Averaging :: enum {
	Avg :: 0;
	Avg2 :: 1;
	Avg3 :: 2;
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


