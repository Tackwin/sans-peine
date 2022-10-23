#scope_file

#import "Basic";
#import "Math";

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
	m := V3d.make(0, 0, 1) * magnet_strength;

	RINGS :: 4;
	STARS :: 8;
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
	law := Normald.make(lawv3.mean.z, lawv3.std.z);
	if true for yi: 0..region.matrix.rows-1 for xi: 0..region.matrix.cols-1 {
		b := V3d.{};
		for probes {
			r := V3d.make(
				(cast(float64) xi) / region.matrix.cols,
				(cast(float64) yi) / region.matrix.rows,
				0
			);
			r.x = r.x - 0.5;
			r.y = r.y - 0.5;
			r.x = region.space_origin.x + r.x * region.space_size.x;
			r.y = region.space_origin.y + r.y * region.space_size.y;
			r += it;
			r.x = r.x - beacon.pos.x;
			r.y = r.y - beacon.pos.y;

			b += magnetic_field(m / probes.count, r);
		}
		p := pdfv3(lawv3, b);
		if !isnan(p) { <<at(region.matrix, xi, yi) *= cast(float32) p; }
	}

	law = Normald.make(cast(float64) length(mag.B) - length(beacon.mean), length(beacon.std));
	if false  for yi: 0..region.matrix.rows-1 for xi: 0..region.matrix.cols-1 {
		r := V3d.make(
			(cast(float64) xi) / region.matrix.cols,
			(cast(float64) yi) / region.matrix.rows,
			-0.04
		);
		r.x = r.x - 0.5;
		r.y = r.y - 0.5;
		r.x = region.space_origin.x + r.x * region.space_size.x;
		r.y = region.space_origin.y + r.y * region.space_size.y;
		r.x = r.x - beacon.pos.x;
		r.y = r.y - beacon.pos.y;

		b := magnetic_field_strength(length(m), length(r));
		p := pdf(law, b);
		if !isnan(p) { <<at(region.matrix, xi, yi) += cast(float32) p; }
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

	RINGS :: 4;
	STARS :: 8;
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

	samples_pos := NewArray(pos_samples.count, V2f, false, temp);
	samples := NewArray(pos_samples.count, float32, false, temp);

	for pos: pos_samples {
		r := V3d.make(pos.x, pos.y, 0);
		r.x /= region.matrix.cols;
		r.y /= region.matrix.rows;
		r.x -= 0.5;
		r.y -= 0.5;
		r.x = region.space_origin.x + r.x * region.space_size.x;
		r.y = region.space_origin.y + r.y * region.space_size.y;
		r.z = -0.04;

		b := V3d.{};
		for probe: probes {
			rr := r + probe;
			rr.x -= beacon.pos.x;
			rr.y -= beacon.pos.y;
			b += magnetic_field(m / probes.count, rr);
		}

		p := pdfv3(lawv3, b);

		samples_pos[it_index] = V2f.make(xx r.x, xx r.y);
		samples[it_index] = xx p;
	}

	regular := inv_interpolate(samples, samples_pos, region, temp);
	for regular.elements region.matrix.elements[it_index] *= it;
}

Averaging :: enum {
	Avg :: 0;
}

estimate_position :: (region: Region) -> V2f {
	sumP: V2f = .{0, 0};
	sumW : float32 = 0;

	for yi: 0..region.matrix.rows-1 for xi: 0..region.matrix.cols-1 {
		p := transform(region, xi, yi);
		w := <<at(region.matrix, xi, yi);

		sumP += w * p;
		sumW += w;
	}

	return sumP / sumW;
}


