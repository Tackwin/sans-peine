

#scope_file
ImGui :: #import "my_imgui";
ImGui_GL :: #import "my_imgui/imgui_impl_gl";
Render :: #import "Render";
Interface :: #import "Interface";
Physics :: #import "Physics";
Data :: #import "Data";
Simulation :: #import "Simulation";
Debug :: #import "xDebug";
Tracy :: #import "Tracy";
CLI :: #import "Command_Line";
UI :: #import "UI";

#import "Random";
#import "File";
#import "File_Watcher";

#import "Window_Creation";
#import "Input";
#import "GL";
#import "Basic";
#import "Math";
#if OS == .WINDOWS #import "Windows";

#scope_export

Sensor_Hertz :: 100;
Sensor_Dt    :: 1000.0 / Sensor_Hertz;
Ui :: struct {
	imgui: ImGui_GL.UI_State;
	window: Window_Type;

	exit := false;
	demo_window_open := true;
	timeline_window_open := true;

	camera_size : Physics.V2f = .{ 0.72, 0.405 };
	camera_pos  : Physics.V2f = .{ 0.0, 0.0 };
	camera_speed := 1.0;

	data_rate := 37.5;
	data_rate_changed := false;

	pause := false;
	interactive : struct {
		open_timeline := false; @trigger;
		clear_beacon := false; @trigger;
		default_beacon := true; @trigger;

		clear_inputs := false; @trigger;
		replay_inputs___ := false;
		replay_loop   := false;

		use_beacon1 := true;
		use_beacon2 := true;

		calibrating := true;

		magnet_strength : float64 = 0.1;

		continuous := true;
		step := false;
		
		averaging_selected: Simulation.Averaging = .Avg2;
		space_rendered: Simulation.Spaces = .Magnetic;

		eps : float32 = 0.00001;
		n_probes := 300;

		record_inputs := false;
		clear_records := false; @trigger;
	};
	using interactive;

	replay_input_cursor := 0;
	inputs : struct {
		mouse_pos : Physics.V2f;
		left_clicked := false;
	};

	Dyn_Proc :: struct {
		input_timeline: (systems: *[..] Data.System_Reading) = null;
		plot: (label: string, conf: UI.PlotConfig) -> UI.PlotStatus = null;
	};
	dyn_proc: Dyn_Proc;
};

State :: struct {
	current_driver: Interface.Driver;
	ui: Ui;

	past_inputs: [..] Interface.Input;
	pending_inputs: [..] Interface.Input;

	systems: [..] Data.System_Reading;

	min_time := S64_MAX;
	max_time := 0;

	dt: float64;
	dt32: float32;

	magnetic_space: Simulation.Region;
	velocity_space: Simulation.Region;
	position_space: Simulation.Region;

	velocity_swarm: [..] Physics.V2f;
	position_swarm: [..] Physics.V2f;

	beacons : [..] Data.Beacon;
	trackers: [..] Data.Tracker;

	pos_samples : [..] Physics.V2f;

	time_record_loop_begin: s64;

	ui_dll_watcher: File_Watcher;
}

#add_context state: State;

app_input :: () {
	using context.state;

	update_mouse_pos :: () {
		mousex, mousey, success := get_mouse_pointer_position(true);
		posx, posy, sizex, sizey, success_ := get_dimensions(context.state.ui.window, true);
		context.state.ui.inputs.mouse_pos.x = (cast(float) mousex) / sizex;
		context.state.ui.inputs.mouse_pos.y = (cast(float) mousey) / sizey;
	};
	update_mouse_pos();
	io := ImGui.GetIO();
	for events_this_frame {
		ImGui_GL.ImGui_Impl_ProcessEvent(*ui.imgui, it);
		if it.type == {
			case .QUIT; ui.exit = true;
			case .KEYBOARD;
				if !io.WantCaptureMouse {
					if it.key_code == .MOUSE_BUTTON_LEFT && it.key_pressed == 0
						ui.inputs.left_clicked = true;
				}
				if io.WantCaptureKeyboard continue;
				if it.key_pressed == 0 continue;
				if it.key_code == .ESCAPE ui.exit = true;
				if it.key_code == #char "Q" ui.camera_pos.x -= ui.camera_speed * dt32;
				if it.key_code == #char "D" ui.camera_pos.x += ui.camera_speed * dt32;
				if it.key_code == #char "Z" ui.camera_pos.y += ui.camera_speed * dt32;
				if it.key_code == #char "S" ui.camera_pos.y -= ui.camera_speed * dt32;
				if it.key_code == #char "P" ui.pause ^= true;
		}
	}
}

handle_input :: (input: Interface.Input) {
	using context.state;
	if ui.pause return;

	array_add(*pending_inputs, input);

	if ui.interactive.record_inputs {
		make_directory_if_it_does_not_exist("records/");

		file, opened := file_open(
			"records/record.dat", for_writing=true, keep_existing_content=true
		);

		if opened {
			cursor := file_length(file);
			file_set_position(*file, cursor);
			file_write(*file, *input, size_of(type_of(input)));
			file_close(*file);
		}
	}
}
pen_input :: () {
	using context.state;

	if context.state.ui.interactive.replay_inputs___ {
		file, openend := read_entire_file("records/record.dat");
		inputs_read: [] Interface.Input;
		inputs_read.data = xx file.data;
		inputs_read.count = file.count / size_of(Interface.Input);

		if ui.replay_input_cursor == 0 {
			time_record_loop_begin = to_nanoseconds(current_time_monotonic());
		}
		
		if ui.replay_input_cursor >= inputs_read.count && ui.interactive.replay_loop {
			dt := to_nanoseconds(current_time_monotonic()) - time_record_loop_begin;
			Debug.append(*Debug.timer_input_replay, dt / 1_000_000_000.0);
			clear_inputs();
		}
		else {
			if ui.replay_input_cursor < inputs_read.count
				handle_input(inputs_read[ui.replay_input_cursor]);
			ui.replay_input_cursor += 1;
		}
	}
	else {
		max_it := 100;
		for 1..max_it {
			new_input, has_new_input := Interface.next(current_driver);
			if has_new_input handle_input(new_input);
			else break;
		}
	}
}

clear_inputs :: () {
	using context.state;
	past_inputs.count = 0;
	pending_inputs.count = 0;
	systems.count = 0;
	min_time = S64_MAX;
	max_time = 0;
	ui.replay_input_cursor = 0;
}

app_update :: () {
	using context.state;
	// if ui.inputs.left_clicked {
	// 	defer ui.inputs.left_clicked = false;

	// 	pos : Physics.V2f;
	// 	pos.x = ui.camera_pos.x - ui.camera_size.x / 2;
	// 	pos.x += ui.inputs.mouse_pos.x * ui.camera_size.x;
	// 	pos.y = ui.camera_pos.y - ui.camera_size.y / 2;
	// 	pos.y += ui.inputs.mouse_pos.y * ui.camera_size.y;

	// 	b: Data.Beacon = .{};
	// 	b.pos = pos;
	// 	array_add(*beacons, b);
	// }

	if ui.interactive.clear_beacon {
		defer ui.interactive.clear_beacon = false;

		beacons.count = 0;
	}

	if ui.interactive.default_beacon {
		defer ui.interactive.default_beacon = false;
		init_beacon_trackers();
	}

	if ui.interactive.clear_inputs {
		defer ui.interactive.clear_inputs = false;
		clear_inputs();
	}

	if ui.data_rate_changed && false {
		defer ui.data_rate_changed = false;

		cmd: Interface.Command;
		cmd.type = Interface.Command.DataRate;
		cmd.data_rate.hz = ui.data_rate;
		Interface.send(current_driver, cmd);
	}

	if ui.interactive.open_timeline {
		defer ui.interactive.open_timeline = false;

		ui.timeline_window_open = true;
	}
}

pen_calibrating :: (input: Interface.Input) {
	using context.state;

	min_sample_count := S64_MAX;
	if input.kind == .Magnetometer for * beacons {
		using Physics;
		B := Cast(input.mag.B, float64);
		it.samples_count += 1;

		if it.samples_count == 1 {
			it.mean = B;
			it.s = .{0, 0, 0};
			it.var = it.s;
		} else {
			it.mean_prev = it.mean;
			it.s_prev    = it.s;
			it.mean = it.mean_prev + (B - it.mean_prev) / cast(float64) it.samples_count;
			it.s.x  = it.s.x     + (B.x - it.mean_prev.x) * cast(float64) (B.x - it.mean.x);
			it.s.y  = it.s.y     + (B.y - it.mean_prev.y) * cast(float64) (B.y - it.mean.y);
			it.s.z  = it.s.z     + (B.z - it.mean_prev.z) * cast(float64) (B.z - it.mean.z);
			it.var  = it.s / (it.samples_count - 1.0);
		}

		it.std.x = sqrt(it.var.x);
		it.std.y = sqrt(it.var.y);
		it.std.z = sqrt(it.var.z);

		min_sample_count = min(min_sample_count, it.samples_count);
	}

	if input.kind == .Accelerometer {
		it := *trackers[input.id];
		using Physics;
		g := Cast(input.acc.g, float64);
		it.acc_samples_count += 1;

		FIRST_BAD_SAMPLES :: 10;

		if it.acc_samples_count == 1 {
			it.acc_mean = g;
			it.acc_s = .{0, 0, 0};
			it.acc_var = it.acc_s;
		} else if it.acc_samples_count > FIRST_BAD_SAMPLES {
			it.acc_mean_prev = it.acc_mean;
			it.acc_s_prev    = it.acc_s;

			n := it.acc_samples_count - FIRST_BAD_SAMPLES;
			it.acc_mean = it.acc_mean_prev + (g - it.acc_mean_prev) / (n + 0.0);
			it.acc_s.x  = it.acc_s.x     + (g.x - it.acc_mean_prev.x) * (g.x - it.acc_mean.x);
			it.acc_s.y  = it.acc_s.y     + (g.y - it.acc_mean_prev.y) * (g.y - it.acc_mean.y);
			it.acc_s.z  = it.acc_s.z     + (g.z - it.acc_mean_prev.z) * (g.z - it.acc_mean.z);
			it.acc_var  = it.acc_s / (n - 1.0);
		}

		it.acc_std.x = sqrt(it.acc_var.x);
		it.acc_std.y = sqrt(it.acc_var.y);
		it.acc_std.z = sqrt(it.acc_var.z);

		min_sample_count = min(min_sample_count, max(it.acc_samples_count, FIRST_BAD_SAMPLES) - FIRST_BAD_SAMPLES);
	}

	if min_sample_count != S64_MAX && min_sample_count > 100 ui.interactive.calibrating = false;
}

pen_update :: (input: Interface.Input) {
	using context.state;
	context_state := *context.state;
	Debug.scope_timer("update");

	if ui.interactive.calibrating {
		pen_calibrating(input);
		return;
	}
	
	last_magnetometers: [2] Data.Reading(Data.Magnetometer);
	last_accelerometers: [2] Data.Reading(Data.Accelerometer);

	for < past_inputs {
		if it.kind == .Magnetometer {
			if !last_magnetometers[it.id].valid {
				last_magnetometers[it.id].payload = it.mag;
				last_magnetometers[it.id].valid = true;
				last_magnetometers[it.id].interpolated_beg = it.time;
			}
		}

		if it.kind == .Accelerometer {
			if !last_accelerometers[it.id].valid {
				last_accelerometers[it.id].payload = it.acc;
				last_accelerometers[it.id].valid = true;
				last_accelerometers[it.id].interpolated_beg = it.time;
			}
		}

		if last_magnetometers[0].valid && last_magnetometers[1].valid
		&& last_accelerometers[0].valid && last_accelerometers[1].valid
			break;
	}


	min_time = min(min_time, input.time);

	add_system :: (t: int) #expand {
		system: Data.System_Reading;
		system.time = t;
		for * system.magnetometer_readings if last_magnetometers[it_index].valid {
			it.interpolated_beg = last_magnetometers[it_index].interpolated_beg;
			it.beg = last_magnetometers[it_index].payload;
		}
		for * system.accelerometer_readings if last_accelerometers[it_index].valid {
			it.interpolated_beg = last_accelerometers[it_index].interpolated_beg;
			it.beg = last_accelerometers[it_index].payload;
		}

		array_add(*systems, system);
	}

	if max_time <= input.time {
		max_time = input.time;
		t := max_time - min_time;
		n_systems := cast(int)ceil(cast(float32) t / Sensor_Dt);

		for systems.count..n_systems add_system(cast(int) (min_time + it * Sensor_Dt));
	}

	system_to_compute := NewArray(systems.count, s64, false, temp);
	system_to_compute.count = 0;
	mark_to_compute :: (idx: int) #expand {
		`system_to_compute.count += 1;
		`system_to_compute[`system_to_compute.count - 1] = idx;
	}
	for *system, it_index: systems {
		if input.kind == .Magnetometer {
			mag := *system.magnetometer_readings[input.id];

			// If the incomming input is a measurement with the exact same time we update it
			if system.time == input.time {
				mag.payload = input.mag;
				mag.valid = true;
				mag.measured = true;
				mag.interpolated_beg = input.time;
				mag.interpolated_end = input.time;
				mark_to_compute(it_index);
			}

			// If the incomming input is within the interpolated range we update and mark it too.
			if (input.time < system.time && input.time > mag.interpolated_beg)
			|| (input.time > system.time && input.time < mag.interpolated_end)
			{
				mag.payload = input.mag;
				mag.valid = true;
				mag.measured = false;
				if input.time > system.time
				{
					mag.interpolated_end = min(mag.interpolated_end, input.time);
					mag.end = input.mag;
				}
				if input.time < system.time
				{
					mag.interpolated_beg = max(mag.interpolated_beg, input.time);
					mag.beg = input.mag;
				}
				mark_to_compute(it_index);
			}
		}

		if input.kind == .Accelerometer {
			acc := *system.accelerometer_readings[input.id];

			if system.time == input.time {
				acc.payload = input.acc;
				acc.valid = true;
				acc.measured = true;
				acc.interpolated_beg = input.time;
				acc.interpolated_end = input.time;
				mark_to_compute(it_index);
			}

			if (input.time < system.time && input.time > acc.interpolated_beg)
			|| (input.time > system.time && input.time < acc.interpolated_end)
			{
				acc.payload = input.acc;
				acc.valid = true;
				acc.measured = false;
				if input.time > system.time {
					acc.interpolated_end = min(acc.interpolated_end, input.time);
					acc.end = input.acc;
				}
				if input.time < system.time {
					acc.interpolated_beg = max(acc.interpolated_beg, input.time);
					acc.beg = input.acc;
				}
				mark_to_compute(it_index);
			}
		}
	}

	for idx: system_to_compute {
		it := *systems[idx];

		// if it isn't valid and the two magnetiometers reading didn't happen to flip to valid this
		// time around, then we don't continue, we wait for both reading to be true.
		if !it.valid
			if !Data.system_is_complete(it) continue;
			// Else we can mark it as valid
			else it.valid = true;
	}

	n := 0;
	for idx: system_to_compute if systems[idx].valid {
		it := *systems[idx];
		update_system(it);
		n += 1;
	}
}

compute_velocity_position_from_cache_at :: (
	velocity_swarm: *[..] Physics.V2f, position_swarm: *[..] Physics.V2f, at: int
) {
	using, except(velocity_swarm, position_swarm) context.state;
	if systems.count == 0 return;

	last_system := -1;
	for < at..0 if !systems[it].dirty {
		last_system = it;
		break;
	}
	if last_system >= systems.count return;
	if last_system < 0 {
		array_reserve(*systems[0].velocity_swarm, 20);
		array_reserve(*systems[0].position_swarm, 20);

		systems[0].velocity_swarm.count = 20;
		systems[0].position_swarm.count = 20;

		for * systems[0].velocity_swarm { <<it = .{0, 0}; }
		for * systems[0].position_swarm { <<it = .{0, 0}; }

		last_system = 0;
	}

	array_copy(velocity_swarm, systems[last_system].velocity_swarm);
	array_copy(position_swarm, systems[last_system].position_swarm);

	for last_system..systems.count-1 {

		Simulation.estimate_to(
			Sensor_Dt, trackers[0], systems[it].accelerometer_readings[0].payload, velocity_swarm
		);
		Simulation.estimate_to(Sensor_Dt, <<velocity_swarm, position_swarm);
	}
}

refresh_cache :: () {
	using context.state;
	
	CACHE_INTERVAL :: 10;

	time_since_last_cache := 0;
	for * systems if it.valid {
		if it.dirty {
			time_since_last_cache += 1;
		} else {
			time_since_last_cache = 0;
		}

		if time_since_last_cache > CACHE_INTERVAL {
			v: [..] Physics.V2f;
			p: [..] Physics.V2f;
			compute_velocity_position_from_cache_at(*v, *p, it_index);
			it.velocity_swarm = v;
			it.position_swarm = p;
			it.dirty = false;
			time_since_last_cache = 0;
		}
	}
}

update_system_mark := "update_system";
update_system :: (system: *Data.System_Reading) {
	if system.computed_time > 0 return;
	defer system.computed_time += 1;
	using context.state;
	Tracy.___tracy_emit_frame_mark_start(update_system_mark.data);
	defer Tracy.___tracy_emit_frame_mark_end(update_system_mark.data);
	Tracy.ZONE("update_system");
	Basis_N :: 5;

	// pos := NewArray(Basis_N*Basis_N, Physics.V2f, false, temp);
	// for x: 0..Basis_N-1 for y: 0..Basis_N-1 {
	// 	p := Physics.V2f.make(
	// 		region.matrix.rows * (x + 0.5) / Basis_N,
	// 		region.matrix.cols * (y + 0.5) / Basis_N
	// 	);
	// 	pos[x + y * Basis_N] = p;
	// }
	
	// pos_samples.count = 0;
	// for pos array_add(*pos_samples, Simulation.transform(region, xx it.x, xx it.y));
	// cdf := Physics.copy(region.matrix);
	// {
	// 	sum := 0.0;
	// 	for cdf.elements sum += it;
	// 	if sum != 0 for * cdf.elements { <<it /= sum; }
	// }
	// pos    = Physics.sample_from(cdf, ui.interactive.n_probes - Basis_N * Basis_N, temp);
	
	// for a: pos {
	// 	at := Simulation.transform(region, xx a.x, xx a.y);
	// 	found := false;
	// 	for b: pos_samples if at.x == b.x && at.y == b.y {
	// 		found = true;
	// 		break;
	// 	}
	// 	if !found array_add(*pos_samples, at);
	// }

	for * magnetic_space.matrix.elements { <<it = 1; }
	for b, b_idx: beacons {
		if !(b.std.x > 0 && b.std.y > 0 && b.std.z > 0) continue;

		reading := system.magnetometer_readings[b_idx];
		mag := system.magnetometer_readings[b_idx].payload;
		if reading.interpolated_beg != S64_MIN && reading.interpolated_end != S64_MAX {
			if !reading.measured {
				t := cast(float) (reading.interpolated_end - reading.interpolated_beg);
				t = (system.time - reading.interpolated_beg) / t;
				mag = Data.interpolate(reading.beg, reading.end, t);
			}
		}
		USE_AVX2 :: false;

		#if USE_AVX2 {
			Simulation.estimate_to_avx2(ui.interactive.magnet_strength, b, mag, *magnetic_space);
		} else {
			Simulation.estimate_to(ui.interactive.magnet_strength, b, mag, *magnetic_space);
		}
	}

	// avg_tracker: Data.Tracker;
	// avg_tracker.acc_mean.x = +trackers[1].acc_mean.z;
	// avg_tracker.acc_mean.y = +trackers[1].acc_mean.y;
	// avg_tracker.acc_mean.z = -trackers[1].acc_mean.x;
	// avg_tracker.acc_std.x = trackers[1].acc_std.z;
	// avg_tracker.acc_std.y = trackers[1].acc_std.y;
	// avg_tracker.acc_std.z = trackers[1].acc_std.x;

	// avg_tracker.acc_mean.x -= trackers[0].acc_mean.z;
	// avg_tracker.acc_mean.y -= trackers[0].acc_mean.y;
	// avg_tracker.acc_mean.z -= trackers[0].acc_mean.x;
	// avg_tracker.acc_std.x += trackers[0].acc_std.z;
	// avg_tracker.acc_std.y += trackers[0].acc_std.y;
	// avg_tracker.acc_std.z += trackers[0].acc_std.x;

	// avg_tracker.acc_mean.x /= 2;
	// avg_tracker.acc_std.x  /= 2;
	// avg_tracker.acc_mean.y /= 2;
	// avg_tracker.acc_std.y  /= 2;
	// avg_tracker.acc_mean.z /= 2;
	// avg_tracker.acc_std.z  /= 2;

	// acc: Data.Accelerometer;
	// acc.g.x = +system.accelerometer_readings[1].payload.g.z;
	// acc.g.y = +system.accelerometer_readings[1].payload.g.y;
	// acc.g.z = -system.accelerometer_readings[1].payload.g.x;
	// acc.g.x -= system.accelerometer_readings[0].payload.g.z;
	// acc.g.y -= system.accelerometer_readings[0].payload.g.y;
	// acc.g.z -= system.accelerometer_readings[0].payload.g.x;


	// if system.time == max_time {
	// 	dt := 1.0 / ui.Sensor_Hertz;

	// 	for t: trackers Simulation.estimate_to(
	// 		dt / trackers.count, t, system.accelerometer_readings[it_index].payload, *velocity_swarm
	// 	);
	// 	Simulation.estimate_to(
	// 		dt, velocity_swarm, *position_swarm
	// 	);
	// }

	ma : float32 = -FLOAT32_MAX;
	mi : float32 = +FLOAT32_MAX;
	for magnetic_space.matrix.elements ma = max(ma, it);
	for magnetic_space.matrix.elements mi = min(mi, it);
	if ma == 0 ma = 1;
	for * magnetic_space.matrix.elements { <<it /= ma; };

	if !ui.interactive.calibrating if ui.interactive.averaging_selected == {
		case .Avg;
			system.estimated_pos = Simulation.estimate_position_avg(magnetic_space);
		case .Avg2;
			system.estimated_pos = Simulation.estimate_position_avg2(magnetic_space);
		case .Avg3;
			system.estimated_pos = Simulation.estimate_position_avg3(magnetic_space);
	}
	
	// refresh_cache();
}

render_plot_mag :: () {
	using context.state;
	ys: [2][]float;
	for * ys { <<it = NewArray(past_inputs.count, float, false, temp); }

	idx1 := 0;
	idx2 := 0;
	for past_inputs if it.kind == .Magnetometer {
		if it.id == 0 {
			ys[0][idx1] = Physics.length(it.mag.B);
			idx1 += 1;
		}
		if it.id == 1 {
			ys[1][idx2] = Physics.length(it.mag.B);
			idx2 += 1;
		}
	}
	ys[0].count = idx1;
	ys[1].count = idx2;

	cys: [2]*float;
	cys[0] = ys[0].data;
	cys[1] = ys[1].data;
	plot: UI.PlotConfig;
	plot.values.count = S64_MAX;
	for ys plot.values.count = min(plot.values.count, it.count);
	plot.values.ys_list = cys.data;
	plot.values.ys_count = cys.count;
	plot.frame_size.x = ImGui.GetWindowSize().x * 0.95;
	plot.frame_size.y = 100;

	plot.scale.min = FLOAT32_MAX;
	plot.scale.max = -FLOAT32_MAX;
	for y: ys for y plot.scale.min = min(plot.scale.min, it - 0.000_000_1);
	for y: ys for y plot.scale.max = max(plot.scale.max, it + 0.000_000_1);

	ui.dyn_proc.plot("Magnet strength", plot);
}

render_plot_acc :: () {
	using context.state;
	
	ys: [2][]float;
	for * ys { <<it = NewArray(past_inputs.count, float, false, temp); }

	idx1 := 0;
	idx2 := 0;
	for past_inputs if it.kind == .Accelerometer {
		if it.id == 0 {
			ys[0][idx1] = Physics.length(it.acc.g);
			idx1 += 1;
		}
		if it.id == 1 {
			ys[1][idx2] = Physics.length(it.acc.g);
			idx2 += 1;
		}
	}
	ys[0].count = idx1;
	ys[1].count = idx2;

	cys: [2]*float;
	cys[0] = ys[0].data;
	cys[1] = ys[1].data;
	plot: UI.PlotConfig;
	plot.values.count = S64_MAX;
	for ys plot.values.count = min(plot.values.count, it.count);
	plot.values.ys_list = cys.data;
	plot.values.ys_count = cys.count;
	plot.frame_size.x = ImGui.GetWindowSize().x * 0.95;
	plot.frame_size.y = 100;

	plot.scale.min = FLOAT32_MAX;
	plot.scale.max = -FLOAT32_MAX;
	for y: ys for y plot.scale.min = min(plot.scale.min, it - 0.000_000_1);
	for y: ys for y plot.scale.max = max(plot.scale.max, it + 0.000_000_1);

	ui.dyn_proc.plot("Acceleration strength", plot);
}

render_imgui :: () {
	using context.state;

	dockspace_flags : ImGui.DockNodeFlags = .PassthruCentralNode;
	window_flags : ImGui.WindowFlags = .NoDocking;
	viewport := ImGui.GetMainViewport();
	ImGui.SetNextWindowPos(viewport.WorkPos);
	ImGui.SetNextWindowSize(viewport.WorkSize);
	ImGui.SetNextWindowViewport(viewport.ID_);
	ImGui.PushStyleVar(.WindowRounding, 0.0);
	ImGui.PushStyleVar(.WindowBorderSize, 0.0);
	window_flags |= ImGui.WindowFlags.NoTitleBar | .NoCollapse;
	window_flags |= ImGui.WindowFlags.NoResize | .NoMove;
	window_flags |= ImGui.WindowFlags.NoBringToFrontOnFocus | .NoNavFocus;
	window_flags |= ImGui.WindowFlags.NoBackground;

	ImGui.SetNextWindowBgAlpha(0);
	ImGui.PushStyleVar(.WindowPadding, .{ 0.0, 0.0 });
	ImGui.Begin("DockSpace", null, window_flags);
	defer ImGui.End();
	ImGui.PopStyleVar(3);

	dockspace_id := ImGui.GetID("MyDockSpace");
	ImGui.DockSpace(dockspace_id, .{ 0.0, 0.0 }, dockspace_flags);

	timeline_ :: () {
		using context.state;
		defer ImGui.End();
		if !ImGui.Begin("Timeline", *ui.timeline_window_open) return;

		if ui.dyn_proc.input_timeline ui.dyn_proc.input_timeline(*systems);
	} timeline_();

	ImGui.Begin("Test");
	defer ImGui.End();

	if ImGui.CollapsingHeader("Interactive", .DefaultOpen)
	{
		ImGui.inspect("State", *context.state.ui.interactive);
		if ImGui.InputFloat("Data rate", *context.state.ui.data_rate, 1)
			context.state.ui.data_rate_changed = true;
	}
	
	if ImGui.CollapsingHeader("Timer", .DefaultOpen)
	{
		text: string;
		text = tprint("Update %, % FPS\0", Debug.timer_update.value, 1.0 / Debug.timer_update.value);
		ImGui.Text(text.data);
		text = tprint("RBF %, % FPS\0", Debug.timer_RBF.value, 1.0 / Debug.timer_RBF.value);
		ImGui.Text(text.data);
		text = tprint("solve %, % FPS\0", Debug.timer_solve.value, 1.0 / Debug.timer_solve.value);
		ImGui.Text(text.data);
		text = tprint("loop %, % FPS\0", Debug.timer_loop.value, 1.0 / Debug.timer_loop.value);
		ImGui.Text(text.data);

		ImGui.Indent();
		for Debug.timer_input_replay.values {
			text = tprint("record loop %, % FPS\0", it, 1.0 / it);
			ImGui.Text(text.data);
		}
		ImGui.Unindent();
	}
	
	if ImGui.CollapsingHeader("Info", .DefaultOpen) {
		ImGui.Separator();
		for beacons {
			using Physics;
			std := it.std * 1000 * 1000;
			mean := it.mean * 1000 * 1000;
			ImGui.Text(tprint("std  (µT): % % %\0", std.x, std.y, std.z).data);
			ImGui.Text(tprint("mean (µT): % % %\0", mean.x, mean.y, mean.z).data);
		}
		ImGui.Separator();
		for trackers {
			using Physics;
			std := it.acc_std * 1000;
			mean := it.acc_mean * 1000;
			ImGui.Text(tprint("std  (mg): % % %\0", std.x, std.y, std.z).data);
			ImGui.Text(tprint("mean (mg): % % %\0", mean.x, mean.y, mean.z).data);
		}
		ImGui.Text(
			tprint("Inputs: Pending % Past %\0", pending_inputs.count, past_inputs.count).data
		);
		ImGui.Text(tprint("Systems: %\0", systems.count).data);
	}
	
	render_plot_mag();
	render_plot_acc();
}

render_magnetic_space :: () {
	using context.state;

	Render.matrix(magnetic_space.matrix, magnetic_space.space_size);
	for beacons {
		Render.circle(it.pos, 0.01);
	}
	{
		rs := NewArray(pos_samples.count, float, false, temp);
		cl := NewArray(pos_samples.count, Physics.V3f, false, temp);
		for * rs { <<it = 0.0005; }
		for * cl { <<it = .{1.0, 1.0, 1.0}; }
		Render.circles(pos_samples, rs, cl);
	}

	Vertex :: struct {
		pos: Physics.V2f #align 1;
		thick : float32 = 0.003 #align 1;
		color := Physics.V3f.{1, 1, 1} #align 1;
	}
	vertices := NewArray(systems.count, Vertex, false, temp);
	vertices.count = 0;
	for systems if it.valid {
		vertices.count += 1;
		t := it_index / (vertices.count - 1.0);
		vertices[vertices.count - 1].pos = it.estimated_pos;
		vertices[vertices.count - 1].thick = 0.0001;
		vertices[vertices.count - 1].color = Render.magma_colormap(t);
	}

	verticesu8: [] u8;
	verticesu8.data = xx vertices.data;
	verticesu8.count = size_of(Vertex) * vertices.count;
	Render.polyline(
		vertices.count,
		verticesu8,
		Render.Attribute.make(0, size_of(Vertex)),
		Render.Attribute.make(8, size_of(Vertex)),
		Render.Attribute.make(12, size_of(Vertex))
	);
}

render_velocity_space :: () {
	using context.state;

	using Physics;

	for * velocity_swarm { it.x = 0; it.y = 0; }
	for * position_swarm { it.x = 0; it.y = 0; }

	for systems {
		Simulation.estimate_to(
			Sensor_Dt, trackers[0], it.accelerometer_readings[0].payload, *velocity_swarm
		);
		Simulation.estimate_to(
			Sensor_Dt, velocity_swarm, *position_swarm
		);
	}

	for trackers {
		Render.circle(it.offset, 0.003, .{0, 1, 0});
	}

	for velocity_swarm {
		Render.circle(it, 0.0005, .{1, 1, 1});
	}
}

render_position_space :: () {
	using context.state;

	for * velocity_swarm { it.x = 0; it.y = 0; }
	for * position_swarm { it.x = 0; it.y = 0; }

	for systems {
		Simulation.estimate_to(
			Sensor_Dt, trackers[0], it.accelerometer_readings[0].payload, *velocity_swarm
		);
		Simulation.estimate_to(
			Sensor_Dt, velocity_swarm, *position_swarm
		);
	}

	using Physics;
	for trackers {
		Render.circle(it.offset, 0.003, .{0, 1, 0});
	}

	for position_swarm {
		Render.circle(it, 0.0005, .{1, 1, 1});
	}
}

render :: () {
	using context.state;
	Render.render_to_framebuffer();
	glViewport(0, 0, 1920, 1080);
	glClearColor(0, 0, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT);
	defer {
		w, h := get_window_size(ui.window);
		glViewport(0, 0, xx w, xx h);
		Render.render_to_screen();
		Render.quad(.{0, 0}, .{2, 2}, Render.framebuffer_texture);
	}

	// compute_velocity_position_from_cache_at(
	// 	*velocity_swarm, *position_swarm, systems.count - 1
	// );

	if ui.interactive.space_rendered == {
		case .Magnetic; render_magnetic_space();
		case .Velocity; render_velocity_space();
		case .Position; render_position_space();
	}
}

Args :: struct {
	run: string;
}

init_beacon_trackers :: () {
	context.state.beacons.count = 0;

	b: Data.Beacon;
	b.pos = Physics.V2f.make(-0.2, +0.07);
	array_add(*context.state.beacons, b);
	b.pos = Physics.V2f.make(-0.2, -0.07);
	array_add(*context.state.beacons, b);

	t: Data.Tracker;
	t.offset = Physics.V2f.make(-0.01, 0);
	array_add(*context.state.trackers, t);
	t.offset = Physics.V2f.make(+0.01, 0);
	array_add(*context.state.trackers, t);
}

make_spaces :: () {
	using context.state;

	magnetic_space.space_origin = .{0, 0};
	magnetic_space.space_size = .{0.3, 0.3};
	magnetic_space.matrix.rows = 100;
	magnetic_space.matrix.cols = 100;
	magnetic_space.matrix.elements = NewArray(
		magnetic_space.matrix.rows * magnetic_space.matrix.cols, float32
	);
	for * magnetic_space.matrix.elements { <<it = (random_get() % 255) / 255.0; }

	velocity_space.space_origin = .{0, 0};
	velocity_space.space_size = .{0.2, 0.2};
	velocity_space.matrix = Physics.Matrixf.make(51, 51, false);

	<<Physics.at(
		velocity_space.matrix,
		velocity_space.matrix.rows / 2,
		velocity_space.matrix.cols / 2
	) = 1;

	position_space.space_origin = .{0, 0};
	position_space.space_size = .{0.2, 0.2};
	position_space.matrix = Physics.Matrixf.make(51, 51, false);

	<<Physics.at(
		position_space.matrix,
		position_space.matrix.rows / 2,
		position_space.matrix.cols / 2
	) = 1;

	array_reserve(*velocity_swarm, 20);
	array_reserve(*position_swarm, 20);
}

free_spaces :: () {
	using context.state;

	free(magnetic_space.matrix.elements.data);
	free(velocity_space.matrix.elements.data);
	free(position_space.matrix.elements.data);

	free(velocity_swarm.data);
	free(position_swarm.data);
}

main :: () {
	using context.state;

	Tracy.set_thread_name("Main");
	success, args, is_set := CLI.parse_arguments(Args);

	make_spaces();
	defer free_spaces();


	if is_set.run {
		init_beacon_trackers();

		file, openend := read_entire_file("records/record.dat");
		inputs_read: [] Interface.Input;
		inputs_read.data = xx file.data;
		inputs_read.count = file.count / size_of(Interface.Input);

		ui.interactive.calibrating = true;
		t0 := to_nanoseconds(current_time_monotonic());
		for inputs_read {
			pen_update(it);
			t1 := to_nanoseconds(current_time_monotonic());
			if (it_index % 10) == 0
				print("Handled % inputs in %s\n", it_index, (t1 - t0) / 1000000000.0);
		}

		return;
	}

	ui.window = create_window(1200, 675, "ImGui OpenGL Test");
	setup_mousemove_hook();

	gl_create_context(ui.window, 4, 5, false, true);
	gl_load(*gl);
	gl_enable_debug_output(true);

	ImGui_GL.ImGui_Impl_CreateContext(*ui.imgui, set_default_allocators=true);
	ImGui_GL.ImGui_Impl_Init(*ui.imgui, ui.window);

	context.state.ui.dyn_proc = load_ui("dll/ui.dll");
	init(*ui_dll_watcher, (watcher: *File_Watcher(), change: *File_Change, user_data: *void) {
		context.state.ui.dyn_proc = load_ui("dll/ui.dll");
	}, null, .MODIFIED);
	defer deinit(ui_dll_watcher);
	add_directories(*ui_dll_watcher, "dll");

	current_time: float64;
	last_time: float64;

	current_driver = load_driver("USB.dll");
	current_driver.memory = alloc(Interface.size(current_driver));
	defer free(current_driver.memory);

	init_res := Interface.init(current_driver);
	if !init_res.success {
		print("Error while opening driver %\n", init_res);
		// return; 
	}
	defer Interface.shut(current_driver);

	Interface.play(current_driver);
	wglSwapIntervalEXT(0);
	while !ui.exit {
		Debug.scope_timer("loop");
		at := current_time_monotonic();
		current_time = to_float64_seconds(at);
		context.state.dt = cast(float64)(current_time - last_time);
		context.state.dt32 = cast(float32)(current_time - last_time);
		defer last_time = current_time;

		process_changes(*ui_dll_watcher);
		update_window_events();
		app_input();
		app_update();

		pen_input();
		for pending_inputs {
			pen_update(it);
			array_add(*past_inputs, it);
		}
		pending_inputs.count = 0;
		
		w, h := get_window_size(ui.window);

		ImGui_GL.ImGui_Impl_NewFrame(*ui.imgui, ui.window, w, h);
		render_imgui();
		ImGui.Render();

		glViewport(0, 0, xx w, xx h);
		glClearColor(0, 0, 0, 1);
		glClear(GL_COLOR_BUFFER_BIT);

		glEnable(GL_MULTISAMPLE);
		render();
		ImGui_GL.ImGui_Impl_RenderDrawData(*ui.imgui, ImGui.GetDrawData());

		dc := GetDC(ui.window);
		SwapBuffers(dc);
		
		reset_temporary_storage();
	}
}

get_window_size :: (window: Window_Type) -> s32, s32 {
	#if OS == .WINDOWS {
		rect: RECT;
		GetClientRect(window, *rect);
		return rect.right, rect.bottom;
	}
}


opengl_debug :: (
	source: GLenum,
	type: GLenum,
	id: GLenum,
	severity: GLenum,
	length: GLsizei,
	msg: *GLchar,
	user: *void
) #c_call {
	debug_source_str :: (source: GLenum) -> string {
		sources :: string.[
		"API",   "Window System", "Shader Compiler", "Third Party", "Application",
		"Other", "Unknown"
		];
		str_idx := min(source - GL_DEBUG_SOURCE_API, sources.count - 1);
		return sources[str_idx];
	}

	debug_type_str :: (type: GLenum) -> string {
		types :: string.[
			"Error",       "Deprecated Behavior", "Undefined Behavior", "Portability",
			"Performance", "Other",               "Unknown"
		];

		str_idx := min(type - GL_DEBUG_TYPE_ERROR, types.count - 1);
		return types[str_idx];
	}

	debug_severity_str :: (severity: GLenum) -> string {
		severities :: string.[ "High", "Medium", "Low", "Unknown" ];

		str_idx := min(severity - GL_DEBUG_SEVERITY_HIGH, severities.count - 1);
		return severities[str_idx];
	}

	message: string;
	message.data = msg;
	message.count = length;
	format :: #string DONE
=========================
OpenGL:
Object Id: %
Severity : %
Type     : %
Source   : %
Message  : %
=========================
DONE;
	new_context: Context;
	push_context new_context {
		to_ignore :: GLenum.[
			131185
		];

		for to_ignore if id == it return;

		print(
			format,
			id,
			debug_severity_str(severity),
			debug_type_str(type),
			debug_source_str(source),
			message
		);
	}
}

#if OS == .WINDOWS {
	// Jai's stdlib lacks mouse move events. Install a Windows "Hook" and catch them.
	//
	// Hopefully in the future it will be easier to provide custom WindowProc functionality
	// through Jai's windows impl.
	HHOOK :: HANDLE;
	HOOKPROC :: #type (code: s32, wparam: WPARAM, lparam: LPARAM) -> LRESULT #c_call;
	WH_MOUSE :: 7;
	user32 :: #system_library "user32";
	kernel32 :: #system_library "kernel32";

	CallNextHookEx    :: (hhk: HHOOK, code: s32, wparam: WPARAM, lparam: LPARAM) -> LRESULT #foreign user32;
	SetWindowsHookExA :: (id_hook: s32, hook_fn: HOOKPROC, hmod: HINSTANCE, thread_id: DWORD) -> HHOOK #foreign user32;

	HookProc :: (code: s32, wparam: WPARAM, lparam: LPARAM) -> LRESULT #c_call {
		return CallNextHookEx(null, code, wparam, lparam);
	}

	hook: HHOOK;

	setup_mousemove_hook :: () {
		if SetWindowsHookExA(WH_MOUSE, HookProc, null, GetCurrentThreadId()) == null {
			print("SetWindowsHookExA error: %\n", GetLastError());
			assert(false);
		}
	}

	load_driver :: ($name: string) -> Interface.Driver {
		dll := LoadLibraryA(name);
		if !dll return .{};

		driver: Interface.Driver;
		driver.size = xx GetProcAddress(dll, "size");
		driver.init = xx GetProcAddress(dll, "init");
		driver.shut = xx GetProcAddress(dll, "shut");
		driver.play = xx GetProcAddress(dll, "play");
		driver.stop = xx GetProcAddress(dll, "stop");
		driver.next = xx GetProcAddress(dll, "next");
		return driver;
	}

	ui_dll: HMODULE;

	load_ui :: ($name: string) -> Ui.Dyn_Proc {
		if ui_dll FreeLibrary(ui_dll);

		make_directory_if_it_does_not_exist("temp");
		CopyFileA(name, "temp/temp.dll", 0);
		ui_dll = LoadLibraryA("temp/temp.dll");
		if !ui_dll return .{};

		ret: Ui.Dyn_Proc;
		ret.input_timeline = xx GetProcAddress(ui_dll, "input_timeline");
		ret.plot = xx GetProcAddress(ui_dll, "plot");
		return ret;
	}
}

#scope_file

USB :: #library "../../USB";
#import "GL";
