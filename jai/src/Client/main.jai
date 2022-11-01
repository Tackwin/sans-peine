
#scope_file
ImGui :: #import "my_imgui";
ImGui_GL :: #import "my_imgui/imgui_impl_gl";
Render :: #import "Render";
Interface :: #import "Interface";
Physics :: #import "Physics";
Data :: #import "Data";
Simulation :: #import "Simulation";
Debug :: #import "xDebug";
#import "Random";
#import "File";

#import "Window_Creation";
#import "Input";
#import "GL";
#import "Basic";
#import "Math";
#if OS == .WINDOWS #import "Windows";

#scope_export

Ui :: struct {
	imgui: ImGui_GL.UI_State;
	window: Window_Type;

	exit := false;
	demo_window_open := true;

	camera_size : Physics.V2f = .{ 0.72, 0.405 };
	camera_pos  : Physics.V2f = .{ 0.0, 0.0 };
	camera_speed := 1.0;

	data_rate := 37.5;
	data_rate_changed := false;

	pause := false;
	interactive : struct {
		clear_beacon := false; @trigger;
		default_beacon := true; @trigger;

		clear_inputs := false; @trigger;
		replay_inputs___ := false;
		replay_loop   := false;

		use_beacon1 := true;
		use_beacon2 := true;

		calibrating := true;

		magnet_strength : float64 = 0.420;

		continuous := true;
		step := false;
		
		averaging_selected: Simulation.Averaging = .Avg;

		eps : float32 = 0.00001;
		n_probes := 300;

		record_inputs := false;
		clear_records := false; @trigger;

		evaluation_rate := 100;
		Sensor_Hertz :: 100;

	};
	using interactive;

	replay_input_cursor := 0;
	inputs : struct {
		mouse_pos : Physics.V2f;
		left_clicked := false;
	};
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

	region: Simulation.Region;
	beacons : [..] Data.Beacon;

	pos_samples : [..] Physics.V2f;

	time_record_loop_begin: s64;
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

pen_input :: () {
	using context.state;
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
		max_it := 1;
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
	if ui.inputs.left_clicked {
		defer ui.inputs.left_clicked = false;

		pos : Physics.V2f;
		pos.x = ui.camera_pos.x - ui.camera_size.x / 2;
		pos.x += ui.inputs.mouse_pos.x * ui.camera_size.x;
		pos.y = ui.camera_pos.y - ui.camera_size.y / 2;
		pos.y += ui.inputs.mouse_pos.y * ui.camera_size.y;

		b: Data.Beacon = .{};
		b.pos = pos;
		array_add(*beacons, b);
	}

	if ui.interactive.clear_beacon {
		defer ui.interactive.clear_beacon = false;

		beacons.count = 0;
	}

	if ui.interactive.default_beacon {
		defer ui.interactive.default_beacon = false;

		beacons.count = 0;

		b: Data.Beacon;
		b.pos = Physics.V2f.make(-0.2, -0.07);
		array_add(*beacons, b);
		b.pos = Physics.V2f.make(-0.2, +0.07);
		array_add(*beacons, b);
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
}

pen_calibrating :: (input: Interface.Input) {
	using context.state;

	if input.kind == .Magnetometer {
		min_sample_count := S64_MAX;
		for * beacons {
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

		if min_sample_count > 100 ui.interactive.calibrating = false;
	}
}

pen_update :: (input: Interface.Input) {
	using context.state;
	context_state := *context.state;
	Debug.scope_timer("update");

	if ui.interactive.calibrating {
		pen_calibrating(input);
		return;
	}
	
	if min_time > input.time min_time = input.time;

	if max_time <= input.time {
		max_time = input.time;
		t := max_time - min_time;

		for systems.count..cast(int) ceil(cast(float32) t * ui.Sensor_Hertz / 100.0) {
			system: Data.System_Reading;
			system.time = (it * 100) / ui.Sensor_Hertz + min_time;
			array_add(*systems, system);
		}
	}

	system_to_compute := NewArray(systems.count, s64, false, temp);
	system_to_compute.count = 0;
	mark_to_compute :: (idx: int) #expand {
		`system_to_compute.count += 1;
		`system_to_compute[`system_to_compute.count - 1] = idx;
	}

	for *system, it_index: systems {
		mag := *system.magnetometer_readings[input.id];

		// If the incomming input is a measurement with the exact same time we update it
		if system.time == input.time {
			mag.payload = input.mag;
			mag.valid = true;
			mag.measured = true;
			mag.interpolated_beg = input.time;
			mag.interpolated_end = input.time;
			mark_to_compute(it_index);
			continue;
		}

		// If the incomming input is within the interpolated range we update and mark it too.
		if (input.time < system.time && input.time > mag.interpolated_beg) ||
		   (input.time > system.time && input.time < mag.interpolated_end)
		{
			mag.payload = input.mag;
			mag.valid = true;
			mag.measured = false;
			if input.time > system.time
				mag.interpolated_end = min(mag.interpolated_end, input.time);
			if input.time < system.time
				mag.interpolated_beg = min(mag.interpolated_beg, input.time);
			mark_to_compute(it_index);
			continue;
		}
	}

	for idx: system_to_compute {
		it := *systems[idx];

		// if it isn't valid and the two magnetiometers reading didn't happen to flip to valid this
		// time around, then we don't continue, we wait for both reading to be true.
		if !it.valid
			if (!it.magnetometer_readings[0].valid || !it.magnetometer_readings[0].valid) continue;
			// Else we can mark it as valid
			else it.valid = true;
	}

	n := 0;
	for idx: system_to_compute if systems[idx].valid {
		it := *systems[idx];
		if n < 1 update_system(it);
		else it.valid = false;
		n += 1;
	}
}

update_system :: (system: *Data.System_Reading) {
	using context.state;

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

	for * region.matrix.elements { <<it = 1; }
	for b, b_idx: beacons {
		if !(b.std.x > 0 && b.std.y > 0 && b.std.z > 0) continue;

		mag := system.magnetometer_readings[b_idx].payload;

		Simulation.estimate_to(ui.interactive.magnet_strength, b, mag, *region);
		// Simulation.estimate_to(
		// 	ui.interactive.magnet_strength,
		// 	b,
		// 	last_mag,
		// 	*region,
		// 	pos_samples
		// );
	}

	ma : float32 = -FLOAT32_MAX;
	mi : float32 = +FLOAT32_MAX;
	for region.matrix.elements ma = max(ma, it);
	for region.matrix.elements mi = min(mi, it);
	for * region.matrix.elements { <<it /= ifx ma == 0 then 1 else ma; };

	if !ui.interactive.calibrating if ui.interactive.averaging_selected == {
		case .Avg;
			system.estimated_pos = Simulation.estimate_position_avg(region);
		case .Avg2;
			system.estimated_pos = Simulation.estimate_position_avg2(region);
	}
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
		for beacons {
			ImGui.Text(tprint("std:  % % %\0", it.std.x, it.std.y, it.std.z).data);
			ImGui.Text(tprint("mean: % % %\0", it.mean.x, it.mean.y, it.mean.z).data);
		}
		ImGui.Text(
			tprint("Inputs: Pending % Past %\0", pending_inputs.count, past_inputs.count).data
		);
		ImGui.Text(tprint("Systems: %\0", systems.count).data);
	}
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
	plot: Render.PlotConfig;
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

	Render.Plot("Magnet strength", plot);
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

	Render.matrix(region.matrix, region.space_size);
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


main :: () {
	using context.state;

	ui.window = create_window(1200, 675, "ImGui OpenGL Test");
	setup_mousemove_hook();

	gl_create_context(ui.window, 4, 5, false, true);
	gl_load(*gl);
	gl_enable_debug_output(true);

	ImGui_GL.ImGui_Impl_CreateContext(*ui.imgui, set_default_allocators=true);
	ImGui_GL.ImGui_Impl_Init(*ui.imgui, ui.window);

	current_time: float64;
	last_time: float64;

	region.space_origin = .{0, 0};
	region.space_size = .{0.3, 0.3};
	region.matrix.rows = 50;
	region.matrix.cols = 50;
	region.matrix.elements = NewArray(region.matrix.rows * region.matrix.cols, float32);
	defer free(region.matrix.elements.data);

	for * region.matrix.elements { <<it = (random_get() % 255) / 255.0; }
 
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
	user32 :: #foreign_system_library "user32";

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
}

#scope_file

USB :: #foreign_library "../../USB";
#import "GL";