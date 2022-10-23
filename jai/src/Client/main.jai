#scope_file
ImGui :: #import "my_imgui";
ImGui_GL :: #import "my_imgui/imgui_impl_gl";
Render :: #import "Render";
Interface :: #import "Interface";
Physics :: #import "Physics";
Data :: #import "Data";
Simulation :: #import "Simulation";
#import "Random";

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

	camera_size : Physics.V2f = .{ 0.8, 0.45 };
	camera_pos  : Physics.V2f = .{ 0.0, 0.0 };
	camera_speed := 1.0;

	pause := false;
	interactive : struct {
		clear_beacon := false; @trigger;
		default_beacon := false; @trigger;

		clear_inputs := false; @trigger;

		use_beacon1 := true;
		use_beacon2 := true;

		calibrating := true;

		magnet_strength : float64 = 1;
		
		averaging_selected: Simulation.Averaging = .Avg;
	};
	using interactive;

	inputs : struct {
		mouse_pos : Physics.V2f;
		left_clicked := false;
	};
};

State :: struct {
	current_driver: Interface.Driver;
	ui: Ui;

	inputs: [..] Interface.Input;
	
	dt: float64;
	dt32: float32;

	region: Simulation.Region;
	beacons : [..] Data.Beacon;

	computed_path : Data.Path;
}

#add_context state: State;

input :: () {
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


update :: () {
	using context.state;

	handle_input :: (input: Interface.Input) {
		using context.state;
		if ui.pause return;

		array_add(*inputs, input);

		if input.kind == .Magnetometer && ui.interactive.calibrating {
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
			}
		}
	}

	max_it := 10;
	for 1..max_it {
		new_input, has_new_input := Interface.next(current_driver);
		if has_new_input handle_input(new_input);
		else break;
	}

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
		computed_path.points.count = 0;
		inputs.count = 0;
	}

	cdf := Physics.copy(region.matrix);
	{
		sum := 0.0;
		for cdf.elements sum += it;
		for * cdf.elements { <<it /= sum; }
		for * cdf.elements { <<it += 0.02 / (cdf.rows * cdf.cols); }
	}
	pos    := Physics.sample_from(cdf, 100, temp);

	for * region.matrix.elements { <<it = 1; }
	for b, b_idx: beacons {
		last_mag: Data.Magnetometer;
		for < inputs if it.id == b_idx { last_mag = it.mag; break; }

		Simulation.estimate_to(ui.interactive.magnet_strength, b, last_mag, *region);
		// Simulation.estimate_to(
		// 	ui.interactive.magnet_strength,
		// 	b,
		// 	last_mag,
		// 	*region,
		// 	pos
		// );
	}

	ma : float32 = -FLOAT32_MAX;
	mi : float32 = +FLOAT32_MAX;
	for region.matrix.elements ma = max(ma, it);
	for region.matrix.elements mi = min(mi, it);
	for * region.matrix.elements { <<it /= ifx ma == 0 then 1 else ma; };

	if !ui.interactive.calibrating if ui.interactive.averaging_selected == {
		case .Avg;
			array_add(*computed_path.points, Simulation.estimate_position(region));
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

	if ImGui.CollapsingHeader("Interactive")
		ImGui.inspect("State", *context.state.ui.interactive);
	
	ys: [2][]float;
	for * ys { <<it = NewArray(inputs.count, float, false, temp); }

	idx1 := 0;
	idx2 := 0;
	for inputs if it.kind == .Magnetometer {
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
	plot.frame_size.x = 100;
	plot.frame_size.y = 100;

	plot.scale.min = FLOAT32_MAX;
	plot.scale.max = -FLOAT32_MAX;
	for y: ys for y plot.scale.min = min(plot.scale.min, it - 0.000_001);
	for y: ys for y plot.scale.max = max(plot.scale.max, it + 0.000_001);

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

	Vertex :: struct {
		pos: Physics.V2f #align 1;
		thick : float32 = 0.003 #align 1;
		color := Physics.V3f.{1, 1, 1} #align 1;
	}
	vertices := NewArray(computed_path.points.count, Vertex, false, temp);
	for * vertices {
		it.pos = computed_path.points[it_index];
		it.thick = 0.001;
		it.color = Render.magma_colormap(it_index / (vertices.count - 1.0));
	}
	verticesu8: [] u8;
	verticesu8.data = xx vertices.data;
	verticesu8.count = size_of(Vertex) * vertices.count;
	Render.polyline(
		computed_path.points.count,
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
	region.space_size = .{0.5, 0.5};
	region.matrix.rows = 300;
	region.matrix.cols = 300;
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

	while !ui.exit {
		at := current_time_monotonic();
		current_time = to_float64_seconds(at);
		context.state.dt = cast(float64)(current_time - last_time);
		context.state.dt32 = cast(float32)(current_time - last_time);
		defer last_time = current_time;

		update_window_events();
		input();

		update();
		
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