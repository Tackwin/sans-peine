#scope_file
ImGui :: #import "my_imgui";
ImGui_GL :: #import "my_imgui/imgui_impl_gl";
Render :: #import "Render";
Interface :: #import "Interface";
Physics :: #import "Physics";
#import "Random";
#import "Reading";

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

	camera_size : V2d = .{ 0.8, 0.45 };
	camera_pos  : V2d = .{ 0.0, 0.0 };
	camera_speed := 1.0;
};

State :: struct {
	current_driver: Interface.Driver;
	ui: Ui;
	inputs: [..] Interface.Input;
	
	dt: float64;

	space_probability: Physics.Matrixf;
}

#add_context state: State;

input :: () {
	using context.state;

	io := ImGui.GetIO();
	for events_this_frame {
		ImGui_GL.ImGui_Impl_ProcessEvent(*ui.imgui, it);
		if it.type == {
			case .QUIT; ui.exit = true;
			case .KEYBOARD;
				if io.WantCaptureKeyboard continue;
				if it.key_pressed == 0 continue;
				if it.key_code == .ESCAPE ui.exit = true;
				if it.key_code == #char "Q" ui.camera_pos.x -= ui.camera_speed * dt;
				if it.key_code == #char "D" ui.camera_pos.x += ui.camera_speed * dt;
				if it.key_code == #char "Z" ui.camera_pos.y += ui.camera_speed * dt;
				if it.key_code == #char "S" ui.camera_pos.y -= ui.camera_speed * dt;
		}
	}
}


update :: () {
	using context.state;

	handle_input :: (input: Interface.Input) {
		using context.state;

		array_add(*inputs, input);
	}

	max_it := 10;
	for 1..max_it {
		new_input, has_new_input := Interface.next(current_driver);
		if has_new_input handle_input(new_input);
		else break;
	}
}


render :: () {
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

	if ui.demo_window_open
		ImGui.ShowDemoWindow(*ui.demo_window_open);
	ImGui.Begin("Test");
	defer ImGui.End();

	// if ImGui.CollapsingHeader("Context")
	// 	ImGui.inspect("State", *context.state);
	
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
	for y: ys for y plot.scale.min = min(plot.scale.min, it - 1);
	for y: ys for y plot.scale.max = max(plot.scale.max, it + 1);

	Render.Plot("Magnet strength", plot);
}


main :: () {
	using context.state;

	current_driver = load_driver("USB.dll");
	current_driver.memory = alloc(Interface.size(current_driver));
	defer free(current_driver.memory);

	init_res := Interface.init(current_driver);
	if !init_res.success {
		print("Error while opening driver %\n", init_res);
		return;
	}
	defer Interface.shut(current_driver);

	Interface.play(current_driver);

	ui.window = create_window(1200, 675, "ImGui OpenGL Test");
	setup_mousemove_hook();

	gl_create_context(ui.window, 3, 2);
	gl_load(*gl);

	ImGui_GL.ImGui_Impl_CreateContext(*ui.imgui, set_default_allocators=true);
	ImGui.StyleColorsLight();
	ImGui_GL.ImGui_Impl_Init(*ui.imgui, ui.window);

	current_time: float64;
	last_time: float64;

	space_probability.rows = 100;
	space_probability.cols = 100;
	space_probability.elements = NewArray(
		space_probability.rows * space_probability.cols,
		float32
	);
	defer free(space_probability.elements.data);

	for * space_probability.elements { <<it = (random_get() % 255) / 255.0; }

	while !ui.exit {
		at := current_time_monotonic();
		current_time = to_float64_seconds(at);
		context.state.dt = cast(float)(current_time - last_time);
		defer last_time = current_time;

		update_window_events();
		input();

		update();
		
		w, h := get_window_size(ui.window);

		ImGui_GL.ImGui_Impl_NewFrame(*ui.imgui, ui.window, w, h);
		render();
		ImGui.Render();

		glViewport(0, 0, xx w, xx h);
		glClearColor(0, 0, 0, 1);
		glClear(GL_COLOR_BUFFER_BIT);

		Render.render(space_probability, .{});

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
		driver.shut = xx GetProcAddress(dll, "shit");
		driver.play = xx GetProcAddress(dll, "play");
		driver.stop = xx GetProcAddress(dll, "stop");
		driver.next = xx GetProcAddress(dll, "next");
		return driver;
	}
}

#scope_file

USB :: #foreign_library "../../USB";