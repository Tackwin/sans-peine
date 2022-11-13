#scope_file

#import "GL";
#import "Basic";
#import "Math";
#import "Shader";
Physics :: #import "Physics";
ImGui :: #import "my_imgui";

#scope_export

Attribute :: struct {
	offset: int;
	stride: int;

	make :: (offset: int, stride: int) -> Attribute {
		a: Attribute = ---;
		a.offset = offset;
		a.stride = stride;
		return a;
	}
}

quad :: (pos: Physics.V2f, size: Physics.V2f, texture: u32) {
	if !quad_gl_initted {
		defer quad_gl_initted = true;

		glGenVertexArrays(1, *quad_vao);
		glGenBuffers(1, *quad_vbo);

		quad_vertices :: float32.[
			-.5, -.5, 0, 0,
			+.5, -.5, 1, 0,
			+.5, +.5, 1, 1,
			-.5, -.5, 0, 0,
			+.5, +.5, 1, 1,
			-.5, +.5, 0, 1
		];

		glBindBuffer(GL_ARRAY_BUFFER, quad_vbo);
		glBufferData(
			GL_ARRAY_BUFFER,
			size_of(float32) * quad_vertices.count,
			quad_vertices.data,
			GL_STATIC_DRAW
		);

		textured_shader = shader_init(TEXTURED_VERT, TEXTURED_FRAG);
	}

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, texture);

	glUseProgram(textured_shader.program);
	set_uniform(textured_shader, "scaling", size);
	set_uniform(textured_shader, "offset", pos);
	set_uniform(textured_shader, "tex", cast(int) 0);

	glBindVertexArray(quad_vao);
	glBindBuffer(GL_ARRAY_BUFFER, quad_vbo);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, size_of(float32) * 4, xx 0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, size_of(float32) * 4, xx 8);
	glDrawArrays(GL_TRIANGLES, 0, 6);
}

matrix :: (mat: Physics.Matrixf, size: Physics.V2f) {
	if !quad_gl_initted {
		defer quad_gl_initted = true;

		glGenVertexArrays(1, *quad_vao);
		glGenBuffers(1, *quad_vbo);

		quad_vertices :: float32.[
			-.5, -.5, 0, 0,
			+.5, -.5, 1, 0,
			+.5, +.5, 1, 1,
			-.5, -.5, 0, 0,
			+.5, +.5, 1, 1,
			-.5, +.5, 0, 1
		];


		glBindBuffer(GL_ARRAY_BUFFER, quad_vbo);
		glBufferData(
			GL_ARRAY_BUFFER,
			size_of(float32) * quad_vertices.count,
			quad_vertices.data,
			GL_STATIC_DRAW
		);

		glGenTextures(1, *mat_tex);
		glBindTexture(GL_TEXTURE_2D, mat_tex);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

		heatmap_shader = shader_init(HEATMAP_VERT, HEATMAP_FRAG);
	}

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, mat_tex);
	glTexImage2D(
		GL_TEXTURE_2D,
		0,
		GL_R32F,
		xx mat.cols,
		xx mat.rows,
		0,
		GL_RED,
		GL_FLOAT,
		mat.elements.data
	);

	glUseProgram(heatmap_shader.program);
	set_uniform(heatmap_shader, "camera_size", context.state.ui.camera_size);
	set_uniform(heatmap_shader, "camera_pos",  context.state.ui.camera_pos);
	set_uniform(heatmap_shader, "scaling", size);

	glBindVertexArray(quad_vao);
	glBindBuffer(GL_ARRAY_BUFFER, quad_vbo);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, size_of(float32) * 4, xx 0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, size_of(float32) * 4, xx 8);
	glDrawArrays(GL_TRIANGLES, 0, 6);
}

circle :: (pos: Physics.V2f, r: float, c: Physics.V3f = .{1, 0, 0}) {
	if !circle_gl_initted {
		defer circle_gl_initted = true;

		glGenVertexArrays(1, *circle_vao);
		glGenBuffers(1, *circle_vbo);

		circle_vertices :: #run -> [20] float32 {
			res : [10 * 2] float32;
			for 0..(res.count/2)-1 {
				t := 2 * PI * it / (res.count/2.0);
				res[it * 2 + 0] = cos(t);
				res[it * 2 + 1] = sin(t);
			}
			return res;
		};

		glBindVertexArray(circle_vao);
		glBindBuffer(GL_ARRAY_BUFFER, circle_vbo);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, size_of(float32) * 2, xx 0);

		glBindBuffer(GL_ARRAY_BUFFER, circle_vbo);
		glBufferData(
			GL_ARRAY_BUFFER,
			size_of(float32) * circle_vertices.count,
			circle_vertices.data,
			GL_STATIC_DRAW
		);

		simple_shader = shader_init(SIMPLE_VERT, SIMPLE_FRAG);
	}

	glUseProgram(simple_shader.program);
	set_uniform(simple_shader, "camera_size", context.state.ui.camera_size);
	set_uniform(simple_shader, "camera_pos",  context.state.ui.camera_pos);
	set_uniform(simple_shader, "scaling", Physics.V2f.make(r, r));
	set_uniform(simple_shader, "color", c);
	set_uniform(simple_shader, "offset", pos);

	glBindVertexArray(circle_vao);
	glBindBuffer(GL_ARRAY_BUFFER, circle_vbo);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, size_of(float32) * 2, xx 0);

	glBindVertexArray(circle_vao);
	glDrawArrays(GL_TRIANGLE_FAN, 0, 10);
}

circles :: (pos: [] Physics.V2f, r: [] float, colors: [] Physics.V3f)
{
	if !circles_gl_initted {
		defer circles_gl_initted = true;

		glGenVertexArrays(1, *circles_vao);
		glGenBuffers(1, *circles_vbo);
		glGenBuffers(1, *circles_ibo);

		circles_vertices :: #run -> [20] float32 {
			res : [10 * 2] float32;
			for 0..(res.count/2)-1 {
				t := 2 * PI * it / (res.count/2.0);
				res[it * 2 + 0] = cos(t);
				res[it * 2 + 1] = sin(t);
			}
			return res;
		};


		glBindBuffer(GL_ARRAY_BUFFER, circles_vbo);
		glBufferData(
			GL_ARRAY_BUFFER,
			size_of(float32) * circles_vertices.count,
			circles_vertices.data,
			GL_STATIC_DRAW
		);

		simple_instanced_shader = shader_init(SIMPLE_INSTANCED_VERT, SIMPLE_INSTANCED_FRAG);
	}

	buffer := NewArray(pos.count * size_of(float) * 7, u8, false, temp);
	view := cast(*float) buffer.data;
	cursor := 0;
	for pos {
		view[cursor] = pos[it_index].x; cursor += 1;
		view[cursor] = pos[it_index].y; cursor += 1;
		view[cursor] = r[it_index]; cursor += 1;
		view[cursor] = r[it_index]; cursor += 1;
		view[cursor] = colors[it_index].x; cursor += 1;
		view[cursor] = colors[it_index].y; cursor += 1;
		view[cursor] = colors[it_index].z; cursor += 1;
	}

	glBindBuffer(GL_ARRAY_BUFFER, circles_ibo);
	glBufferData(GL_ARRAY_BUFFER, buffer.count, buffer.data, GL_DYNAMIC_DRAW);

	glUseProgram(simple_instanced_shader.program);
	set_uniform(simple_instanced_shader, "camera_size", context.state.ui.camera_size);
	set_uniform(simple_instanced_shader, "camera_pos",  context.state.ui.camera_pos);

	glBindVertexArray(circles_vao);
	glBindBuffer(GL_ARRAY_BUFFER, circles_vbo);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, size_of(float32) * 2, xx 0);
	glVertexAttribDivisor(0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, circles_ibo);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, size_of(float32) * 7, xx 0);
	glVertexAttribDivisor(1, 1);

	glEnableVertexAttribArray(2);
	glVertexAttribPointer(
		2, 2, GL_FLOAT, GL_FALSE, size_of(float32) * 7, xx (size_of(float32) * 2)
	);
	glVertexAttribDivisor(2, 1);

	glEnableVertexAttribArray(3);
	glVertexAttribPointer(
		3, 3, GL_FLOAT, GL_FALSE, size_of(float32) * 7, xx (size_of(float32) * 4)
	);
	glVertexAttribDivisor(3, 1);

	glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, 10, xx pos.count);
}

polyline :: (
	n: int, buffer: [] u8, pos: Attribute, thickness: Attribute, colors: Attribute
) {
	if !polyline_gl_initted {
		defer polyline_gl_initted = true;

		glGenVertexArrays(1, *polyline_vao);
		glGenBuffers(1, *polyline_vbo);

		glBindVertexArray(polyline_vao);
		glBindBuffer(GL_ARRAY_BUFFER, polyline_vbo);

		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, size_of(float32) * 12, xx 0);
		glVertexAttribDivisor(0, 0);

		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, size_of(float32) * 12, xx 8);
		glVertexAttribDivisor(1, 0);
		
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, size_of(float32) * 12, xx 24);
		glVertexAttribDivisor(2, 0);
		
		glEnableVertexAttribArray(3);
		glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, size_of(float32) * 12, xx 32);
		glVertexAttribDivisor(3, 0);

		glEnableVertexAttribArray(4);
		glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, size_of(float32) * 12, xx 36);
		glVertexAttribDivisor(4, 0);

		polyline_shader = shader_init(POLYLINE_VERT, POLYLINE_FRAG);
	}

	glBindVertexArray(polyline_vao);
	glBindBuffer(GL_ARRAY_BUFFER, polyline_vbo);

	end_buffer := NewArray(n * 16 * 2, float32, false, temp);

	cursor := 0;
	prev := Physics.V2f.{0, 0};
	curr := Physics.V2f.{0, 0};
	s0 : float64 = 0;

	for 0..n {
		next := curr;
		defer {
			prev = curr;
			curr = next;
		}

		pos_cursor := pos.offset + it * pos.stride;

		if it < n next = Physics.V2f.make(
			(cast(*float32) (buffer.data + pos_cursor))[0],
			(cast(*float32) (buffer.data + pos_cursor))[1]
		);

		if it == 0 continue;
		if it == n next = curr;

		thick_cursor := thickness.offset + (it - 1) * thickness.stride;
		color_cursor := colors.offset + (it - 1) * colors.stride;

		using Physics;
		defer s0 += cast(float64) Physics.length(curr - prev);

		end_buffer[cursor] = prev.x; cursor += 1;
		end_buffer[cursor] = prev.y; cursor += 1;
		end_buffer[cursor] = curr.x; cursor += 1;
		end_buffer[cursor] = curr.y; cursor += 1;
		end_buffer[cursor] = +1; cursor += 1;
		end_buffer[cursor] = cast(float) s0; cursor += 1;
		end_buffer[cursor] = next.x; cursor += 1;
		end_buffer[cursor] = next.y; cursor += 1;
		end_buffer[cursor] = (cast (*float32) (buffer.data + thick_cursor))[0]; cursor += 1;
		end_buffer[cursor] = (cast (*float32) (buffer.data + color_cursor))[0]; cursor += 1;
		end_buffer[cursor] = (cast (*float32) (buffer.data + color_cursor))[1]; cursor += 1;
		end_buffer[cursor] = (cast (*float32) (buffer.data + color_cursor))[2]; cursor += 1;

		end_buffer[cursor] = prev.x; cursor += 1;
		end_buffer[cursor] = prev.y; cursor += 1;
		end_buffer[cursor] = curr.x; cursor += 1;
		end_buffer[cursor] = curr.y; cursor += 1;
		end_buffer[cursor] = -1; cursor += 1;
		end_buffer[cursor] = cast(float) s0; cursor += 1;
		end_buffer[cursor] = next.x; cursor += 1;
		end_buffer[cursor] = next.y; cursor += 1;
		end_buffer[cursor] = (cast (*float32) (buffer.data + thick_cursor))[0]; cursor += 1;
		end_buffer[cursor] = (cast (*float32) (buffer.data + color_cursor))[0]; cursor += 1;
		end_buffer[cursor] = (cast (*float32) (buffer.data + color_cursor))[1]; cursor += 1;
		end_buffer[cursor] = (cast (*float32) (buffer.data + color_cursor))[2]; cursor += 1;
	}

	glBufferData(
		GL_ARRAY_BUFFER, size_of(float32) * end_buffer.count, end_buffer.data, GL_DYNAMIC_DRAW
	);

	glUseProgram(polyline_shader.program);
	set_uniform(polyline_shader, "camera_size", context.state.ui.camera_size);
	set_uniform(polyline_shader, "camera_pos",  context.state.ui.camera_pos);
	set_uniform(polyline_shader, "antialias",   0.0001);
	set_uniform(polyline_shader, "arclength",   s0);

	glDrawArrays(GL_TRIANGLE_STRIP, 0, xx (2 * n));
}

render_to_framebuffer :: () {
	if !framebuffer_initted {
		defer framebuffer_initted = true;

		glGenFramebuffers(1, *framebuffer);
		glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
		glGenTextures(1, *framebuffer_texture);
		glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, framebuffer_texture);
		glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, 16, GL_RGBA, 1920, 1027, GL_TRUE);

		glFramebufferTexture2D(
			GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, framebuffer_texture, 0
		);

		textured_shader = shader_init(TEXTURED_VERT, TEXTURED_FRAG);
	}
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
}
render_to_screen :: () {
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

framebuffer_texture: u32 = 0;

magma_colormap :: (t: float) -> Physics.V3f {
	_magma_table :: Physics.V3f.[
		.{0.001, 0.000, 0.013},.{0.002, 0.001, 0.018},.{0.003, 0.002, 0.023},.{0.004, 0.003, 0.029},
		.{0.005, 0.004, 0.037},.{0.007, 0.006, 0.044},.{0.009, 0.008, 0.052},.{0.011, 0.009, 0.060},
		.{0.013, 0.011, 0.068},.{0.016, 0.013, 0.076},.{0.018, 0.016, 0.084},.{0.021, 0.018, 0.092},
		.{0.024, 0.020, 0.100},.{0.028, 0.023, 0.108},.{0.031, 0.025, 0.116},.{0.035, 0.028, 0.125},
		.{0.039, 0.031, 0.133},.{0.043, 0.033, 0.141},.{0.048, 0.036, 0.150},.{0.052, 0.039, 0.158},
		.{0.056, 0.042, 0.167},.{0.060, 0.044, 0.176},.{0.065, 0.047, 0.184},.{0.069, 0.049, 0.193},
		.{0.074, 0.052, 0.202},.{0.078, 0.054, 0.211},.{0.083, 0.056, 0.220},.{0.088, 0.058, 0.229},
		.{0.092, 0.059, 0.239},.{0.097, 0.061, 0.248},.{0.102, 0.063, 0.257},.{0.107, 0.064, 0.267},
		.{0.113, 0.065, 0.276},.{0.118, 0.066, 0.286},.{0.123, 0.067, 0.295},.{0.129, 0.067, 0.305},
		.{0.135, 0.068, 0.315},.{0.140, 0.068, 0.324},.{0.146, 0.068, 0.334},.{0.152, 0.068, 0.343},
		.{0.159, 0.068, 0.352},.{0.165, 0.067, 0.361},.{0.171, 0.067, 0.370},.{0.178, 0.066, 0.379},
		.{0.184, 0.065, 0.387},.{0.191, 0.064, 0.396},.{0.198, 0.063, 0.404},.{0.204, 0.062, 0.411},
		.{0.211, 0.061, 0.418},.{0.218, 0.061, 0.425},.{0.225, 0.060, 0.431},.{0.232, 0.059, 0.437},
		.{0.238, 0.059, 0.443},.{0.245, 0.059, 0.448},.{0.252, 0.059, 0.453},.{0.258, 0.059, 0.457},
		.{0.265, 0.060, 0.461},.{0.271, 0.060, 0.465},.{0.278, 0.061, 0.469},.{0.284, 0.063, 0.472},
		.{0.291, 0.064, 0.475},.{0.297, 0.066, 0.478},.{0.304, 0.067, 0.480},.{0.310, 0.069, 0.483},
		.{0.316, 0.071, 0.485},.{0.322, 0.073, 0.487},.{0.329, 0.075, 0.489},.{0.335, 0.078, 0.491},
		.{0.341, 0.080, 0.492},.{0.347, 0.082, 0.494},.{0.353, 0.085, 0.495},.{0.359, 0.087, 0.496},
		.{0.366, 0.090, 0.497},.{0.372, 0.092, 0.499},.{0.378, 0.095, 0.500},.{0.384, 0.097, 0.501},
		.{0.390, 0.100, 0.501},.{0.396, 0.102, 0.502},.{0.402, 0.105, 0.503},.{0.408, 0.107, 0.504},
		.{0.414, 0.110, 0.504},.{0.420, 0.112, 0.505},.{0.426, 0.115, 0.505},.{0.432, 0.117, 0.506},
		.{0.439, 0.120, 0.506},.{0.445, 0.122, 0.506},.{0.451, 0.125, 0.507},.{0.457, 0.127, 0.507},
		.{0.463, 0.129, 0.507},.{0.469, 0.132, 0.507},.{0.475, 0.134, 0.507},.{0.481, 0.136, 0.507},
		.{0.488, 0.139, 0.508},.{0.494, 0.141, 0.507},.{0.500, 0.143, 0.507},.{0.506, 0.145, 0.507},
		.{0.512, 0.148, 0.507},.{0.519, 0.150, 0.507},.{0.525, 0.152, 0.507},.{0.531, 0.154, 0.506},
		.{0.537, 0.156, 0.506},.{0.544, 0.159, 0.506},.{0.550, 0.161, 0.505},.{0.556, 0.163, 0.505},
		.{0.562, 0.165, 0.504},.{0.569, 0.167, 0.504},.{0.575, 0.169, 0.503},.{0.581, 0.171, 0.502},
		.{0.588, 0.173, 0.502},.{0.594, 0.175, 0.501},.{0.600, 0.177, 0.500},.{0.607, 0.179, 0.499},
		.{0.613, 0.181, 0.498},.{0.620, 0.183, 0.497},.{0.626, 0.185, 0.496},.{0.632, 0.187, 0.495},
		.{0.639, 0.189, 0.494},.{0.645, 0.191, 0.492},.{0.652, 0.193, 0.491},.{0.658, 0.196, 0.490},
		.{0.664, 0.198, 0.488},.{0.671, 0.200, 0.487},.{0.677, 0.202, 0.485},.{0.684, 0.204, 0.484},
		.{0.690, 0.206, 0.482},.{0.697, 0.208, 0.480},.{0.703, 0.210, 0.479},.{0.709, 0.212, 0.477},
		.{0.716, 0.214, 0.475},.{0.722, 0.217, 0.473},.{0.729, 0.219, 0.471},.{0.735, 0.221, 0.469},
		.{0.742, 0.224, 0.467},.{0.748, 0.226, 0.464},.{0.754, 0.228, 0.462},.{0.761, 0.231, 0.460},
		.{0.767, 0.233, 0.457},.{0.773, 0.236, 0.455},.{0.779, 0.238, 0.452},.{0.786, 0.241, 0.450},
		.{0.792, 0.244, 0.447},.{0.798, 0.247, 0.444},.{0.804, 0.249, 0.442},.{0.810, 0.252, 0.439},
		.{0.816, 0.255, 0.436},.{0.822, 0.259, 0.433},.{0.828, 0.262, 0.430},.{0.834, 0.265, 0.427},
		.{0.840, 0.268, 0.424},.{0.846, 0.272, 0.421},.{0.852, 0.276, 0.418},.{0.857, 0.279, 0.415},
		.{0.863, 0.283, 0.412},.{0.868, 0.287, 0.409},.{0.874, 0.291, 0.406},.{0.879, 0.296, 0.403},
		.{0.884, 0.300, 0.400},.{0.889, 0.305, 0.397},.{0.894, 0.309, 0.393},.{0.899, 0.314, 0.391},
		.{0.904, 0.319, 0.388},.{0.908, 0.324, 0.385},.{0.913, 0.330, 0.382},.{0.917, 0.335, 0.379},
		.{0.921, 0.341, 0.377},.{0.925, 0.346, 0.374},.{0.929, 0.352, 0.372},.{0.933, 0.358, 0.370},
		.{0.937, 0.364, 0.368},.{0.940, 0.371, 0.366},.{0.944, 0.377, 0.365},.{0.947, 0.384, 0.363},
		.{0.950, 0.390, 0.362},.{0.953, 0.397, 0.361},.{0.955, 0.404, 0.360},.{0.958, 0.411, 0.360},
		.{0.960, 0.418, 0.359},.{0.963, 0.425, 0.359},.{0.965, 0.432, 0.359},.{0.967, 0.439, 0.359},
		.{0.969, 0.446, 0.360},.{0.971, 0.454, 0.361},.{0.973, 0.461, 0.361},.{0.975, 0.468, 0.363},
		.{0.976, 0.476, 0.364},.{0.978, 0.483, 0.366},.{0.979, 0.491, 0.367},.{0.981, 0.498, 0.369},
		.{0.982, 0.505, 0.371},.{0.983, 0.513, 0.374},.{0.984, 0.520, 0.376},.{0.985, 0.528, 0.379},
		.{0.986, 0.535, 0.382},.{0.987, 0.543, 0.385},.{0.988, 0.550, 0.388},.{0.989, 0.557, 0.391},
		.{0.990, 0.565, 0.395},.{0.990, 0.572, 0.398},.{0.991, 0.580, 0.402},.{0.992, 0.587, 0.406},
		.{0.992, 0.594, 0.410},.{0.993, 0.602, 0.414},.{0.993, 0.609, 0.418},.{0.994, 0.616, 0.422},
		.{0.994, 0.624, 0.427},.{0.995, 0.631, 0.431},.{0.995, 0.639, 0.436},.{0.995, 0.646, 0.441},
		.{0.996, 0.653, 0.446},.{0.996, 0.660, 0.451},.{0.996, 0.668, 0.456},.{0.996, 0.675, 0.461},
		.{0.996, 0.682, 0.466},.{0.997, 0.690, 0.471},.{0.997, 0.697, 0.477},.{0.997, 0.704, 0.482},
		.{0.997, 0.711, 0.488},.{0.997, 0.719, 0.493},.{0.997, 0.726, 0.499},.{0.997, 0.733, 0.505},
		.{0.997, 0.740, 0.510},.{0.997, 0.747, 0.516},.{0.997, 0.755, 0.522},.{0.997, 0.762, 0.528},
		.{0.996, 0.769, 0.534},.{0.996, 0.776, 0.541},.{0.996, 0.783, 0.547},.{0.996, 0.791, 0.553},
		.{0.996, 0.798, 0.559},.{0.995, 0.805, 0.566},.{0.995, 0.812, 0.572},.{0.995, 0.819, 0.579},
		.{0.995, 0.827, 0.585},.{0.994, 0.834, 0.592},.{0.994, 0.841, 0.598},.{0.994, 0.848, 0.605},
		.{0.993, 0.855, 0.612},.{0.993, 0.862, 0.619},.{0.993, 0.870, 0.626},.{0.992, 0.877, 0.633},
		.{0.992, 0.884, 0.640},.{0.992, 0.891, 0.647},.{0.991, 0.898, 0.654},.{0.991, 0.905, 0.661},
		.{0.990, 0.912, 0.668},.{0.990, 0.920, 0.675},.{0.990, 0.927, 0.682},.{0.989, 0.934, 0.690},
		.{0.989, 0.941, 0.697},.{0.989, 0.948, 0.704},.{0.988, 0.955, 0.712},.{0.988, 0.962, 0.719},
		.{0.988, 0.970, 0.727},.{0.987, 0.977, 0.734},.{0.987, 0.984, 0.742},.{0.987, 0.991, 0.749},
	];

	return _magma_table[clamp(cast(int) (255 * t), 0, 255)];
}

#scope_file

framebuffer_initted := false;
framebuffer: u32 = 0;
textured_shader: Shader;

quad_gl_initted := false;
quad_vao : u32 = 0;
quad_vbo : u32 = 0;
quad_tex : u32;
mat_tex : u32;
heatmap_shader : Shader;


circle_gl_initted := false;
circle_vao : u32 = 0;
circle_vbo : u32 = 0;
circle_tex : u32;
simple_shader : Shader;

circles_gl_initted := false;
circles_vao : u32 = 0;
circles_vbo : u32 = 0;
circles_ibo : u32 = 0;
circles_instance_data_cpu : [..] u8;
simple_instanced_shader : Shader;

polyline_gl_initted := false;
polyline_vao : u32 = 0;
polyline_vbo : u32 = 0;
polyline_shader: Shader;

HEATMAP_FRAG :: #string FRAG
#version 410 core
in vec2 frag_uv;

uniform sampler2D heatmap;

out vec4 frag_color;

vec3 viridis(float t) {
	const vec3 c0 = vec3(0.2777273272234177, 0.005407344544966578, 0.3340998053353061);
	const vec3 c1 = vec3(0.1050930431085774, 1.404613529898575, 1.384590162594685);
	const vec3 c2 = vec3(-0.3308618287255563, 0.214847559468213, 0.09509516302823659);
	const vec3 c3 = vec3(-4.634230498983486, -5.799100973351585, -19.33244095627987);
	const vec3 c4 = vec3(6.228269936347081, 14.17993336680509, 56.69055260068105);
	const vec3 c5 = vec3(4.776384997670288, -13.74514537774601, -65.35303263337234);
	const vec3 c6 = vec3(-5.435455855934631, 4.645852612178535, 26.3124352495832);

	return c0+t*(c1+t*(c2+t*(c3+t*(c4+t*(c5+t*c6)))));
}

void main() {
	frag_color = vec4(viridis(texture(heatmap, frag_uv).r), 1);
}
FRAG

HEATMAP_VERT :: #string VERT
#version 410 core
layout (location = 0) in vec2 pos;
layout (location = 1) in vec2 uv;

uniform vec2 camera_size;
uniform vec2 camera_pos;

uniform vec2 scaling;

out vec2 frag_uv;

void main() {
	frag_uv = uv;
	gl_Position = vec4(
		(scaling.x * pos.x - camera_pos.x) / camera_size.x,
		(scaling.y * pos.y - camera_pos.y) / camera_size.y,
		0,
		1.0
	);
	gl_Position.x *= 2;
	gl_Position.y *= 2;
}
VERT
SIMPLE_FRAG :: #string FRAG
#version 410 core

uniform vec3 color;

out vec4 frag_color;

void main() {
	frag_color = vec4(color, 1);
}
FRAG

SIMPLE_VERT :: #string VERT
#version 410 core
layout (location = 0) in vec2 pos;

uniform vec2 camera_size;
uniform vec2 camera_pos;
uniform vec2 scaling;
uniform vec2 offset;

void main() {
	gl_Position = vec4(
		(scaling.x * pos.x - camera_pos.x + offset.x) / camera_size.x,
		(scaling.y * pos.y - camera_pos.y + offset.y) / camera_size.y,
		0,
		1.0
	);
	gl_Position.x *= 2;
	gl_Position.y *= 2;
}
VERT
TEXTURED_FRAG :: #string FRAG
#version 410 core
layout (location = 0) in vec2 uv;

uniform sampler2DMS tex;

out vec4 frag_color;

vec4 lookup(sampler2DMS t, vec2 uv) {
	vec4 sum = vec4(0);
	ivec2 size = textureSize(t);
	for (int i = 0; i < 16; i += 1) {
		sum += texelFetch(t, ivec2(size.x * uv.x, size.y * uv.y), i);
	}
	return sum / 16;
}

void main() {
	frag_color = lookup(tex, uv);
}
FRAG

TEXTURED_VERT :: #string VERT
#version 410 core
layout (location = 0) in vec2 pos;
layout (location = 1) in vec2 uv;

layout (location = 0) out vec2 out_uv;

uniform vec2 camera_size;
uniform vec2 camera_pos;
uniform vec2 scaling;
uniform vec2 offset;

void main() {
	gl_Position = vec4(
		scaling.x * pos.x + offset.x,
		scaling.y * pos.y + offset.y,
		0,
		1.0
	);
	out_uv = uv;
}
VERT

SIMPLE_INSTANCED_FRAG :: #string FRAG
#version 410 core

in vec3 color;

out vec4 frag_color;

void main() {
	frag_color = vec4(1);
}
FRAG

SIMPLE_INSTANCED_VERT :: #string VERT
#version 410 core
layout (location = 0) in vec2 pos;
layout (location = 1) in vec2 offset;
layout (location = 2) in vec2 scaling;
layout (location = 3) in vec3 in_color;

uniform vec2 camera_size;
uniform vec2 camera_pos;

out vec3 out_color;

void main() {
	gl_Position = vec4(
		(scaling.x * pos.x - camera_pos.x + offset.x) / camera_size.x,
		(scaling.y * pos.y - camera_pos.y + offset.y) / camera_size.y,
		0.5,
		1.0
	);
	gl_Position.x *= 2;
	gl_Position.y *= 2;
	out_color = in_color;
}
VERT

POLYLINE_FRAG :: #string FRAG
#version 410 core
layout (location = 0) in float in_thickness;
layout (location = 1) in vec2  in_uv;
layout (location = 2) in vec3  in_color;

out vec4 out_color;

uniform float arclength;
uniform float antialias;

void main() {
	float w = in_thickness / 2 - antialias;

	float d = 0;
	     if (in_uv.x < 0)          d = length(in_uv) - w;
	else if (in_uv.x >= arclength) d = length(in_uv - vec2(arclength, 0)) - w;
	else                           d = abs(in_uv.y) - w;

	d = max(d, 0);
	d /= antialias;
	out_color = vec4(in_color, exp(-d*d));
}
FRAG

POLYLINE_VERT :: #string VERT
#version 410 core
layout (location = 0) in vec2  prev;
layout (location = 1) in vec4  curr;
layout (location = 2) in vec2  next;
layout (location = 3) in float curr_thick;
layout (location = 4) in vec3  curr_color;

layout (location = 0) out float out_thickness;
layout (location = 1) out vec2  out_uv;
layout (location = 2) out vec3  out_color;

uniform vec2 camera_size;
uniform vec2 camera_pos;
uniform float arclength;
uniform float antialias;

void main() {
	out_thickness = curr_thick;
	out_color     = curr_color;

	float w = curr_thick / 2 + antialias;
	
	vec2 t0 = normalize(curr.xy - prev.xy);
	if (curr.xy == prev.xy) t0 = vec2(1, 0);
	vec2 t1 = normalize(next.xy - curr.xy);
	if (curr.xy == next.xy) t1 = vec2(1, 0);
	vec2 n0 = vec2(-t0.y, t0.x);
	vec2 n1 = vec2(-t1.y, t1.x);

	vec2 uv;
	vec2 p;

	if (prev.xy == curr.xy) {
		uv = vec2(-w, curr.z * w);
		p = curr.xy - w * t1 + curr.z * w * n1;
	}
	else if (curr.xy == next.xy) {
		uv = vec2(arclength + w, curr.z * w);
		p = curr.xy + w * t0  + curr.z * w * n0;
	}
	else {
		vec2 m = normalize(n0 + n1);
		float dy = w / max(dot(m, n1), 0.5);
		uv = vec2(curr.w, curr.z * w);
		p = curr.xy + dy * curr.z * m;
	}

	out_uv = uv;
	gl_Position = vec4(
		2 * (p.x - camera_pos.x) / camera_size.x,
		2 * (p.y - camera_pos.y) / camera_size.y,
		0.0,
		1.0
	);
}
VERT
