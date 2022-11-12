#scope_file

#import "GL";
#import "Basic";
#import "Math";
#import "Shader";
Physics :: #import "Physics";
ImGui :: #import "my_imgui";

#scope_export


// Use this structure to pass the plot data and settings into the Plot function
PlotConfig :: struct {
	Values :: struct {
		// if necessary, you can provide x-axis values
		xs : *float = null;
		// array of y values. If null, use ys_list (below)
		ys : *float = null;
		// the number of values in each array
		count := 0;
		// at which offset to start plotting.
		// Warning: count+offset must be <= length of array!
		offset := 0;
		// Plot color. If 0, use ImGuiCol_PlotLines.
		color : u32 = 0;

		// in case you need to draw multiple plots at once, use this instead of ys
		ys_list: **float = null;
		// the number of plots to draw
		ys_count := 0;
		// colors for each plot
		colors : *u32 = null;
	}
	values: Values;
	Scale :: struct {
		// Minimum plot value
		min: float;
		// Maximum plot value
		max: float;
		Kind :: enum {
			Linear :: 0;
			Log10 :: 1;
		}
		// How to scale the x-axis
		type : Kind = .Linear;
	}
	scale: Scale;

	Tooltip :: struct {
		show := false;
		ys_names : *string = null;
		format := "%g: %8.4g";
	}
	tooltip : Tooltip;

	Grid :: struct {
		show := false;
		size := 100.0; // at which intervals to draw the grid
		subticks := 10; // how many subticks in each tick
	}
	grid_x : Grid;
	grid_y : Grid;

	Selection :: struct {
		show := false;
		start : *u32 = null;
		length : *u32 = null;
		// "Sanitize" function. Give it selection length, and it will return
		// the "allowed" length. Useful for FFT, where selection must be
		// of power of two
		sanitize_fn : (u32) -> u32 = null;
	}
	selection : Selection;

	VerticalLines :: struct {
		show := false;
		indices : *int = null; // at which indices to draw the lines
		count := 0;
	}
	v_lines: VerticalLines;

	frame_size := ImGui.ImVec2.{0, 0};
	line_thickness := 1.0;
	skip_small_lines := true;
	overlay_text: string;
}

PlotStatus :: enum {
	nothing;
	selection_updated;
}



#scope_file
using ImGui;
log10 :: log;

ImLerp :: (a: $T, b: T, t: T) -> T { return a + (b - a) * t; }
ImLerp :: (a: ImVec2, b: ImVec2, t: float) -> ImVec2 {
	return make_ImVec2(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t);
}
ImSaturate :: (x: float) -> float { return ifx x < 0 then 0 else (ifx x < 1 then x else 1); }
ImMin :: (a: $T, b: T) -> T { if a < b return a; else return b; }
ImMax :: (a: $T, b: T) -> T { if a > b return a; else return b; }

Contains :: (rect: ImRect, v: ImVec2) -> bool {
	return
		rect.Min.x < v.x && v.x < rect.Max.x &&
		rect.Min.y < v.y && v.y < rect.Max.y;
}
// [0..1] -> [0..1]
rescale :: (t: float, min: float, max: float, type: PlotConfig.Scale.Kind) -> float {
    if type == {
    case .Linear; return t;
    case .Log10;  return cast(float) (log10(ImLerp(min, max, t) / min) / log10(max / min));
    }
    return 0;
}

// [0..1] -> [0..1]
rescale_inv :: (t: float, min: float, max: float, type: PlotConfig.Scale.Kind) -> float {
    if type == {
    case .Linear; return t;
    case .Log10;  return cast(float)((pow(max/min, t) * min - min) / (max - min));
    }
    return 0;
}

cursor_to_idx :: (pos: ImVec2, bb: ImRect, conf: PlotConfig, x_min: float, x_max: float) -> int {
    t := ImSaturate((pos.x - bb.Min.x) / (bb.Max.x - bb.Min.x));
    inv := rescale_inv(t, x_min, x_max, conf.scale.type) * (conf.values.count - 1);
    v_idx := cast(int) (inv + 0.5);
    // IM_ASSERT(v_idx >= 0 && v_idx < conf.values.count);

    return v_idx;
}
#scope_export

imgui :: (data: *$T) {
	info := type_info(T);

	ImGui.PushID(xx data);
	defer ImGui.PopID();

	if #complete info.type == {
		case .INTEGER;
		case .FLOAT;
		case .BOOL;
			
		case .STRING;
		case .POINTER;
		case .PROCEDURE;
		case .VOID;
		case .CODE;
		case .STRUCT;
		case .ARRAY;
		case .ANY;
		case .ENUM;
		case .POLYMORPHIC_VARIABLE;
		case .TYPE;
		case .OVERLOAD_SET;
		case .VARIANT;
	}
}

Plot :: (label: string, conf: PlotConfig) -> PlotStatus{
	using ImGui;
	status : PlotStatus = .nothing;

	ys_list := conf.values.ys_list;
	ys_count := conf.values.ys_count;
	colors := conf.values.colors;
	if (conf.values.ys != null) { // draw only a single plot
		ys_list = *conf.values.ys;
		ys_count = 1;
		colors = *conf.values.color;
	}

	style := <<ImGui.GetStyle();
	io := ImGui.GetIO();
	id := ImGui.GetID(label);
	drawlist := ImGui.GetForegroundDrawList();

	frame_bb: ImRect;
	frame_bb.Min = ImGui.GetWindowPos() + ImGui.GetCursorPos();
	frame_bb.Max = frame_bb.Min + conf.frame_size;

	inner_bb: ImRect;
	inner_bb.Min = frame_bb.Min + style.FramePadding;
	inner_bb.Max = frame_bb.Max - style.FramePadding;

	total_bb := frame_bb;
	ImGui.ItemSize(*total_bb, style.FramePadding.y);
	if (!ItemAdd(*total_bb, 0, *frame_bb))
		return status;
	hovered := ItemHoverable(*frame_bb, id);

	RenderFrame(
		frame_bb.Min,
		frame_bb.Max,
		GetColorU32(.FrameBg),
		true,
		style.FrameRounding
	);

	if conf.values.count > 0 {
		res_w: int;
		if conf.skip_small_lines
			res_w = ImMin(cast(int) conf.frame_size.x, conf.values.count);
		else
			res_w = conf.values.count;
		res_w -= 1;
		item_count := conf.values.count - 1;

		x_min := cast(float) conf.values.offset;
		x_max := cast(float) conf.values.offset + conf.values.count - 1.0;
		if conf.values.xs {
			x_min = conf.values.xs[cast(int) x_min];
			x_max = conf.values.xs[cast(int) x_max];
		}

		// Tooltip on hover
		v_hovered := -1;
		if conf.tooltip.show && hovered && Contains(inner_bb, io.MousePos) {
			v_idx := cursor_to_idx(io.MousePos, inner_bb, conf, x_min, x_max);
			data_idx := conf.values.offset + (v_idx % conf.values.count);
			x0 := ifx conf.values.xs then conf.values.xs[data_idx] else 1.0 * v_idx;

			format : string;

			temp_buffer : string;
			for i: 0..ys_count-1 {
				if conf.tooltip.ys_names {
					temp_buffer = tprint(
						conf.tooltip.format,
						conf.tooltip.ys_names[i],
						x0,
						ys_list[i][data_idx]
					);
				}
				else {
					temp_buffer = tprint(conf.tooltip.format, x0, ys_list[i][data_idx]);
				}
				format = tprint("%\n%", format, temp_buffer);
			}
			format = tprint("%\0", format);

			SetTooltip("%s", format.data);
			v_hovered = v_idx;
		}

		t_step := 1.0 / cast(float)res_w;
		inv_scale := 0.0;
		if conf.scale.min != conf.scale.max inv_scale = 1.0 / (conf.scale.max - conf.scale.min);

		if conf.grid_x.show {
			y0 := cast(int)inner_bb.Min.y;
			y1 := cast(int)inner_bb.Max.y;
			if conf.scale.type == {
			case .Linear; {
				cnt := conf.values.count / (conf.grid_x.size / conf.grid_x.subticks);
				inc := 1.0 / cnt;
				for i: 0..cast(int)cnt {
					x0 := cast(int)ImLerp(inner_bb.Min.x, inner_bb.Max.x, i * inc);
					ImDrawList.AddLine(
						drawlist,
						make_ImVec2(cast(float)x0, cast(float)y0),
						make_ImVec2(cast(float)x0, cast(float)y1),
						make_col32(200, 200, 200, xx ifx i % conf.grid_x.subticks then 128 else 255));
				}
			}
			case .Log10; {
				start := 1.0;
				while start < x_max {
					for i: 0..9 {
						x := start * i;
						if x < x_min continue;
						if x > x_max break;
						t := cast(float)(log10(x / x_min) / log10(x_max / x_min));
						x0 := cast(int)ImLerp(inner_bb.Min.x, inner_bb.Max.x, t);
						ImDrawList.AddLine(
							drawlist,
							make_ImVec2(cast(float)x0, cast(float)y0),
							make_ImVec2(cast(float)x0, cast(float)y1),
							make_col32(200, 200, 200, xx ifx i > 1 then 128 else 255));
					}
					start *= 10.0;
				}
			}
			}
		}
		if conf.grid_y.show {
			x0 := cast(int)inner_bb.Min.x;
			x1 := cast(int)inner_bb.Max.x;
			cnt := (conf.scale.max - conf.scale.min) / (conf.grid_y.size / conf.grid_y.subticks);
			inc := 1.0 / cnt;
			for i: 0..cast(int)cnt {
				y0 := cast(int)ImLerp(inner_bb.Min.y, inner_bb.Max.y, i * inc);
				ImDrawList.AddLine(
					drawlist,
					make_ImVec2(cast(float)x0, cast(float)y0),
					make_ImVec2(cast(float)x1, cast(float)y0),
					make_col32(0, 0, 0, xx ifx i % conf.grid_y.subticks then 16 else 64));
			}
		}

		col_hovered := GetColorU32(.PlotLinesHovered);
		col_base := GetColorU32(.PlotLines);


		for i: 0..ys_count-1 {
			if colors {
				if   colors[i] col_base = colors[i];
				else col_base = GetColorU32(.PlotLines);
			}
			v0 := ys_list[i][conf.values.offset];
			t0 := 0.0;
			// Point in the normalized space of our target rectangle
			tp0 := make_ImVec2(t0, 1.0 - ImSaturate((v0 - conf.scale.min) * inv_scale));

			for n: 0..res_w
			{
				t1 := t0 + t_step;
				v1_idx := cast(int)(t0 * item_count + 0.5);
				// IM_ASSERT(v1_idx >= 0 && v1_idx < conf.values.count);
				v1 := ys_list[i][conf.values.offset + (v1_idx + 1) % conf.values.count];
				tp1 := make_ImVec2(
					rescale(t1, x_min, x_max, conf.scale.type),
					1.0 - ImSaturate((v1 - conf.scale.min) * inv_scale)
				);

				// NB: Draw calls are merged together by the DrawList system.
				// Still, we should render our batch are lower level to save a bit of CPU.
				pos0 := ImLerp(inner_bb.Min, inner_bb.Max, tp0);
				pos1 := ImLerp(inner_bb.Min, inner_bb.Max, tp1);

				if v1_idx == v_hovered ImDrawList.AddCircleFilled(drawlist, pos0, 3, col_hovered);

				if (n < res_w) ImDrawList.AddLine(
					drawlist,
					pos0,
					pos1,
					col_base,
					conf.line_thickness
				);

				t0 = t1;
				tp0 = tp1;
			}
		}

		if conf.v_lines.show {
			for i: 0..conf.v_lines.count-1 {
				idx := conf.v_lines.indices[i];
				t1 := rescale(idx * t_step, x_min, x_max, conf.scale.type);
				pos0 := ImLerp(inner_bb.Min, inner_bb.Max, make_ImVec2(t1, 0));
				pos1 := ImLerp(inner_bb.Min, inner_bb.Max, make_ImVec2(t1, 1));
				ImDrawList.AddLine(drawlist, pos0, pos1, make_col32(0xff, 0, 0, 0x88));
			}
		}

		if conf.selection.show {
			if hovered {
				if io.MouseClicked[0] {
					// SetActiveID(id, window);
					// FocusWindow(window);

					v_idx := cursor_to_idx(io.MousePos, inner_bb, conf, x_min, x_max);
					start := conf.values.offset + (v_idx % conf.values.count);
					end := start;
					if conf.selection.sanitize_fn
						end = conf.selection.sanitize_fn(xx (end - start)) + start;
					if end < conf.values.offset + conf.values.count {
						<<conf.selection.start = xx start;
						<<conf.selection.length = xx (end - start);
						status = .selection_updated;
					}
				}
			}

			// if g.ActiveId == id {
				if io.MouseDown[0] {
					v_idx := cursor_to_idx(io.MousePos, inner_bb, conf, x_min, x_max);
					start := <<conf.selection.start;
					end := conf.values.offset + (v_idx % conf.values.count);
					if end > start {
						if conf.selection.sanitize_fn
							end = conf.selection.sanitize_fn(xx (end - start)) + start;
						if end < conf.values.offset + conf.values.count {
							<<conf.selection.length = xx (end - start);
							status = .selection_updated;
						}
					}
				} else {
					// ClearActiveID();
				}
			// }
			fSelectionStep := 1.0 / item_count;
			pos0 := ImLerp(inner_bb.Min, inner_bb.Max,
			               make_ImVec2(fSelectionStep * <<conf.selection.start, 0.0));
			pos1 := ImLerp(
				inner_bb.Min,
				inner_bb.Max,
				make_ImVec2(fSelectionStep * (<<conf.selection.start + <<conf.selection.length), 1.0)
			);
			ImDrawList.AddRectFilled(drawlist, pos0, pos1, make_col32(128, 128, 128, 32));
			ImDrawList.AddRect(drawlist, pos0, pos1, make_col32(128, 128, 128, 128));
		}
	}

	// Text overlay
	// if conf.overlay_text ImDrawList.RenderTextClipped(
	// 	make_ImVec2(frame_bb.Min.x, frame_bb.Min.y + style.FramePadding.y),
	// 	frame_bb.Max, conf.overlay_text,
	// 	null,
	// 	null,
	// 	make_ImVec2(0.5,0.0)
	// );

	return status;
}

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

circle :: (pos: Physics.V2f, r: float) {
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
	set_uniform(simple_shader, "color", Physics.V3f.make(1, 0, 0));
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
