
#scope_file

Data :: #import "Data";
Interface :: #import "Interface";
ImGui :: #import "my_imgui";
#import "Math";
#import "Basic";

#scope_export

#program_export "input_timeline"
input_timeline :: (systems: *[..] Data.System_Reading) {
	c: type_of(context);
	push_context c {
		using ImGui;
		ImGui.Text("Input timeline");

		style := <<ImGui.GetStyle();
		frame_bb: ImRect;
		frame_bb.Min = ImGui.GetWindowPos();
		frame_bb.Max = frame_bb.Min + ImGui.GetWindowSize() + ImGui.GetWindowPos() - frame_bb.Min;

		v2 :: ImGui.make_ImVec2;
		inner_bb: ImRect;
		inner_bb.Min = frame_bb.Min + style.FramePadding;
		inner_bb.Max = frame_bb.Max - v2(3, 3) * style.FramePadding;

		pos := ImGui.GetCursorPos();

		rect :: (p: ImGui.ImVec2, s: ImGui.ImVec2, c: u32) #expand {
			ImGui.ImDrawList.AddRectFilled(drawlist, inner_bb.Min + p, inner_bb.Min + p + s, c);
		}

		item_size :: ImGui.ImVec2.{10, 30};
		item_col_invalid  : u32 : 0xFFFF0000;
		item_col_\ valid  : u32 : 0xFF00FF00;
		item_col_partial  : u32 : 0xFFFFFF00;

		item_col_measured     : u32 : 0xFF00FF00;
		item_col_interpolated : u32 : 0xFFFFFF00;
		item_padding :: 2;

		beg := 0;
		end := systems.count;
		// if end > 100 {
		// 	end = systems.count;
		// 	beg = end - 100;
		// }

		ImGui.Checkbox("Autoscroll", *auto_scroll);

		ImGui.BeginChild(
			"Systems",
			v2((inner_bb.Max.x - inner_bb.Min.x), 0),
			flags=.HorizontalScrollbar
		);
		drawlist := ImGui.GetWindowDrawList();
		ImGui.ImDrawList.PushClipRect(drawlist, inner_bb.Min, inner_bb.Max, true);
		if auto_scroll ImGui.SetScrollX(ImGui.GetScrollMaxX());
		scroll_x := ImGui.GetScrollX();


		handle_system :: (system: Data.System_Reading, idx: int) #expand {
			ImGui.SetCursorPos(pos);
			ImGui.Dummy(item_size);
			if ImGui.IsItemHovered() {
				ImGui.BeginTooltip();
				defer ImGui.EndTooltip();
				format_mag :: #string DONE
System %1
time: %2
mag[0] = t: %3 [%4 - %5]
mag[1] = t: %6 [%7 - %8]
DONE;
				formatted := tprint(
					format_mag,
					idx,
					system.time,
					system.magnetometer_readings[0].valid,
					system.magnetometer_readings[0].interpolated_beg,
					system.magnetometer_readings[0].interpolated_end,
					system.magnetometer_readings[1].valid,
					system.magnetometer_readings[1].interpolated_beg,
					system.magnetometer_readings[1].interpolated_end
				);
				ImGui.Text(tprint("%\0", formatted).data);

				format_acc :: #string DONE
acc[0] = t: % [% - %]
acc[1] = t: % [% - %]
DONE;
				formatted = tprint(
					format_acc,
					system.accelerometer_readings[0].valid,
					system.accelerometer_readings[0].interpolated_beg,
					system.accelerometer_readings[0].interpolated_end,
					system.accelerometer_readings[1].valid,
					system.accelerometer_readings[1].interpolated_beg,
					system.accelerometer_readings[1].interpolated_end
				);
				ImGui.Text(tprint("%\0", formatted).data);
			}

			beg_x := pos.x + item_size.x - scroll_x;
			end_x := pos.x - scroll_x;

			color := item_col_invalid;
			if system.valid color = item_col_valid;
			else {
				for system.magnetometer_readings if it.valid {
					color = item_col_partial;
					break;
				}
				for system.accelerometer_readings if it.valid {
					color = item_col_partial;
					break;
				}
			}

			if 0 < beg_x && end_x < inner_bb.Max.x - inner_bb.Min.x
			{
				p := GetItemRectMin() - inner_bb.Min;
				s := GetItemRectSize();
				rect(p, v2(s.x, s.y * 0.7), color);

				for system.magnetometer_readings {
					pp := p + v2(it_index * s.x / 2, s.y * 0.7);
					ss := v2(s.x / 2, 0.125 * s.y);

					color := item_col_interpolated;
					if it.measured color = item_col_measured;
					rect(pp, ss, color);
				}

				rect(p + v2(0, s.y * (0.65)), v2(s.x, 0.05 * s.y), 0xFF000000);
				rect(p + v2(0, s.y * (0.7 + 0.125)), v2(s.x, 0.05 * s.y), 0xFF000000);
				
				for system.accelerometer_readings {
					pp := p + v2(it_index * s.x / 2, s.y * (0.7 + 0.125 + 0.05));
					ss := v2(s.x / 2, 0.125 * s.y);

					color := item_col_interpolated;
					if it.measured color = item_col_measured;
					rect(pp, ss, color);
				}

				rect(p + v2(0.45 * s.x, s.y * 0.7), v2(0.1 * s.x, 0.3 * s.y), 0xFF000000);
			}
			pos.x += item_size.x + item_padding;

		}

		for idx: beg..end-1 {
			system := systems.data[idx];
			handle_system(system, idx);
		}
		ImGui.ImDrawList.PopClipRect(drawlist);
		ImGui.EndChild();
	}
}


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

#program_export "plot"
plot :: (label: string, conf: PlotConfig) -> PlotStatus {
	c: type_of(context);
	push_context c {

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
		drawlist := ImGui.GetWindowDrawList();
		scroll_y := ImGui.GetScrollY();

		frame_bb: ImRect;
		frame_bb.Min = ImGui.GetWindowPos() + ImGui.GetCursorPos();
		frame_bb.Min.y -= scroll_y;
		frame_bb.Max = frame_bb.Min + conf.frame_size;

		inner_bb: ImRect;
		inner_bb.Min = frame_bb.Min + style.FramePadding;
		inner_bb.Max = frame_bb.Max - style.FramePadding;

		total_bb := frame_bb;
		ImGui.ItemSize(*total_bb, style.FramePadding.y);
		if (!ItemAdd(*total_bb, 0, *frame_bb)) {
			return status;
		}
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
}

auto_scroll := true;