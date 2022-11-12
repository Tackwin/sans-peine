
#scope_file

Data :: #import "Data";
Interface :: #import "Interface";
ImGui :: #import "my_imgui";
#import "Basic";

#scope_export

#program_export "input_timeline"
input_timeline :: (systems: *[..] Data.System_Reading) {
	c: Context;
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

		item_size :: ImGui.ImVec2.{10, 20};
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
				format :: #string DONE
System %1
time: %2
mag[0] = t: %3 [%4 - %5]
mag[1] = t: %6 [%7 - %8]
DONE;
				formatted := tprint(
					format,
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
			}

			beg_x := pos.x + item_size.x - scroll_x;
			end_x := pos.x - scroll_x;

			color := item_col_invalid;
			if system.valid color = item_col_valid;
			else for system.magnetometer_readings if it.valid {
				color = item_col_partial;
				break;
			}

			if 0 < beg_x && end_x < inner_bb.Max.x - inner_bb.Min.x
			{
				p := GetItemRectMin() - inner_bb.Min;
				s := GetItemRectSize();
				rect(p, v2(s.x, s.y * 0.8), color);

				for system.magnetometer_readings {
					pp := p + v2(it_index * s.x / 2, s.y * 0.8);
					ss := v2(s.x / 2, 0.2 * s.y);

					color := item_col_interpolated;
					if it.measured color = item_col_measured;
					rect(pp, ss, color);
				}
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

auto_scroll := true;