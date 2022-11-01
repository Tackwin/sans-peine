#scope_file

#import "Basic";


#scope_export
append :: (timer: *Timer($Is_Series), x: float64) {
	#if Is_Series array_add(*timer.values, x);
	else          timer.value = x;
}

scope_timer :: (name: $string) #expand {
	t0 := to_nanoseconds(current_time_monotonic());
	`defer {
		t1 := to_nanoseconds(current_time_monotonic());

		#insert #run tprint("append(*timer_%, (t1 - t0) / 1_000_000_000.0);\n", name);
	}
}

Timer :: struct(Is_Series: bool) {
	#if Is_Series {
		values: [..] float64;
	}
	else {
		value: float64;
	}
}

timer_update: Timer(false);
timer_RBF: Timer(false);
timer_solve: Timer(false);
timer_loop: Timer(false);
timer_input_replay: Timer(true);