#import "Compiler";
#import "Basic";

build_usb :: () -> Workspace {
	w := compiler_create_workspace("USB_Driver");

	import : [..] string;
	opts := get_build_options(w);
	array_add(*import, ..opts.import_path);
	array_add(*import, "./src/");
	array_add(*import, "./src/Driver/");
	opts.output_executable_name = "USB";
	opts.import_path = import;
	opts.text_output_flags = 0;
	opts.output_type = .DYNAMIC_LIBRARY;
	set_build_options(opts, w);

	return w;
}

build_client :: () -> Workspace {
	w := compiler_create_workspace("Client");
	import : [..] string;
	opts := get_build_options(w);
	array_add(*import, ..opts.import_path);
	array_add(*import, "./src/");
	array_add(*import, "./src/Client/");
	opts.output_executable_name = "Client";
	opts.import_path = import;
	opts.text_output_flags = 0;
	set_build_options(opts, w);
	return w;
}

build :: () {
	usb    := build_usb();
	client := build_client();

	compiler_begin_intercept(client);
	compiler_begin_intercept(usb);

	add_build_file("src/Client/main.jai", client);
	add_build_file("src/Driver/main.jai", usb);
	
	while true {
		message := compiler_wait_for_message();

		if message.kind == {
			case .COMPLETE; break;
		}
	}

	compiler_end_intercept(client);
	compiler_end_intercept(usb);

	
	set_build_options_dc(.{ do_output=false });
}


#run build();
