#scope_export

___tracy_source_location_data :: struct {
    name     : *c.char;
    function : *c.char;
    file     : *c.char;
    line     : u32;
    color    : u32;
}

___tracy_c_zone_context :: struct {
    id     : u32;
    active : c.int;
}

TracyCZoneCtx :: ___tracy_c_zone_context;

___tracy_set_thread_name                   :: (name: *c.char)                                                                                                                   #foreign tracy;
___tracy_alloc_srcloc                      :: (line: u32, source: *c.char, sourceSz: c.size_t, function: *c.char, functionSz: c.size_t) -> u64                                  #foreign tracy;
___tracy_alloc_srcloc_name                 :: (line: u32, source: *c.char, sourceSz: c.size_t, function: *c.char, functionSz: c.size_t, name: *c.char, nameSz: c.size_t) -> u64 #foreign tracy;
___tracy_emit_zone_begin                   :: (srcloc: *___tracy_source_location_data, active: c.int) -> TracyCZoneCtx                                                          #foreign tracy;
___tracy_emit_zone_begin_callstack         :: (srcloc: *___tracy_source_location_data, depth: c.int, active: c.int) -> TracyCZoneCtx                                            #foreign tracy;
___tracy_emit_zone_begin_alloc             :: (srcloc: u64, active: c.int) -> TracyCZoneCtx                                                                                     #foreign tracy;
___tracy_emit_zone_begin_alloc_callstack   :: (srcloc: u64, depth: c.int, active: c.int) -> TracyCZoneCtx                                                                       #foreign tracy;
___tracy_emit_zone_end                     :: (ctx: TracyCZoneCtx)                                                                                                              #foreign tracy;
___tracy_emit_zone_text                    :: (ctx: TracyCZoneCtx, txt: *c.char, size: c.size_t)                                                                                #foreign tracy;
___tracy_emit_zone_name                    :: (ctx: TracyCZoneCtx, txt: *c.char, size: c.size_t)                                                                                #foreign tracy;
___tracy_emit_zone_color                   :: (ctx: TracyCZoneCtx, color: u32)                                                                                                  #foreign tracy;
___tracy_emit_zone_value                   :: (ctx: TracyCZoneCtx, value: u64)                                                                                                  #foreign tracy;
___tracy_emit_memory_alloc                 :: (ptr: *void, size: c.size_t, secure: c.int)                                                                                       #foreign tracy;
___tracy_emit_memory_alloc_callstack       :: (ptr: *void, size: c.size_t, depth: c.int, secure: c.int)                                                                         #foreign tracy;
___tracy_emit_memory_free                  :: (ptr: *void, secure: c.int)                                                                                                       #foreign tracy;
___tracy_emit_memory_free_callstack        :: (ptr: *void, depth: c.int, secure: c.int)                                                                                         #foreign tracy;
___tracy_emit_memory_alloc_named           :: (ptr: *void, size: c.size_t, secure: c.int, name: *c.char)                                                                        #foreign tracy;
___tracy_emit_memory_alloc_callstack_named :: (ptr: *void, size: c.size_t, depth: c.int, secure: c.int, name: *c.char)                                                          #foreign tracy;
___tracy_emit_memory_free_named            :: (ptr: *void, secure: c.int, name: *c.char)                                                                                        #foreign tracy;
___tracy_emit_memory_free_callstack_named  :: (ptr: *void, depth: c.int, secure: c.int, name: *c.char)                                                                          #foreign tracy;
___tracy_emit_message                      :: (txt: *c.char, size: c.size_t, callstack: c.int)                                                                                  #foreign tracy;
___tracy_emit_messageL                     :: (txt: *c.char, callstack: c.int)                                                                                                  #foreign tracy;
___tracy_emit_messageC                     :: (txt: *c.char, size: c.size_t, color: u32, callstack: c.int)                                                                      #foreign tracy;
___tracy_emit_messageLC                    :: (txt: *c.char, color: u32, callstack: c.int)                                                                                      #foreign tracy;
___tracy_emit_frame_mark                   :: (name: *c.char)                                                                                                                   #foreign tracy;
___tracy_emit_frame_mark_start             :: (name: *c.char)                                                                                                                   #foreign tracy;
___tracy_emit_frame_mark_end               :: (name: *c.char)                                                                                                                   #foreign tracy;
___tracy_emit_frame_image                  :: (image: *void, w: u16, h: u16, offset: u8, flip: c.int)                                                                           #foreign tracy;
___tracy_emit_plot                         :: (name: *c.char, val: c.double)                                                                                                    #foreign tracy;
___tracy_emit_message_appinfo              :: (txt: *c.char, size: c.size_t)                                                                                                    #foreign tracy;

#if OS == .WINDOWS {
    tracy :: #library,no_dll "../lib/TracyClient";
} else #if OS == .LINUX {
    // tracy :: #foreign_library "linux/tracy";
} else #if OS == .MACOS {
    // tracy :: #foreign_library "mac/libtracy";
}

#load "c.jai";
#load "wrapper.jai";