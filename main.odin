package main

import "base:runtime"
import "core:fmt"
import "core:strconv"
import "core:time"

import "full_simd"
import "simple_simd"


run_simple_simd :: #force_inline proc() {
    args := runtime.args__
    n := strconv.parse_int(auto_cast args[1]) or_else 1000

    bodies := simple_simd.init_system_auto_soa()
    defer delete_soa(bodies)
    simple_simd.offset_momentum_auto_soa(&bodies)

    fmt.printf("Simple simd   energy:  %.9f\n", simple_simd.energy_auto_soa(&bodies))

    start:= time.tick_now()
    simple_simd.advance_auto_soa(&bodies, 0.01, n)
    duration := time.tick_since(start)
    fmt.eprintf("Simple simd time: %.3f µs\n", f64(duration) / f64(time.Microsecond))

    fmt.printf("Simple simd   energy:  %.9f\n", simple_simd.energy_auto_soa(&bodies))
}

run_full_simd :: #force_inline proc() {
    args := runtime.args__
    n := strconv.parse_int(auto_cast args[1]) or_else 1000

    bodies := full_simd.init_system_soa_stack()
    full_simd.offset_momentum_soa_stack_simd(&bodies)

    fmt.printf("Full simd   energy:  %.9f\n", full_simd.energy_soa_stack_simd(&bodies))

    start:= time.tick_now()
    full_simd.advance_stack_simd(&bodies, 0.01, n)
    duration := time.tick_since(start)
    fmt.eprintf("Full simd time: %.3f µs\n", f64(duration) / f64(time.Microsecond))

    fmt.printf("Full simd   energy:  %.9f\n", full_simd.energy_soa_stack_simd(&bodies))
}

main :: proc() {
    run_simple_simd()
    run_full_simd()
}