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

    sys := simple_simd.init_system_soa_stack()
    simple_simd.offset_momentum_soa_stack(&sys)

    fmt.printf("Simple simd   energy:  %.9f\n", simple_simd.energy_soa_stack(&sys))

    start:= time.tick_now()
    simple_simd.advance_stack(&sys, 0.01, n)
    duration := time.tick_since(start)
    fmt.eprintf("Simple simd time: %.3f µs\n", f64(duration) / f64(time.Microsecond))

    fmt.printf("Simple simd   energy:  %.9f\n", simple_simd.energy_soa_stack(&sys))
}

run_full_simd :: #force_inline proc() {
    args := runtime.args__
    n := strconv.parse_int(auto_cast args[1]) or_else 1000

    sys := full_simd.init_system_soa_stack()
    full_simd.offset_momentum_soa_stack_simd(&sys)

    fmt.printf("Full simd   energy:  %.9f\n", full_simd.energy_soa_stack_simd(&sys))

    start:= time.tick_now()
    full_simd.advance_stack_simd(&sys, 0.01, n)
    duration := time.tick_since(start)
    fmt.eprintf("Full simd time: %.3f µs\n", f64(duration) / f64(time.Microsecond))

    fmt.printf("Full simd   energy:  %.9f\n", full_simd.energy_soa_stack_simd(&sys))
}

main :: proc() {
    run_simple_simd()
    run_full_simd()
}