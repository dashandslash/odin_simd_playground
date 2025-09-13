package main

import "core:fmt"

Test :: struct {
    x, y: f64
}

TestSOA :: #soa[]Test

main :: proc() {
    soa := make_soa(TestSOA, 3)
    defer delete_soa(soa)

    soa[0].x = 1.0
    soa[0].y = 2.0
    soa[1].x = 3.0
    soa[1].y = 4.0
    soa[2].x = 5.0
    soa[2].y = 6.0

    // Print memory layout
    fmt.printf("SOA base address: %p\n", &soa)
    fmt.printf("Element 0 address: %p\n", &soa[0])
    fmt.printf("Element 0.x address: %p\n", &soa[0].x)
    fmt.printf("Element 0.y address: %p\n", &soa[0].y)
    fmt.printf("Element 1.x address: %p\n", &soa[1].x)
    fmt.printf("Element 1.y address: %p\n", &soa[1].y)
    fmt.printf("Element 2.x address: %p\n", &soa[2].x)
    fmt.printf("Element 2.y address: %p\n", &soa[2].y)

    // Print values
    fmt.printf("Values: [%.1f, %.1f], [%.1f, %.1f], [%.1f, %.1f]\n",
               soa[0].x, soa[0].y, soa[1].x, soa[1].y, soa[2].x, soa[2].y)
}