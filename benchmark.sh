#!/bin/bash

echo "=== Odin SIMD N-Body Benchmark ==="
echo "Running 100 iterations for each step count..."
echo

run_benchmark() {
    local steps=$1
    local iterations=100
    local simple_total=0
    local full_total=0
    local simple_count=0
    local full_count=0

    echo "Testing with $steps simulation steps:"

    for i in $(seq 1 $iterations); do
        # Clear system cache
        sudo purge > /dev/null 2>&1

        output=$(./odin_simd_playground $steps 2>&1)

        simple_time=$(echo "$output" | grep "Simple simd time:" | sed 's/.*: \([0-9.]*\) µs/\1/')
        if [[ -n "$simple_time" ]]; then
            simple_total=$(echo "$simple_total + $simple_time" | bc -l)
            ((simple_count++))
        fi

        full_time=$(echo "$output" | grep "Full simd time:" | sed 's/.*: \([0-9.]*\) µs/\1/')
        if [[ -n "$full_time" ]]; then
            full_total=$(echo "$full_total + $full_time" | bc -l)
            ((full_count++))
        fi
    done

    if [[ $simple_count -gt 0 ]]; then
        simple_avg=$(echo "scale=3; $simple_total / $simple_count" | bc -l)
        echo "  Simple SIMD: $simple_avg µs (avg of $simple_count runs)"
    fi

    if [[ $full_count -gt 0 ]]; then
        full_avg=$(echo "scale=3; $full_total / $full_count" | bc -l)
        echo "  Full SIMD:   $full_avg µs (avg of $full_count runs)"
    fi

    if [[ $simple_count -gt 0 && $full_count -gt 0 ]]; then
        speedup=$(echo "scale=2; $simple_total / $full_total" | bc -l)
        echo "  Speedup:     ${speedup}x"
    fi

    echo
}

echo "Building program..."
odin build . -o:aggressive -microarch:native -no-bounds-check -out:odin_simd_playground
echo "Build complete."
echo

step_counts=(1000 10000 1000000 10000000)

for steps in "${step_counts[@]}"; do
    run_benchmark $steps
done

echo "Benchmark complete!"