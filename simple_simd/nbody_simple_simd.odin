package simple_simd

import "base:runtime"
import "core:fmt"
import "core:math"
import "core:simd"
import "core:strconv"
import "core:time"
import "base:intrinsics"


PI :: 3.141592653589793
SOLAR_MASS :: 4.0 * PI * PI
DAYS_PER_YEAR :: 365.24
TOTAL_BODIES :: 5

Body :: struct {
    pos_x, pos_y, pos_z: f64,
    vel_x, vel_y, vel_z: f64,
    mass: f64,
}

Bodies :: #soa[]Body

init_system_auto_soa :: #force_inline proc() -> Bodies {
    bodies := make_soa(Bodies, TOTAL_BODIES)
    
    // Sun
    bodies[0] = Body{
        pos_x = 0.0, pos_y = 0.0, pos_z = 0.0,
        vel_x = 0.0, vel_y = 0.0, vel_z = 0.0,
        mass = SOLAR_MASS,
    }
    
    // Jupiter
    bodies[1] = Body{
        pos_x = 4.84143144246472090e+00,
        pos_y = -1.16032004402742839e+00,
        pos_z = -1.03622044471123109e-01,
        vel_x = 1.66007664274403694e-03 * DAYS_PER_YEAR,
        vel_y = 7.69901118419740425e-03 * DAYS_PER_YEAR,
        vel_z = -6.90460016972063023e-05 * DAYS_PER_YEAR,
        mass = 9.54791938424326609e-04 * SOLAR_MASS,
    }
    
    // Saturn
    bodies[2] = Body{
        pos_x = 8.34336671824457987e+00,
        pos_y = 4.12479856412430479e+00,
        pos_z = -4.03523417114321381e-01,
        vel_x = -2.76742510726862411e-03 * DAYS_PER_YEAR,
        vel_y = 4.99852801234917238e-03 * DAYS_PER_YEAR,
        vel_z = 2.30417297573763929e-05 * DAYS_PER_YEAR,
        mass = 2.85885980666130812e-04 * SOLAR_MASS,
    }
    
    // Uranus
    bodies[3] = Body{
        pos_x = 1.28943695621391310e+01,
        pos_y = -1.51111514016986312e+01,
        pos_z = -2.23307578892655734e-01,
        vel_x = 2.96460137564761618e-03 * DAYS_PER_YEAR,
        vel_y = 2.37847173959480950e-03 * DAYS_PER_YEAR,
        vel_z = -2.96589568540237556e-05 * DAYS_PER_YEAR,
        mass = 4.36624404335156298e-05 * SOLAR_MASS,
    }
    
    // Neptune
    bodies[4] = Body{
        pos_x = 1.53796971148509165e+01,
        pos_y = -2.59193146099879641e+01,
        pos_z = 1.79258772950371181e-01,
        vel_x = 2.68067772490389322e-03 * DAYS_PER_YEAR,
        vel_y = 1.62824170038242295e-03 * DAYS_PER_YEAR,
        vel_z = -9.51592254519715870e-05 * DAYS_PER_YEAR,
        mass = 5.15138902046611451e-05 * SOLAR_MASS,
    }
    
    return bodies
}

advance_auto_soa :: #force_inline proc(bodies: ^Bodies, dt: f64, n: int) {
    for _ in 0..<n {
        #unroll for i in 0..<TOTAL_BODIES {
            for j in (i+1)..<TOTAL_BODIES {
                dx := bodies.pos_x[i] - bodies.pos_x[j]
                dy := bodies.pos_y[i] - bodies.pos_y[j]
                dz := bodies.pos_z[i] - bodies.pos_z[j]
                
                distance_square := dx*dx + dy*dy + dz*dz
                distance := math.sqrt_f64(distance_square)
                mag := dt / (distance * distance_square)
                
                fx := dx * mag
                fy := dy * mag
                fz := dz * mag
                
                bodies.vel_x[i] -= fx * bodies.mass[j]
                bodies.vel_y[i] -= fy * bodies.mass[j]
                bodies.vel_z[i] -= fz * bodies.mass[j]
                
                bodies.vel_x[j] += fx * bodies.mass[i]
                bodies.vel_y[j] += fy * bodies.mass[i]
                bodies.vel_z[j] += fz * bodies.mass[i]
            }
        }
        
        #unroll for i in 0..<TOTAL_BODIES {
            bodies.pos_x[i] += bodies.vel_x[i] * dt
            bodies.pos_y[i] += bodies.vel_y[i] * dt
            bodies.pos_z[i] += bodies.vel_z[i] * dt
        }
    }
}

energy_auto_soa :: #force_inline proc(bodies: ^Bodies) -> f64 {
    e := 0.0
    #unroll for i in 0..<TOTAL_BODIES {
        speed2 := bodies.vel_x[i]*bodies.vel_x[i] + bodies.vel_y[i]*bodies.vel_y[i] + bodies.vel_z[i]*bodies.vel_z[i]
        e += 0.5 * bodies.mass[i] * speed2

        for j in (i+1)..<TOTAL_BODIES {
            dx := bodies.pos_x[i] - bodies.pos_x[j]
            dy := bodies.pos_y[i] - bodies.pos_y[j]
            dz := bodies.pos_z[i] - bodies.pos_z[j]
            distance := math.sqrt_f64(dx*dx + dy*dy + dz*dz)
            e -= bodies.mass[i] * bodies.mass[j] / distance
        }
    }
    return e
}

offset_momentum_auto_soa :: #force_inline proc(bodies: ^Bodies) {
    px, py, pz := 0.0, 0.0, 0.0
    #unroll for i in 0..<TOTAL_BODIES {
        px -= bodies.vel_x[i] * bodies.mass[i]
        py -= bodies.vel_y[i] * bodies.mass[i]
        pz -= bodies.vel_z[i] * bodies.mass[i]
    }
    bodies.vel_x[0] = px / SOLAR_MASS
    bodies.vel_y[0] = py / SOLAR_MASS
    bodies.vel_z[0] = pz / SOLAR_MASS
}
