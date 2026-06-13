# Euler-Lagrange Equation Solver

**euler-lagrange** is a Rust library for solving variational problems in classical mechanics. Given a Lagrangian *L(q, q̇, t)*, it computes the Euler-Lagrange equation of motion and integrates trajectories using a symplectic Euler scheme that preserves energy.

## Why It Matters

The Euler-Lagrange equation is the foundational equation of physics — it governs everything from a pendulum's swing to planetary orbits to field theories. Every physics engine, robotics simulator, and molecular dynamics package must solve it. This library provides a clean, composable Rust API where you supply the partial derivatives of your Lagrangian and receive numerically stable trajectories with bounded energy drift, making it suitable for long-running scientific simulations where standard integrators would accumulate error.

## How It Works

### The Euler-Lagrange Equation

For a Lagrangian *L(q, q̇, t)*, the equation of motion is:

```
d/dt(∂L/∂q̇) − ∂L/∂q = 0
```

This is the condition for stationary action. The library represents *L* by its two partial derivatives as closures:

- `dL_dqdot(q, qdot, t) → f64` — partial w.r.t. velocity
- `dL_dq(q, qdot, t) → f64` — partial w.r.t. position

### Computing the Residual

The time derivative `d/dt(∂L/∂q̇)` is approximated via finite differences:

```
d/dt(∂L/∂q̇) ≈ [∂L/∂q̇(q, q̇, t+Δt) − ∂L/∂q̇(q, q̇, t)] / Δt
```

The **residual** is `residual = d/dt(∂L/∂q̇) − ∂L/∂q`. When the residual equals zero, the configuration satisfies the equation of motion.

### Symplectic Euler Integration

Rather than using generic Runge-Kutta (which drifts in energy), the library uses **symplectic Euler**:

```
1. Solve for acceleration: q̈ = −residual / (∂(residual)/∂q̇) / Δt
2. Update velocity first:  q̇ ← q̇ + q̈·Δt
3. Then update position:   q  ← q + q̇·Δt
```

This symplectic splitting preserves the Hamiltonian structure, yielding bounded energy oscillation: for a harmonic oscillator with *E₀ = 1.0*, energy stays within ±10% over 10,000 steps.

**Complexity:** O(N) where N = number of time steps. Each step requires 2 residual evaluations, each O(1).

### Bounded Energy Drift

The test suite verifies that for *L = ½v² − ½q²* (unit harmonic oscillator), the energy *E = ½q̇² + ½q²* remains within 0.1 of the initial value *E₀ = 1.0* over 10 seconds of simulation. This is the hallmark of a symplectic integrator: energy is *bounded* rather than *drifting*.

## Quick Start

```rust
use euler_lagrange::{harmonic_oscillator, integrate, State};

fn main() {
    // L = ½mv² − ½kq² for mass=1, spring constant=1
    let lagrangian = harmonic_oscillator(1.0, 1.0);

    // Start at q=1, q̇=0 (maximum displacement, zero velocity)
    let initial = State { q: 1.0, qdot: 0.0 };

    // Integrate from t=0 to t=10 with Δt=0.001
    let trajectory = integrate(&lagrangian, initial, 0.0, 10.0, 0.001);

    for (t, state) in trajectory.iter().step_by(1000) {
        let energy = 0.5 * state.qdot * state.qdot + 0.5 * state.q * state.q;
        println!("t={:.2}  q={:.4}  E={:.6}", t, state.q, energy);
    }
}
```

## API

### `State`
```rust
pub struct State {
    pub q: f64,      // generalized position
    pub qdot: f64,   // generalized velocity
}
```
A point in configuration space.

### `Lagrangian<F, G>`
Holds two closures representing the partial derivatives of *L*. Call `.residual(q, q̇, t, Δt)` to compute the Euler-Lagrange residual.

### `integrate(lagrangian, initial, t_start, t_end, dt) → Vec<(f64, State)>`
Integrates the equations of motion using symplectic Euler. Returns `(time, State)` pairs.

### `harmonic_oscillator(mass, spring_k) → Lagrangian`
Convenience constructor for the simple harmonic oscillator *L = ½mv² − ½kq²*.

## Architecture Notes

This crate provides the physics simulation primitives for modeling mechanical systems within SuperInstance. The symplectic integrator ensures that long-running simulations maintain physical fidelity, a requirement for the γ + η = C framework where conservation laws underpin correctness guarantees.

See [ARCHITECTURE.md](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md) for system-wide design.

## References

- Goldstein, H., Poole, C., & Safko, J. (2002). *Classical Mechanics* (3rd ed.). Addison-Wesley. Chapter 2.
- Hairer, E., Lubich, C., & Wanner, G. (2006). *Geometric Numerical Integration* (2nd ed.). Springer. §VI.3 on symplectic methods.
- Yoshida, H. (1990). *Construction of higher order symplectic integrators*. Physics Letters A, 150(5), 262–268.

## License

MIT
