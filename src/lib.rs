//! Euler-Lagrange equation solver.
//!
//! Solves variational problems by computing the Euler-Lagrange equation:
//! d/dt(∂L/∂q̇) - ∂L/∂q = 0

/// A point in configuration space: position and velocity.
#[derive(Debug, Clone, Copy)]
pub struct State {
    pub q: f64,
    pub qdot: f64,
}

/// Represents a Lagrangian L(q, qdot, t) via closures.
pub struct Lagrangian<F, G>
where
    F: Fn(f64, f64, f64) -> f64,
    G: Fn(f64, f64, f64) -> f64,
{
    /// ∂L/∂q̇ (partial derivative w.r.t. velocity)
    pub dL_dqdot: F,
    /// ∂L/∂q (partial derivative w.r.t. position)
    pub dL_dq: G,
}

impl<F, G> Lagrangian<F, G>
where
    F: Fn(f64, f64, f64) -> f64,
    G: Fn(f64, f64, f64) -> f64,
{
    /// Compute the Euler-Lagrange residual: d/dt(∂L/∂q̇) - ∂L/∂q
    /// using finite differences for the time derivative.
    pub fn residual(&self, q: f64, qdot: f64, t: f64, dt: f64) -> f64 {
        let dL_dqdot_now = (self.dL_dqdot)(q, qdot, t);
        let dL_dqdot_next = (self.dL_dqdot)(q, qdot, t + dt);
        let d_dt_dL_dqdot = (dL_dqdot_next - dL_dqdot_now) / dt;
        let dL_dq = (self.dL_dq)(q, qdot, t);
        d_dt_dL_dqdot - dL_dq
    }
}

/// Integrate equations of motion from the Euler-Lagrange equation using symplectic Euler.
/// Returns a vector of (t, State) pairs.
pub fn integrate<F, G>(
    lagrangian: &Lagrangian<F, G>,
    initial: State,
    t_start: f64,
    t_end: f64,
    dt: f64,
) -> Vec<(f64, State)>
where
    F: Fn(f64, f64, f64) -> f64,
    G: Fn(f64, f64, f64) -> f64,
{
    let mut results = Vec::new();
    let mut state = initial;
    let mut t = t_start;
    let eps = 1e-8;

    results.push((t, state));

    while t < t_end - eps {
        // Approximate acceleration from EL equation residual = 0
        // Use Newton-like step: adjust qdot to drive residual to zero
        let r = lagrangian.residual(state.q, state.qdot, t, dt);
        let r_plus = lagrangian.residual(state.q, state.qdot + eps, t, dt);
        let dr_dqdot = (r_plus - r) / eps;

        let qdotdot = if dr_dqdot.abs() > 1e-14 {
            -r / dr_dqdot / dt
        } else {
            0.0
        };

        // Symplectic Euler: update velocity first, then position
        state.qdot += qdotdot * dt;
        state.q += state.qdot * dt;
        t += dt;

        results.push((t, state));
    }

    results
}

/// Convenience: create a Lagrangian for simple harmonic motion L = ½mv² - ½kq²
pub fn harmonic_oscillator(mass: f64, spring_k: f64) -> Lagrangian<impl Fn(f64, f64, f64) -> f64, impl Fn(f64, f64, f64) -> f64> {
    Lagrangian {
        dL_dqdot: move |_, qdot, _| mass * qdot,
        dL_dq: move |q, _, _| -spring_k * q,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_harmonic_oscillator_conservation() {
        let lang = harmonic_oscillator(1.0, 1.0);
        let initial = State { q: 1.0, qdot: 0.0 };
        let traj = integrate(&lang, initial, 0.0, 10.0, 0.001);
        // Energy: E = ½qdot² + ½q² should stay ~1.0
        for (_, s) in traj.iter().step_by(100) {
            let e = 0.5 * s.qdot * s.qdot + 0.5 * s.q * s.q;
            assert!((e - 1.0).abs() < 0.1, "Energy drift too large: {e}");
        }
    }
}
