// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

// Procedures for configurational bias Monte-Carlo (CBMC)
//
// This implementation - up to now - deals only with non-branched molecules (i.e. united atom).
//
// The probability of creating a trial position {b} is:
// p({b}) d{b} = p(l, theta, phi) l^2 sin(theta) d(l) d(theta) d(phi)
//             = l^2 p(l) dl
//               sin(theta) p(theta) d(theta)
//               p(phi) d(phi)
//
// where
// - {b}: trial position (vector)
// - l: bond length
// - theta: bond angle
// - phi: dihedral angle
// - p(.): probability density function (non normalized)
//
// CBMC attempts to create trial positions by sampling values for
// l, theta and phi yielding high acceptances with respect to the
// Boltzmann factor of the energy of the trial positions.
//
// # Note:
//
// For the initial implementation we will assume that bond stretching and
// angle bending potentials are harmonic functions.
// The resulting Boltzmann factor is _similar_ to a Gaussian distribution.
// We can use this similarity to improve the sampling - using rejection sampling -
// of the bond length as well as the bond angle.
// Naive sampling (random sampling on a sphere) is slower by a factor ~ 2-10x.
// See the work of Vlugt et al. "Improving the efficiency of the configurational-bias
// Monte Carlo algorithm", Molecular Physics, 1998, Vol. 94, No. 4, 727-733
//
// In the future, we might implement a different sampling scheme for arbitrary potentials.

use std::f64::consts::PI;
use rand::distributions::{Sample, Normal, Range};
use rand::Rng;

use sys::{System, Molecule};
use energy::{AnglePotential, DihedralPotential, BondPotential, Harmonic, Potential};
use types::{Vector3D, Matrix3};

/// Returns a position given an array of positions and weights.
///
/// Computes the CDF (cumulative density function) of the weights,
/// then uses binary search to find the position.
// pub fn select_position<T: Rng>(
//     trialpos: &[Vector3D],
//     trial_weights: &[f64],
//     rng: &mut Box<T>
// ) -> Vector3D {
//     let length = trial_position.len()
//     let range = Range::new(0.0, 1.0);

//     // compute the CDF
//     let cumsum: Vec<f64> = trial_weights.iter()
//         .scan(0.0, |sum, &x| {
//             *sum += x;
//             Some(*sum)
//         })
//         .collect();

//     // elements must not be NaN else we have a problem...
//     // see this discussion:
//     // http://stackoverflow.com/questions/28247990/how-to-do-a-binary-search-on-a-vec-of-floats?noredirect=1&lq=1
//     let seek = range.sample(rng);
//     // unwrap will fail with NaN
//     let selection = cum_sum.binary_search_by(|weight| weight.partial_cmp(&seek).unwrap());
//     let select_id = match selection { Err(x) | Ok(x) => x };
//     trial_position[select_id]
// }

// Note: linear search should be sufficient for small cb_steps
/// Returns a position given an array of tiral positions and weights.
///
/// Uses linear search to make selection.
///
/// # Note
///
/// Does not check if weights are normalized, i.e
/// 0 <= w_i <= 1.
pub fn select_position<T: Rng>(
    trialpos: &[Vector3D],
    trial_weights: &[f64],
    rng: &mut T
) -> Vector3D {

    let mut range = Range::new(0.0, 1.0);
    let mut R = range.sample(rng);
    for (i, &weight) in trial_weights.iter().enumerate() {
        if R <= weight { return trialpos[i] };
        R -= weight;
    }
    *trialpos.last().unwrap()
}

/// Creates the position of a segment of the molecule.
pub fn trial_position<T: Rng>(
    system: &System,
    molecule: &Molecule,
    trialpos: &[Vector3D],
    beta: f64,
    rng: &mut T
) -> Vector3D {

    // the number of already grown atoms == length
    let l = trialpos.len();
    // index of the particle to grow
    let trial_pid = molecule.start() + l;
    match l {
        // for now, don't regrow the first particle
        0 => system[trial_pid].position.clone(),
        1 => {
            // get bond potential between the 0'th and the 1'th segment
            let bond_potential = system.bond_potentials(trial_pid-1, trial_pid);
            // check if the potential is a harmonic potential
            let length = match bond_potential[0].as_any().downcast_ref::<Harmonic>() {
                Some(harmonic) => create_bond_length_harmonic(beta, harmonic, rng),
                None => fatal_error!(
                    "CBMC is currently only implemented for harmonic bond and angle potentials.")
            };
            position_from_length(&trialpos[0], length, rng)
        },
        2 => {
            // get bond potential between 1st and 2nd segment
            // get angle potential between 0th, 1th and 2nd segment
            // check if the potential is a harmonic potential
            let bond_potential = system.bond_potentials(trial_pid-1, trial_pid);
            let length = match bond_potential[0].as_any().downcast_ref::<Harmonic>() {
                Some(harmonic) => create_bond_length_harmonic(beta, harmonic, rng),
                None => fatal_error!(
                    "CBMC is currently only implemented for harmonic bond and angle potentials.")
            };

            let angle_potential = system.angle_potentials(trial_pid-2, trial_pid-1, trial_pid);
            let theta = match angle_potential[0].as_any().downcast_ref::<Harmonic>() {
                Some(harmonic) => create_bond_angle_harmonic(beta, harmonic, rng),
                None => fatal_error!(
                    "CBMC is currently only implemented for harmonic bond and angle potentials.")
            };
            position_from_angle(&trialpos, length, theta, rng)
        },
        _ => {
            // all contributions
            let bond_potential = system.bond_potentials(trial_pid-1, trial_pid);
            let length = match bond_potential[0].as_any().downcast_ref::<Harmonic>() {
                Some(harmonic) => create_bond_length_harmonic(beta, harmonic, rng),
                None => fatal_error!(
                    "CBMC is currently only implemented for harmonic bond and angle potentials.")
            };

            let angle_potential = system.angle_potentials(trial_pid-2, trial_pid-1, trial_pid);
            let dihedral_potential = system.dihedral_potentials(trial_pid-3, trial_pid-2, trial_pid-1, trial_pid);
            let (theta, phi) = match angle_potential[0].as_any().downcast_ref::<Harmonic>() {
                Some(harmonic) => create_bond_and_dihedral_angle(beta, harmonic, &dihedral_potential[0], rng),
                None => fatal_error!(
                    "CBMC is currently only implemented for harmonic bond and angle potentials.")
            };
            // use only the last three positions
            let (_, last_three_positions) = trialpos.split_at(l-3);
            position_from_angle_and_dihedral(last_three_positions, length, theta, phi, rng)
        }
    }
}

/// Creates a bond length given a harmonic bond potential.
///
/// The probability of the length is proportional to
/// the Boltzmann distribution of a harmonic bond potential
/// times the squared length:
///
/// p(l) ~ l^2 * exp[-beta * u_bond(l)] d(l)
fn create_bond_length_harmonic<T: Rng>(
    beta: f64,
    harmonic_potential: &Harmonic,
    rng: &mut T
) -> f64 {
    // Why is dynamic dispatch not working? I.e. putting rng: &mut Box<Rng> as fn arg??

    // standard deviation of the gaussian distribution
    let sigma = (1f64 / (beta * harmonic_potential.k)).sqrt();
    // compute factor M for rejection sampling so that
    // p(l) <= normal(l) * M
    // taken from Frenkel, Smit who took it from Numerical Recipes.
    let M = (harmonic_potential.x0 + 3.0 * sigma).powi(2);
    // create bond length according to normal distribution
    let mut normal = Normal::new(harmonic_potential.x0, sigma);

    // use rejection sampling to sample from non-gaussian probability density
    let mut length = 0f64;
    loop {
        length = normal.sample(rng);
        // accept/reject
        if rng.next_f64() <= length * length / M { break }
    }
    length // return length
}

/// Creates an angle given a harmonic bond angle potential.
///
/// p(theta) ~ sin(theta) * exp[-beta * u_angle(theta)] d(theta)
///
/// The probability for the angle is proportional to
/// the Boltzmann distribution of the angle potential times the sine.
fn create_bond_angle_harmonic<T: Rng>(
    beta: f64,
    harmonic_potential: &Harmonic,
    rng: &mut T
) -> f64 {
    // standard deviation of the gaussian distribution
    let sigma = (1f64 / (beta * harmonic_potential.k)).sqrt();
    // create bond length according to normal distribution
    let mut normal = Normal::new(harmonic_potential.x0, sigma);

    // use rejection sampling to sample from non-gaussian probability density
    let mut theta = 0f64;
    loop {
        theta = normal.sample(rng);
        // accept/reject
        if rng.next_f64() <= theta.sin() { break }
    }
    theta // return angle
}

/// Creates a bond angle `theta` and a dihedral angle `phi`.
fn create_bond_and_dihedral_angle<T: Rng>(
    beta: f64,
    harmonic_potential: &Harmonic,
    dihedral_potential: &Box<DihedralPotential>,
    rng: &mut T
) -> (f64, f64) {

    // Sampling range for dihedral angle
    let mut phi_range = Range::new(0.0, 2.0 * PI);
    let mut phi = 0f64; // dihedral angle
    let mut theta = 0f64; // bonding angle
    loop {
        // Create angle according to Boltzmann factor of harmonic potential
        theta = create_bond_angle_harmonic(beta, harmonic_potential, rng);
        let angle_energy = harmonic_potential.energy(theta);
        // Create point on a cone obeying angle theta,
        // that is - an angle from an uniform distribution about the connection axis.
        phi = phi_range.sample(rng);
        let dihedral_energy = dihedral_potential.energy(phi);
        if rng.next_f64() <= f64::exp(-beta * (angle_energy + dihedral_energy)) { break }
    }
    (theta, phi) // return angles
}

/// Returns a new position given the previous position with target bond length.
fn position_from_length<T: Rng>(
    position: &Vector3D,
    length: f64,
    rng: &mut T
) -> Vector3D {
    let mut normal = Normal::new(0.0, 1.0);
    // Create vector on unit sphere,
    // scale to desired length and translate.
    Vector3D::new(
            normal.sample(rng),
            normal.sample(rng),
            normal.sample(rng)
    ).normalized() * length + position
}

/// Returns a new position given the previous two positions in a chain molecule.
///
/// The resulting position will have a bond angle `theta`
/// w.r.t to the previous two positions and bond length `length` to the
/// previous position.
///
/// # Construction method (TODO make this more efficient/elegant)
/// 1. Compute a random vector lying on a cone surface.
///     * the random angle is uniformly picked from [0, 2PI]
///     * the cone angle is `theta` with respect to the z-axis
///     * the length of the vector is `length`
/// 2. Compute the rotation matrix that rotates the z-axis so that it aligns with
///    the vector that connects the previous positions.
/// 3. Rotate vector
/// 4. Translate vector to align with the rest of the molecule
fn position_from_angle<T: Rng>(
    positions: &[Vector3D],
    length: f64,
    theta: f64,
    rng: &mut T
) -> Vector3D {

    // checks that two positions are handed.
    assert!(positions.len() == 2);

    // Create a random point on a cone around pole with distance r.
    // The cone angle is the desired bond angle.
    let mut range = Range::new(0.0, 2.0 * PI);
    let phi = range.sample(rng);
    let r12 = Vector3D::new(
        length * f64::sin(theta) * f64::cos(phi),
        length * f64::sin(theta) * f64::sin(phi),
        length * f64::cos(theta)
    );

    let r10 = (positions[0] - positions[1]).normalized();
    // Compute the axis about which to rotate so that r10 and the pole
    // align.
    let e3 = Vector3D::new(0.0, 0.0, 1.0);
    let axis = r10 ^ e3;
    let angle = f64::acos(r10 * e3);
    let rotation = Matrix3::rotation(&axis, angle);
    // Perform rotation and shift position
    rotation * r12 + positions[1]
}

/// Returns a new position given the previous three positions in a chain molecule.
///
/// The new position is computed from a bonding length, bonding angle and
/// a dihedral angle.
fn position_from_angle_and_dihedral<T: Rng>(
    positions: &[Vector3D],
    length: f64,
    theta: f64,
    phi: f64,
    rng: &mut T
) -> Vector3D {

    assert!(positions.len() == 3);

    // The lazy way: create random position, obeying bond angle
    // we can replace this by constructing new_position obeying both angles
    let new_position = position_from_angle(&positions[1..3], length, theta, rng);

    // compute dihedral angle
    let r01 = positions[1] - positions[0];
    let r12 = positions[2] - positions[1];
    let r23 = new_position - positions[2];
    // assert!(r23.norm() == length);
    // assert!(f64::acos(-r12.normalized() * r23.normalized()) == theta);

    let u = r01 ^ r12;
    let v = r12 ^ r23;
    let angle = f64::atan2(r12.norm() * v * r01, u * v);
    // target angle: phi = angle + delta_phi
    let delta_phi = phi - angle;
    let rotation = Matrix3::rotation(& -r12, delta_phi);
    rotation * r23 + positions[2]
}

// Note: I use random (not really, since seed is always the same) input for tests.
// Maybe use quickcheck or fuzzing?
#[cfg(test)]
mod tests {
    use super::*;
    use types::Vector3D;
    use sys::{System,Particle};
    use energy::Harmonic;
    use std::f64;
    use std::f64::consts::PI;
    use rand::{XorShiftRng, SeedableRng, Rng};
    use consts::K_BOLTZMANN;

    #[test]
    fn test_create_bond_length_harmonic() {
        let mut rng = Box::new(XorShiftRng::new_unseeded());
        rng.reseed([2015u32, 42u32, 3u32, 12u32]);
        let beta = 1.0 / (K_BOLTZMANN * 50.0);
        let harmonic_potential = Harmonic{k: 62500. * K_BOLTZMANN, x0: 1.540};
        let iterations = 5000;
        let l: f64 = (0..iterations)
            .map(|_| create_bond_length_harmonic(beta, &harmonic_potential, &mut rng))
            .sum();
        assert!(f64::abs((l/iterations as f64 - harmonic_potential.x0)/harmonic_potential.x0) < 1e-3)
    }

    #[test]
    fn test_create_bond_angle_harmonic() {
        let mut rng = Box::new(XorShiftRng::new_unseeded());
        rng.reseed([2015u32, 42u32, 3u32, 12u32]);
        let beta = 1.0 / (K_BOLTZMANN * 50.0);
        let harmonic_potential = Harmonic{k: 62500. * K_BOLTZMANN, x0: f64::to_radians(114.0)};
        let iterations = 1000;
        let theta: f64 = (0..iterations).map(|_| create_bond_angle_harmonic(beta, &harmonic_potential, &mut rng)).sum();
        assert!(f64::abs((theta/iterations as f64 - harmonic_potential.x0)/harmonic_potential.x0) < 1e-3)
    }

    #[test]
    fn test_position_from_length() {
        let mut rng = Box::new(XorShiftRng::new_unseeded());
        rng.reseed([25u32, 42u32, 3u32, 12u32]);

        let position = Vector3D::new(1.0, 2.0, 3.0);
        // Create random bonding length
        let length = rng.next_f64();

        let new_position = position_from_length(&position, length, &mut rng);
        let r01 = new_position - position;
        assert_relative_eq!(r01.norm(), length, epsilon = 1.0e-14);
    }

    #[test]
    fn test_position_from_angle() {
        let mut rng = Box::new(XorShiftRng::new_unseeded());
        rng.reseed([25u32, 42u32, 3u32, 12u32]);

        let positions = vec!(
            Vector3D::new(17.0, 1.0, -3.7),
            Vector3D::new(-12.0, 0.6, 1.0)
        );

        // Create a random bonding angle.
        let mut range = Range::new(0.0, PI);
        let angle = range.sample(&mut rng);
        // Create random bonding length
        let length = rng.next_f64();

        // Create a new position from angle, length and old positions.
        let new_position = position_from_angle(&positions, length, angle, &mut rng);
        let r10 = positions[0] - positions[1];
        let r12 = new_position - positions[1];
        assert_relative_eq!(r12.norm(), length, epsilon = 1.0e-14);
        assert_relative_eq!(
            f64::acos(r10.normalized() * r12.normalized()),
            angle,
            epsilon = 1.0e-14
        );
    }

    #[test]
    fn test_position_from_angle_and_dihedral() {
        let mut rng = Box::new(XorShiftRng::new_unseeded());
        rng.reseed([25u32, 42u32, 3u32, 12u32]);

        let positions = vec!(
            Vector3D::new(1.2, 2.3, 3.4),
            Vector3D::new(2.4, 4.6, 6.8),
            Vector3D::new(3.6, 6.9, 9.0)
        );

        // Create a random bonding angle.
        let mut range = Range::new(0.0, PI);
        let theta = range.sample(&mut rng);
        let phi = range.sample(&mut rng);
        // Create random bonding length
        let length = rng.next_f64();
        // let length = 1.0;

        // Create a new position from angle, length and dihedral.
        let new_position = position_from_angle_and_dihedral(&positions, length, theta, phi, &mut rng);
        let r01 = positions[1] - positions[0];
        let r21 = positions[1] - positions[2];
        let r23 = new_position - positions[2];
        // Check length
        assert_relative_eq!(r23.norm(), length, epsilon = 1.0e-14);
        // Check bond angle
        assert_relative_eq!(
            f64::acos(r21.normalized() * r23.normalized()),
            theta,
            epsilon = 1.0e-14
        );
        // Check dihedral
        let u = r01 ^ (-r21);
        let v = (-r21) ^ r23;
        assert_relative_eq!(
            f64::atan2((-r21).norm() * v * r01, u * v),
            phi,
            epsilon = 1.0e-14
        );
    }

    // not yet a good test, but making it crash will print reasonable results.
    //
    // Some checks I did:
    // low temperature -> low derivation from set lengths/angles
    // high temperature -> high derivation from set lengths/angles
    // low force constant -> high derivation from set lengths/angles
    // high force constant -> low derivation from set lengths/angles
    // angles of 0°, 90°, 180° seem to work -> no numerical issues with trigonometric functions used.
    #[test]
    /// Builds a molecule, given an initial bead.
    fn test_trial_position() {
        let mut rng = Box::new(XorShiftRng::new_unseeded());
        rng.reseed([21u32, 42u32, 3u32, 12u32]);

        let beta = 1.0 / (K_BOLTZMANN * 300.0);

        let mut system = System::new();
        for _ in 0..5 {
            system.add_particle(Particle::new("C"))
        }
        for i in 0..4 {
            let _ = system.add_bond(i, i+1);
        }
        // set the first particle to some position
        system[0].position = Vector3D::new(1.0, 1.0, 1.0);
        assert_eq!(system.molecules().len(), 1);

        let (length, theta, phi) = (1.540, 114.0, 0.0);
        system.interactions_mut().add_bond(
            "C", "C", Box::new(Harmonic{k: 62500. * K_BOLTZMANN, x0: length}));
        system.interactions_mut().add_angle(
            "C", "C", "C", Box::new(Harmonic{k: 62500. * K_BOLTZMANN, x0: f64::to_radians(theta)}));
        system.interactions_mut().add_dihedral(
            "C", "C", "C", "C", Box::new(Harmonic{k: 62500. * K_BOLTZMANN, x0: f64::to_radians(phi)}));

        // Create trial positions from first bead
        let mut trial: Vec<Vector3D> = Vec::with_capacity(5);
        trial.push(system[0].position);
        for _ in 0..4 {
            let new_trial = trial_position(&system, system.molecule(0), &trial, beta, &mut rng);
            trial.push(new_trial);
        }

        // print all bonding lengths
        for i in 0..4 {
            println!("bond length: {:?}", (trial[i] - trial[i+1]).norm())
        }
        // print all bonding angles
        for i in 0..3 {
            let rji = (trial[i+1] - trial[i]).normalized();
            let rjk = (trial[i+1] - trial[i+2]).normalized();
            println!("bond angle: {:?}", f64::to_degrees(f64::acos(rji * rjk)))
        }
        // print all dihedral angles
        for i in 0..2 {
            let r01 = trial[i+1] - trial[i];
            let r12 = trial[i+2] - trial[i+1];
            let r23 = trial[i+3] - trial[i+2];
            let u = r01 ^ r12;
            let v = r12 ^ r23;
            println!("dihedral angle: {:?}", f64::to_degrees(f64::atan2(r12.norm() * v * r01, u * v)))
        }
        println!("built molecule: {:?}", trial);
        //assert_eq!(trial[0], trial[1]);
    }

    #[test]
    /// Performs 5000 (let n_step = 5000) selections and compares to weights.
    /// Fails if selections and weights differ in 1.0e-2.
    fn test_select_position() {
        let mut rng = Box::new(XorShiftRng::new_unseeded());
        rng.reseed([24u32, 42u32, 3u32, 12u32]);

        let mut counts = vec!(0, 0, 0);

        let a = Vector3D::new(1.2, 2.3, 3.4);
        let b = Vector3D::new(2.4, 4.6, 6.8);
        let c = Vector3D::new(3.6, 6.9, 9.0);
        let trialpos = vec!(a, b, c);
        let trialweights = vec!(0.2, 0.3, 0.5);
        let n_steps = 5000;
        for _ in 0..n_steps {
            let mut selection = select_position(&trialpos, &trialweights, &mut rng);
            if selection == a {
                counts[0] += 1
            } else if selection == b {
                counts[1] += 1
            } else if selection == c {
                counts[2] += 1
            }
        }
        println!("{:?}", counts);
        assert_relative_eq!(counts[0] as f64 / n_steps as f64, trialweights[0], epsilon = 1.0e-2);
        assert_relative_eq!(counts[1] as f64 / n_steps as f64, trialweights[1], epsilon = 1.0e-2);
        assert_relative_eq!(counts[2] as f64 / n_steps as f64, trialweights[2], epsilon = 1.0e-2)
    }
}

