//! Matrix diffusion with the time domain random walk (DFNTrans
//! `TimeDomainRW`; Hyman et al. 2019, Geophysical Research Letters 46,
//! 13785-13795).

use crate::network::PI;
use crate::random::Drand48;

/// Parameters of the matrix diffusion in the units of the run
#[derive(Clone, Copy, Debug)]
pub struct MatrixDiffusion {
    pub porosity: f64,
    /// m² per time unit
    pub diffusivity: f64,
}

/// The time which a particle stays in the matrix during an advective time
/// `advective_time` in a fracture of aperture `aperture`:
/// (φ·√D·t / (b · erfc⁻¹(1 − r)))², with a random number r
pub fn diffusion_time(
    parameters: MatrixDiffusion,
    aperture: f64,
    advective_time: f64,
    random: &mut Drand48,
) -> f64 {
    let rate = parameters.porosity * parameters.diffusivity.sqrt() / aperture;
    let z = 1.0 - random.next_f64();
    (rate * advective_time / inverse_erfc(z)).powi(2)
}

/// The series of the inverse complementary error function of DFNTrans:
/// the inverse error function of `z` (erfc⁻¹(1 − z) = erf⁻¹(z)), to the
/// power 11
pub fn inverse_erfc(z: f64) -> f64 {
    0.5 * PI.sqrt()
        * (z + (PI / 12.0) * z.powi(3)
            + (7.0 * PI.powi(2) / 480.0) * z.powi(5)
            + (127.0 * PI.powi(3) / 40320.0) * z.powi(7)
            + (4369.0 * PI.powi(4) / 5_806_080.0) * z.powi(9)
            + (34807.0 * PI.powi(5) / 182_476_800.0) * z.powi(11))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// erf⁻¹ of tabulated values; the error of the series grows with z
    #[test]
    fn series_gives_the_inverse_error_function() {
        for (z, expected, tolerance) in [
            (0.1, 0.088_855_990_494_257_7, 1e-12),
            (0.5, 0.476_936_276_204_469_9, 1e-5),
        ] {
            assert!((inverse_erfc(z) - expected).abs() < tolerance, "{z}");
        }
    }

    #[test]
    fn diffusion_time_grows_with_the_square_of_the_advective_time() {
        let parameters = MatrixDiffusion {
            porosity: 0.02,
            diffusivity: 1e-13,
        };
        let short = diffusion_time(parameters, 1e-4, 1.0, &mut Drand48::new(5));
        let long = diffusion_time(parameters, 1e-4, 2.0, &mut Drand48::new(5));
        assert!((long / short - 4.0).abs() < 1e-12);
    }
}
