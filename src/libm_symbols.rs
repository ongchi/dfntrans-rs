//! The C math functions, from the `libm` crate.
//!
//! The math functions of the platforms round differently in the last digits,
//! and the times of the particles add up such differences. `f64::sin` and the
//! like call these C symbols, and the linker takes them from here and not from
//! the system. Thus each platform gives the output of WebAssembly, where the
//! Rust standard library has the same functions.

macro_rules! define_one_argument_functions {
    ($($name:ident),*) => {
        $(
            #[unsafe(no_mangle)]
            pub extern "C" fn $name(value: f64) -> f64 {
                libm::$name(value)
            }
        )*
    };
}

macro_rules! define_two_argument_functions {
    ($($name:ident),*) => {
        $(
            #[unsafe(no_mangle)]
            pub extern "C" fn $name(first: f64, second: f64) -> f64 {
                libm::$name(first, second)
            }
        )*
    };
}

define_one_argument_functions!(
    sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh, exp, exp2, expm1, log,
    log2, log10, log1p, cbrt
);
define_two_argument_functions!(atan2, pow, hypot);

/// `f32::acos` of the velocity reconstruction.
#[unsafe(no_mangle)]
pub extern "C" fn acosf(value: f32) -> f32 {
    libm::acosf(value)
}

/// LLVM joins `sin` and `cos` of one value into this function on Linux.
///
/// # Safety
///
/// `sine` and `cosine` must be valid for writes.
#[cfg(target_os = "linux")]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn sincos(value: f64, sine: *mut f64, cosine: *mut f64) {
    let (sine_value, cosine_value) = libm::sincos(value);
    unsafe {
        *sine = sine_value;
        *cosine = cosine_value;
    }
}

#[cfg(target_os = "macos")]
#[repr(C)]
pub struct SineCosine {
    sine: f64,
    cosine: f64,
}

/// LLVM joins `sin` and `cos` of one value into this function on macOS.
#[cfg(target_os = "macos")]
#[unsafe(no_mangle)]
pub extern "C" fn __sincos_stret(value: f64) -> SineCosine {
    let (sine, cosine) = libm::sincos(value);
    SineCosine { sine, cosine }
}
