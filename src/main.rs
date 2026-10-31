use std::path::PathBuf;
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

use clap::Parser;
use dfntrans::{DfntransError, read_input, track_particles, write_output_to_folder};
use tracing::{error, info};

mod libm_symbols;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to input file (TOML)
    input_file: PathBuf,

    /// Path to output folder
    output_folder: PathBuf,
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    init_tracing();

    info!("Starting dfntrans-rs");

    match run(&cli) {
        Ok(()) => {
            info!("dfntrans-rs - Complete");
            ExitCode::SUCCESS
        }
        Err(e) => {
            error!("{}", e);
            ExitCode::FAILURE
        }
    }
}

/// Reads the input, solves the flow with fehm-rs, tracks the particles and
/// writes the output files.
fn run(cli: &Cli) -> Result<(), DfntransError> {
    let mut input = read_input(&cli.input_file)?;
    if input.seed == 0 {
        input.seed = seed_from_clock();
        info!("Seed from the clock: {}", input.seed);
    }

    info!("Flow input: {}", input.flow_file.display());
    let flow = fehm::read_input(&input.flow_file)?;
    let solution = fehm::solve(&flow)?;
    info!(
        "Flow solved in {} Newton iterations",
        solution.newton_iterations
    );

    let transport = track_particles(&input, &flow, &solution)?;

    write_output_to_folder(&input, &transport, &cli.output_folder)?;

    Ok(())
}

/// Seconds since 1970 (as `time(0)` of DFNTrans), not 0
fn seed_from_clock() -> u32 {
    let seconds = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(1, |duration| duration.as_secs());
    (seconds as u32).max(1)
}

/// Logs to stdout. The level is `info`, unless the environment variable `RUST_LOG` gives another.
fn init_tracing() {
    let level_filter = tracing_subscriber::EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info"));

    tracing_subscriber::fmt()
        .compact()
        .with_env_filter(level_filter)
        .init();
}
