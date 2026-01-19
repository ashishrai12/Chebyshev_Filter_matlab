use chebyshev_engine::{Biquad, BiquadCascade, BiquadCoefficients, SignalProcessor};
use hound;
use std::env;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        println!("Usage: chebyshev-cli <input.wav> <output.wav>");
        return Ok(());
    }

    let input_path = &args[1];
    let output_path = &args[2];

    println!("Reading {}...", input_path);
    let mut reader = hound::WavReader::open(input_path)?;
    let spec = reader.spec();

    // Setup a 4th order Chebyshev Type I Low-pass filter (0.5 dB ripple)
    // Cutoff frequency at 0.1 * Fs (normalized)
    let sample_rate = spec.sample_rate as f32;
    let cutoff = sample_rate * 0.1; 
    let order = 4;
    let ripple = 0.5;

    println!("Designing filter: Order={}, Cutoff={}Hz, Ripple={}dB, Fs={}Hz", order, cutoff, ripple, sample_rate);

    let all_coeffs = chebyshev_engine::ChebyshevDesigner::design_lowpass(
        order,
        ripple,
        cutoff,
        sample_rate
    );
     
    // Create cascade from the first 2 sections (4th order = 2 biquads)
    // The designer returns fixed size 8 array, we take what we need.
    // Note: We need to properly initialize the sections array for the Cascade.
    // Since BiquadCascade<T, 2> expects [Biquad<T>; 2], we construct it manually.
    
    let section1 = Biquad::new(all_coeffs[0]);
    let section2 = Biquad::new(all_coeffs[1]);

    let mut cascade = BiquadCascade::new([section1, section2]);

    println!("Processing audio...");
    let mut writer = hound::WavWriter::create(output_path, spec)?;

    match spec.sample_format {
        hound::SampleFormat::Float => {
            for sample in reader.samples::<f32>() {
                let s = sample?;
                let filtered = cascade.process(s);
                writer.write_sample(filtered)?;
            }
        }
        hound::SampleFormat::Int => {
            let max_val = (1 << (spec.bits_per_sample - 1)) as f32;
            for sample in reader.samples::<i32>() {
                let s = sample? as f32 / max_val;
                let filtered = cascade.process(s);
                writer.write_sample((filtered * max_val) as i32)?;
            }
        }
    }

    writer.finalize()?;
    println!("Done! Output saved to {}", output_path);

    Ok(())
}
