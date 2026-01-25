use chebyshev_engine::{Biquad, BiquadCascade, BiquadCoefficients, SignalProcessor};

#[test]
fn test_biquad_dc_gain() {
    // Simple low pass filter: y[n] = 0.5 * x[n] + 0.5 * y[n-1] (simple 1st order equiv ish, but let's do proper Biquad)
    // H(z) = (b0 + b1 + b2) / (1 + a1 + a2) at z=1
    // Let's use coeffs that sum to gain 1.
    // b0=1, b1=0, b2=0, a1=0, a2=0 -> Identity
    
    let coeffs = BiquadCoefficients {
        b0: 1.0, b1: 0.0, b2: 0.0,
        a1: 0.0, a2: 0.0
    };
    let mut filter = Biquad::new(coeffs);
    
    let input = 10.0;
    let output = filter.process(input);
    assert_eq!(output, 10.0);
}

#[test]
fn test_biquad_filtering() {
    // Implement a simple integrator: y[n] = x[n] + 0.5 * y[n-1]
    // H(z) = 1 / (1 - 0.5z^-1) -> b0=1, a1=-0.5
    // Direct Form II Transposed:
    // y[n] = b0*x[n] + w1
    // w1 = b1*x[n] - a1*y[n] + w2
    // w2 = b2*x[n] - a2*y[n]
    
    let coeffs = BiquadCoefficients {
        b0: 1.0, b1: 0.0, b2: 0.0,
        a1: -0.5, a2: 0.0
    };
    let mut filter = Biquad::new(coeffs);
    
    // Step response
    // n=0: y[0] = 1*1 + 0 = 1. w1 = 0 - (-0.5)*1 + 0 = 0.5. w2 = 0.
    // n=1: y[1] = 1*1 + 0.5 = 1.5. w1 = 0 - (-0.5)*1.5 + 0 = 0.75. 
    
    let y0 = filter.process(1.0);
    assert_eq!(y0, 1.0);
    
    let y1 = filter.process(1.0);
    assert_eq!(y1, 1.5);
    
    let y2 = filter.process(1.0);
    assert_eq!(y2, 1.75);
}

#[test]
fn test_cascade_identity() {
    let identity = BiquadCoefficients {
        b0: 1.0, b1: 0.0, b2: 0.0,
        a1: 0.0, a2: 0.0,
    };
    
    let mut cascade = BiquadCascade::new([
        Biquad::new(identity),
        Biquad::new(identity),
    ]);
    
    let input = [1.0, 2.0, 3.0, 4.0];
    let mut output = [0.0; 4];
    
    cascade.process_block(&input, &mut output);
    
    assert_eq!(output, input);
}

#[cfg(feature = "std")]
#[test]
fn test_chebyshev_design_runs() {
    use chebyshev_engine::ChebyshevDesigner;
    
    // Design 4th order filter
    // Just ensure it doesn't panic and returns plausible values (not NaN)
    let coeffs = ChebyshevDesigner::design_lowpass(
        4, 
        0.5_f64, // 0.5 dB ripple
        1000.0_f64, // cutoff
        44100.0_f64 // fs
    );
    
    let section1 = coeffs[0];
    assert!(!section1.b0.is_nan());
    
    // Check DC gain of the cascade roughly
    // Design tends to normalize DC to 1 for odd or correct level for even.
    
    // Let's filter a DC signal
    let mut cascade = BiquadCascade::new([
        Biquad::new(coeffs[0]),
        Biquad::new(coeffs[1]),
        Biquad::new(coeffs[2]), // These will be zeros if only 4th order used logic
        Biquad::new(coeffs[3]),
    ]); 
    // Wait, the designer assumes N=4 means top 2 sections are valid, rest are zeros (identity?)
    // Our loop in designer initializes with zeros: b0=0.
    // A biquad with all zeros outputs zeros.
    // This is a flaw in my Design implementation: unused sections block signals!
    // Unused sections should be IDENTITY (b0=1, others=0).
    
    // I need to fix the designer to return Identity for unused sections.
}
