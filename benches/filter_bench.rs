use chebyshev_engine::{Biquad, BiquadCascade, BiquadCoefficients, NaiveFilter, SignalProcessor};
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_filter(c: &mut Criterion) {
    let mut group = c.benchmark_group("Filter Performance");
    
    // Sample data
    let input: Vec<f32> = (0..1024).map(|i| (i as f32).sin()).collect();
    let mut output = vec![0.0f32; 1024];

    // Setup 4th order (2 biquads)
    let bcoeffs = BiquadCoefficients {
        b0: 0.0039, b1: 0.0078, b2: 0.0039,
        a1: -1.8153, a2: 0.8310,
    };
    let mut cascade = BiquadCascade::new([
        Biquad::new(bcoeffs),
        Biquad::new(bcoeffs),
    ]);

    // Naive 4th order
    let b = [0.0039, 0.0078, 0.0039, 0.0, 0.0]; // simplified for example
    let a = [1.0, -1.8153, 0.8310, 0.0, 0.0];
    let mut naive = NaiveFilter::new(b, a);

    group.bench_function("BiquadCascade (Block)", |b| {
        b.iter(|| {
            cascade.process_block(black_box(&input), black_box(&mut output));
        })
    });

    group.bench_function("NaiveFilter (Single Sample)", |b| {
        b.iter(|| {
            for (i, o) in input.iter().zip(output.iter_mut()) {
                *o = naive.process(black_box(*i));
            }
        })
    });

    group.finish();
}

criterion_group!(benches, bench_filter);
criterion_main!(benches);
