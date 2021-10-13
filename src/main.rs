#![feature(array_from_fn)]#![allow(uncommon_codepoints,mixed_script_confusables)]
use std::array::from_fn as map;
fn sq(x: f64) -> f64 { x*x }

fn main() {
	let cₛ = f64::sqrt(1./3.);
	let c = [-1., 0., 1.];
	let w = [1./6., 4./6., 1./6.];
	const N : usize = 800;
	let ν = 0.06;
	let τ = ν/sq(cₛ);
	let β = 1./(2.*τ+1.);
	let density = |x| if x<N/2 { 1.5 } else { 1.0 };
	let mut populations = map(|i| map(|x| w[i]*density(x)));
	let density = |populations:[[f64; N]; 3]| map(|x| populations.iter().map(|f| f[x]).sum());
	const T: usize = 500;
	for _ in 0..T {
		populations = map(|i| {
			let dx = i as isize - 1;
			map(|x| {
				let sx = x as isize + dx;
				if sx<0 || sx >= N as isize { let dx = -dx; let i = (1+dx) as usize; populations[i][x] } // Bounce back
				else { populations[i][sx as usize] } // Stream
			})
		});
		let density : [f64; N] = density(populations);
		let momentum : [f64; N] = map(|x| populations.iter().enumerate().map(|(i,f)| c[i]*f[x]).sum());
		populations = map(|i| map(|x| {
			let f = populations[i][x];
			f + 2.*β*(w[i]*(density[x] + c[i]/cₛ*momentum[x] + (sq(c[i])-sq(cₛ))/(2.*f64::powi(cₛ,4))*sq(momentum[x])/density[x]) - f)
		}));
	}
	assert_eq!(density(populations), [0.; N]);
}
