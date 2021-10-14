#![feature(array_from_fn)]#![allow(uncommon_codepoints,mixed_script_confusables)]
use std::array::from_fn as map;
fn sq(x: f64) -> f64 { x*x }

fn main() -> Result<(), impl std::fmt::Debug> {
	let cₛ2 = 1./3.;
	let c = [-1., 0., 1.];
	let w = [1./6., 4./6., 1./6.];
	let local_equilibrium = |i, density, momentum| w[i]*(density + c[i]/cₛ2*momentum + (sq(c[i])-cₛ2)/(2.*sq(cₛ2))*sq(momentum)/density);
	const N : usize = 800;
	//let ν = 0.06;
	let ν = 5e-8;
	let τ = ν/cₛ2;
	let β = 1./(2.*τ+1.);
	//let density = |x| if x<N/2 { 1.5 } else { 1.0 };
	let density: [f64; N] = map(|x| 1. + 1./2. * f64::exp(-5000.*sq(x as f64/N as f64 - 1./4.)));
	let momentum: [f64; N] = map(|x| 0.1*density[x]);
	let mut populations = map(|i| map(|x| local_equilibrium(i, density[x], momentum[x])));
	let density = |populations:[[f64; N]; 3]| map(|x| populations.iter().map(|f| f[x]).sum());
	const T: usize = 64;
	for _ in 0..T {
		populations = map(|i| {
			let dx = i as isize - 1;
			map(|x| {
				let sx = x as isize - dx;
				//if sx<0 || sx >= N as isize { let dx = -dx; let i = (1+dx) as usize; populations[i][x] } // Bounce back
				if sx<0 { populations[i][N-1] } // Periodic
				else if sx as usize>=N { populations[i][0] } // Periodic
				else { populations[i][sx as usize] } // Stream
			})
		});
		let density: [f64; N] = density(populations);
		let momentum: [f64; N] = map(|x| populations.iter().enumerate().map(|(i,f)| c[i]*f[x]).sum());
		populations = map(|i| map(|x| {
			let f = populations[i][x];
			f + 2.*β*(local_equilibrium(i, density[x], momentum[x]) - f)
		}));
	}
	let density = density(populations);
	let values: [_; N] = map(|x| (x as f64, [[density[x]].into()].into()));
	ui::app::run(ui::plot::Plot::new(&[&["ρ"]], &values))
}
