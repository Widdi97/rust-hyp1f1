// based on https://people.maths.ox.ac.uk/porterm/papers/hypergeometric-final.pdf
// use statrs::function::gamma;



// lanczos approximation of the gamma function as shown in https://en.wikipedia.org/wiki/Lanczos_approximation


fn gamma (z:f64) -> f64 {
    let p: [f64;8] = [676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
    const PI: f64 = std::f64::consts::PI;
    let y: f64;
    if z < 0.5 {
        y = PI / ( (PI * z).sin() * gamma(1.0 - z) );  // Reflection formula
    } else {
        let z2: f64 = z - 1.0;
        let mut x: f64 = 0.99999999999980993;
        for ii in 0..8 {
            x += p[ii] / (z2 + ii as f64 + 1.0);
        }
            
        let t: f64 = z2 + 8.0 - 0.5;
        y = f64::powf(2.0 * PI, 0.5) * f64::powf(t , z2 + 0.5) * f64::powf(std::f64::consts::E,-t) * x;
    }
    y
}

fn poch (z:f64,k:i64) -> f64 {
    if k == 0 {
        1.0
    } else {
        gamma(z + k as f64) / gamma(z)
    }
}

fn hyp1f1_slow (a:f64,b:f64,z:f64) -> f64 { // inefficient
    let mut result = 0.0;
    const NMAX: i64 = 100;
    const EPS: f64 = 2.2e-16;
    let mut a_j: f64;
    for j in 0..=NMAX {
        a_j = poch(a,j) / poch(b,j) / gamma((j+1) as f64) * f64::powi(z, j as i32);
        result += a_j;
        if (a_j.abs() / result.abs()) < EPS {
            //println!("{}", j);
            break;
        }
    }
    result
}

fn hyp1f1_series (a:f64,b:f64,z:f64) -> f64 { // taken from scipy source code https://github.com/scipy/scipy/blob/main/scipy/special/_hypergeometric.pxd
    let mut k: f64;
    let mut term: f64 = 1.0;
    let mut result = 1.0;
    const EPS: f64 = 2.2e-16;
    let mut kfloat;
    for k in 0..500 {
        kfloat = k as f64;
        term *= (a + kfloat) * z / (b + kfloat) / (kfloat + 1.0);
        result += term;
        if term.abs() <= EPS * result.abs() {
            break;
        }
    }

    result
}




fn main() {
    println!("{}",hyp1f1(2.5,5.33,6.4));

    let now = std::time::Instant::now();

    for _ii in 0..10000 {
        hyp1f1(2.5,5.33,6.4);
    }

    let elapsed_time = now.elapsed();
    println!("Running slow_function() took {} ms.", elapsed_time.as_millis());
}
