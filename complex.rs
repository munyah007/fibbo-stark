#[derive(Debug , Copy , Clone, PartialEq)]
pub struct Comp { pub real: f64 ,
                  pub imag: f64 }
#[derive(Debug , Copy , Clone)]
pub struct CompPolar { mod_: f64 ,
                       theta: f64 }
pub fn normalize (x: f64) -> f64 {
    let eps = 1e-12 ;
    if (x - 1.0).abs() < eps {
        1.0}
    else if x.abs() < eps{
        0.0 }
    else {
        x}}

pub fn new(r: f64 ,i: f64) ->Comp {
    Comp{real: normalize(r) ,
         imag: normalize(i) ,
    }}

pub fn add (comp_1: &Comp , comp_2: &Comp) -> Comp{
    Comp{ real: (comp_1.real + comp_2.real),
          imag: (comp_1.imag + comp_2.imag ),
    }}

pub fn sub (comp_1: &Comp , comp_2: &Comp) -> Comp{
    Comp{ real: (comp_1.real - comp_2.real),
          imag: (comp_1.imag - comp_2.imag ),
    }}





pub fn mul (comp_1: &Comp , comp_2: &Comp) -> Comp{
    let a1a2 =  comp_1.real * comp_2.real ;
    let b1b2 =  comp_1.imag * comp_2.imag ;

    let real = a1a2 - b1b2 ;

    // imaginary part components 
    let a1b2 =  comp_1.real * comp_2.imag;
    let a2b1 =  comp_1.imag * comp_2.real ;
    let imag = a1b2 + a2b1 ;

   new(real , imag ) 
}

pub fn mod_( comp: &Comp) -> f64 {
    (comp.real.powi(2) + comp.imag.powi(2)).sqrt() 
}

pub fn div (comp_1: &Comp , comp_2: &Comp) -> Comp{
    let mod_den = mod_(comp_2) ;
    let a1a2 =  comp_1.real * comp_2.real ;
    let b1b2 =  comp_1.imag * comp_2.imag ;

    let real = (a1a2 + b1b2)/mod_den ;

    // imaginary part components 
    let a1b2 =  comp_1.real * comp_2.imag;
    let a2b1 =  comp_1.imag * comp_2.real ;
    let imag =( a2b1 - a1b2)/mod_den ;

   new(real , imag ) 
}




pub fn conj(comp: &Comp)-> Comp {
    new(comp.real , -1f64 * comp.imag )
      }

pub fn inv(comp: &Comp)-> Comp {
    new(-1f64*comp.real , -1f64 * comp.imag )
      }

pub fn to_polar (comp: &Comp ) -> CompPolar {
let theta = comp.imag.atan2(comp.real);
let  mod_ = mod_(comp) ;
         CompPolar{ mod_: mod_ ,
           theta: theta, 
             } }

pub fn to_cart (comp: &CompPolar ) -> Comp {
// e^{itheta} = cos(theta) + i sin(theta)
let real =  comp.theta.cos()*comp.mod_ ;
let imag = comp.theta.sin()*comp.mod_ ;
         Comp{ real:  real ,
               imag:  imag, 
             } }

pub fn polar_mul (comp1: &CompPolar , comp2: &CompPolar ) -> CompPolar {

         CompPolar{ mod_: comp1.mod_ * comp2.mod_ ,
           theta: comp1.theta + comp2.theta, 
             } }
pub fn polar_pow_n (comp1: &CompPolar , n: i32 ) -> CompPolar {

         CompPolar{ mod_: comp1.mod_.powi(n) ,
           theta: comp1.theta* n as f64 , 
             } }

pub fn polar_root_nth (comp1: &CompPolar , n: i32 ) -> CompPolar {
        let root = 1f64/n as f64 ;
// note below , the other roots have theta as the current calculated theta + root*k*(2pi) running
// all the way from k =  0 to k = n - 1 , after that the cycle starts again 
// all complext numbers with the same modulus and angles theta + k*2pi are identical because they
// are cyclic in a mod 2pi group , the above are called roots of unity because we set the modulus
// to be equal to 1 
         CompPolar{ mod_: comp1.mod_.powf(root) ,
           theta: comp1.theta*root , 
             } }
//  Complex vector spaces functionality 
