
use lagrangian_interpolation::lagrange_interpolate ;
#[derive(Debug)]
pub struct Omega {
   pub coef: i32,
   pub power: u32 , }


#[derive(Debug)]
pub struct Poly {
    coef: i32 ,    power: u32, }

const FIELD: u32  =  97 ;//99991 ;
const PRIM_ROOT: u32 = 5 ; //6 ;

// length of trace should always be equal - ie a power of 2 length 
// for now i will work with 8 
const TRACE_LENGTH:  u32 = 8 ; // this is k 
pub fn  newOmega (power: u32 , coef: i32) -> Omega {
       Omega{ power , coef }
}

pub fn newPoly(power:u32 , coef: i32) -> Poly {
        Poly{ power , coef }
}

pub fn change_poly_to_omega (poly: &Poly , k: u32) -> Omega {
        Omega{ power: (poly.power/k) , coef: poly.coef }
}
// functionality for power 
pub mod powers {
trait One {
    fn one() -> Self ;
}
impl One for u32 {
    fn one () -> Self {1}}

impl One for f64 { 
    fn one() -> Self {1.0}}

pub fn pow_ <T> (mut base: T, mut exp: u32) -> T 
                                          where 
                                            T: Copy + std::ops::Mul<Output = T> + One,
{ let mut result = T::one();
    while exp > 0 {
        if exp % 2 == 1 {
              result = result * base ;
        }
        base = base * base ;
        exp /=2 ;
    }
result 

}
}

// later do fn mod_pow (base: u64 , exp: u64 , modulus: u64) -> u64 
// for cryptography ie here needed 
//
// Produce xs 


////////////////////
pub fn eval_Omega (om: &Omega , order: u32) -> u32 {
    let subprim: u32 =  (PRIM_ROOT.pow((FIELD - 1)/order)) % FIELD
;
     let coef_p: u32 = if om.coef < 0 { (om.coef * -1)as u32 } else {om.coef as u32 } ;
     let valPen: u32 =  (coef_p * (subprim.pow( om.power)) ) % FIELD;
     let val: u32 = if om.coef < 0 {
          FIELD -  valPen }
           else { valPen  };
     val }


pub fn prod_sub_gr (order: u32 , start: u32) -> Vec<Omega>{
    let mut cont: Vec<Omega> = Vec::new();
    if start == 0 {
        for x in 0..order {
            cont.push(Omega{power: x as u32 ,coef: 1})}}
        else { for x in 1..order+1 {cont.push(Omega{power: x , coef: 1})}}

            cont }


pub fn row_x_omega_eval (sub: &Vec<Omega>, coefs: &Vec<i32> ,k: u32 ,rown: u32) ->i32 {
    let track: Vec<usize> = (0..=sub.len() - 1).collect() ;
    let subprim: u32 = FIELD % (powers::pow_(PRIM_ROOT , (FIELD - 1)/k));
    let row: Vec<Omega> = (sub.iter()).map(|x| {Omega{coef: x.coef,
                                                      power: x.power*rown}}).map(|x| x.process_bigger_than_kov2(k)).collect();
    let mut value: i32 = 0 ;
    for x in track{
         value = value + coefs[x]*((powers::pow_(subprim , row[x].power)) as i32)*(row[x].coef) }
             value }
pub fn process (  k: u32 , cont: &mut Vec<Vec<Omega>>) {
    
     let steps = (( k as f64).log2()) as usize ;
     let proc: Vec<usize> = (1..=steps).collect();
      for x in proc {
         let mut imCont: Vec<Omega> = Vec::new() ;
         
          for y in &cont[x-1]{
                  let length: usize = cont[x-1].len();
                  y.sqrt(&mut imCont , k);
                  
          }
          cont.push(imCont)
      }

}

impl Omega {
// process at the very last 
pub fn process_bigger_than_kov2 (&self , k: u32) -> Omega {
    if self.power < (k/2) {
        Omega{ power: self.power , coef: self.coef *(1) , }
           } 
    else if self.power == (k/2) {
         Omega{ power: 0 , coef:  self.coef*(-1) , }
    }
    else { Omega{ power:  (self.power - (k/2)) , coef: self.coef*(-1) }}}

pub fn sqrt (&self ,cont: &mut Vec<Omega> , k: u32) -> (){
    
    if self.power == 0 {
        cont.push(Omega{ coef: 1 ,
                         power: 0}) ;
        cont.push(Omega{ coef: 1,
                         power:  ((self.power/2) + (k/2) )});}
        else {  cont.push(Omega{ coef: 1 ,
                         power: (self.power / 2 )}) ;
        cont.push(Omega{ coef: 1,
                         power: ((self.power/2) + (k/2) )});
    } ;

        
}

}   

