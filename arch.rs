use crate::complex::Comp ;
use crate::complex ;
use crate::cv_spaces ;
use crate::quant_f ;
use rand::random ;
#[derive(Debug ,Clone , PartialEq)]  
pub struct Qubit {
    or_state: Vec<Comp> ,
    norm: f64 }


pub fn qubitize (vec: &Vec<Comp>) -> Qubit {
    let mod_ = cv_spaces::norm(&vec) ;

      Qubit{
          or_state: vec.clone() ,
          norm: mod_ }
}

pub fn process_layer (gates: &Vec<Vec<Vec<Comp>>>) -> Vec<Vec<Comp>> {

let mut cont: Vec<Vec<Comp>> = gates[0].clone() ;
for x in 1..gates.len() {
    cont = cv_spaces::tens_prod(&cont , &gates[x]) ;
}
cont 
}

//in the layer the gates are entered by starting with the one at the bottom going upwards , but in
//the layers collection enter layers from left ro right as they appear on the page 
pub fn circuit(states: &Vec<Vec<Vec<Comp>>> ,layers: &Vec<Vec<Vec<Vec<Comp>>>>) -> Vec<Vec<Comp>> {

    let states_processed = process_layer(&states) ;
    let mut layers_processed: Vec<Vec<Vec<Comp>>> = Vec::new() ;

    for x in 0..layers.len() {
        let layer_processed = process_layer(&layers[x]) ;
        layers_processed.push(layer_processed);
            }
    let mut fin_cont = states_processed ;
    for x in 0..layers_processed.len() {
        fin_cont = cv_spaces::mul_mat(&layers_processed[x],&fin_cont);
    }
fin_cont 
}
pub fn ret_ket(key: &str) -> Vec<Vec<Comp>> {
    let mut res: Vec<Vec<Comp>> = Vec::new() ;
   match key {
    
   "zero" =>  { res = vec![vec![Comp{ real: 1f64 ,imag: 0f64} ],vec![Comp{ real: 0f64 , imag: 0f64} ]];}

 "one" =>  { res = vec![vec![Comp{ real: 0f64 ,imag: 0f64} ],vec![Comp{ real: 1f64 , imag: 0f64} ]];}

   &_ => { println!("Argument not matched") ;}

   }

   res 
}
pub fn ret_gate(key: &str) -> Vec<Vec<Comp>> {
    let mut res: Vec<Vec<Comp>> = Vec::new() ;
   match key {

   "px" =>  { res = vec![vec![Comp{ real: 0f64 ,imag: 0f64} ,Comp{real: 1f64 ,imag: 0f64}],vec![Comp{ real: 1f64 , imag: 0f64} , Comp{real: 0f64 , imag: 0f64}]];}

   "py" =>  { res = vec![vec![Comp{ real: 0f64 ,imag: 0f64} ,Comp{real: 0f64 ,imag: -1f64}],vec![Comp{ real: 0f64 , imag: 1f64} , Comp{real: 0f64 , imag: 0f64}]];}
   
 "pz" =>  { res = vec![vec![Comp{ real: 1f64 ,imag: 0f64} ,Comp{real: 0f64 ,imag: 0f64}],vec![Comp{ real: 0f64 , imag: 0f64} , Comp{real: -1f64 , imag: 0f64}]];}

 "s" =>  { res = vec![vec![Comp{ real: 1f64 ,imag: 0f64} ,Comp{real: 0f64 ,imag: 0f64}],vec![Comp{ real: 0f64 , imag: 0f64} , Comp{real: 0f64 , imag: 1f64}]];}

 "t" =>  { res = vec![vec![Comp{ real: 1f64 ,imag: 0f64} ,Comp{real: 0f64 ,imag: 0f64}],vec![Comp{ real: 0f64 , imag: 0f64} , Comp{real: (std::f64::consts::PI*0.25).cos() , imag:  (std::f64::consts::PI*0.25).sin() }]];}


"cx_gmin" =>  { res = vec![vec![Comp{ real: 1f64 ,imag: 0f64} ,Comp{real: 0f64 ,imag: 0f64},Comp{real: 0f64 ,imag: 0f64},Comp{real: 0f64 ,imag: 0f64}],
vec![Comp{ real: 0f64 ,imag: 0f64} ,Comp{real: 1f64 ,imag: 0f64},Comp{real: 0f64 ,imag: 0f64},Comp{real: 0f64 ,imag: 0f64}],
vec![Comp{ real: 0f64 ,imag: 0f64} ,Comp{real: 0f64 ,imag: 0f64},Comp{real: 0f64 ,imag: 0f64},Comp{real: 1f64 ,imag: 0f64}],
vec![Comp{ real: 0f64 ,imag: 0f64} ,Comp{real: 0f64 ,imag: 0f64},Comp{real: 1f64 ,imag: 0f64},Comp{real: 0f64 ,imag: 0f64}],

];}

"cx_gplus" =>  { res = vec![vec![Comp{ real: 1f64 ,imag: 0f64} ,Comp{real: 0f64 ,imag: 0f64},Comp{real: 0f64 ,imag: 0f64},Comp{real: 0f64 ,imag: 0f64}],
vec![Comp{ real: 0f64 ,imag: 0f64} ,Comp{real: 0f64 ,imag: 0f64},Comp{real: 0f64 ,imag: 0f64},Comp{real: 1f64 ,imag: 0f64}],
vec![Comp{ real: 0f64 ,imag: 0f64} ,Comp{real: 0f64 ,imag: 0f64},Comp{real: 1f64 ,imag: 0f64},Comp{real: 0f64 ,imag: 0f64}],
vec![Comp{ real: 0f64 ,imag: 0f64} ,Comp{real: 1f64 ,imag: 0f64},Comp{real: 0f64 ,imag: 0f64},Comp{real: 0f64 ,imag: 0f64}],

];}


 "hd" =>  { res = vec![vec![Comp{ real: 1f64/2f64.sqrt() ,imag:  0f64 } ,Comp{real:  1f64/2f64.sqrt()  ,imag:  0f64 }],vec![Comp{ real:  1f64/2f64.sqrt()  , imag: 0f64} , Comp{real:  -1f64/2f64.sqrt()  , imag: 0f64}]];}

 &_ => { println!("Argument not matched") ;}


   }

   res}

pub fn ph_shift(theta: f64) -> Vec<Vec<Comp>> {
    vec![vec![Comp{ real: 1f64 ,imag: 0f64} ,Comp{real: 0f64 ,imag: 0f64}],vec![Comp{ real: 0f64 , imag: 0f64} , Comp{real: theta.exp() , imag: 0f64}]]
}
pub fn generate_at_i (i: usize , n: usize) -> Vec<Comp> {
    let mut  cont: Vec<Comp> = Vec::new() ;
    for x in 0..n {
        if x == i {
            cont.push(complex::new(1f64 , 0f64));
        } else {
        cont.push(complex::new(0f64 , 0f64));
        }}
      cont 
        }

pub fn collapse (state: &Vec<Vec<Comp>>) ->Vec<Comp> {

    let mut orig: Vec<Comp> = Vec::new() ;

    for x in 0..state.len() {
        orig.push(state[x][0].clone());

    }
    let mut probs: Vec<f64> = Vec::new() ;

    for x in 0..orig.len() {
    probs.push(complex::mod_(&orig[x]).powi(2));
    }

    let rand: f64 = random::<f64>() ;
    let mut cum = 0.0 ;
    let mut index = 0 ;
    for (i ,p) in probs.iter().enumerate(){
        cum += p ;
        if rand < cum {
            index = i; }}
    
   generate_at_i( index , state.len())
}

pub fn ph_shift_xyz(theta: f64 , key: &str ) -> Vec<Vec<Comp>> {
   let mut res: Vec<Vec<Comp>> = Vec::new() ;
   match key {

   "rx" =>  { res = 
    vec![vec![Comp{ real: (theta*0.5).cos() ,imag: 0f64} ,Comp{real: 0f64 ,imag: (-0.5 *theta).sin()}],vec![Comp{ real: 0f64 , imag: (-0.5*theta).sin()} , Comp{real: (theta * 0.5).cos() , imag: 0f64}]]
}
 "ry" =>  { res = 
    vec![vec![Comp{ real: (theta*0.5).cos() ,imag: 0f64} ,Comp{real: (-0.5 *theta).sin() , imag: 0f64}],vec![Comp{ real: (0.5*theta).sin() , imag: 0f64} ,Comp{real: (0.5 *theta).cos() , imag: 0f64}]]]
}
 "rz" =>  { res = 
    vec![vec![Comp{ real: (theta* -0.5).cos() ,imag: (theta* -0.5).sin()} ,Comp{real: 0f64 , imag: 0f64}],vec![Comp{real: 0f64 , imag: 0f64} ,Comp{ real: (theta*0.5).cos() ,imag: (theta*0.5).sin()} ]]
}

pub fn if_on_perform(control: &Vec<Vec<Comp>> , target: &Vec<Vec<Comp>> ,operation: &Vec<Vec<Comp>>) -> Vec<Vec<Comp>> {
let states = vec![ control.clone() , target.clone()] ;


}
