use crate::complex::Comp ;
use crate::complex ;
use crate::cv_spaces ;
#[derive(Debug ,Clone , PartialEq)]  
pub struct States_N_Prob {
    fr_state: Vec<Comp> ,
    to_states: Vec<Vec<Comp>> ,
    probs:  Vec<Comp> ,
    omega: Vec<Vec<Comp>> ,//Omega_psi
    mean: f64 ,
}

pub fn n_clicks (curr_state: &Vec<Comp> , trans_mat: &Vec<Vec<Comp>> , clicks: usize) -> Vec<Comp> {
    let mut matricised: Vec<Vec<Comp>> = Vec::new() ;
    for x in 0..curr_state.len() {
        let new_row: Vec<Comp> = vec![curr_state[x].clone()] ;
        matricised.push(new_row) ;
    }

    let mut new_trans_mat: Vec<Vec<Comp>> = trans_mat.clone() ;
    for x in 1..clicks{
        new_trans_mat = cv_spaces::mul_mat(&new_trans_mat, &trans_mat) ;
    }

    let final_state = cv_spaces::mul_mat(&new_trans_mat ,&matricised) ;
    let mut dematricised: Vec<Comp> = Vec::new() ;
    for x in 0..final_state.len() {
         dematricised.push(final_state[x][0].clone()) ;
    }

    dematricised 
}

pub fn obs_prob(curr_ket: &Vec<Comp> , index: usize) -> f64{

    let mod_at_index = complex::mod_(&curr_ket[index]).powi(2) ;
    let norm_ket = cv_spaces::norm(&curr_ket).powi(2) ;
    let prob = mod_at_index/norm_ket ;

    prob }


pub fn trans_prob2(curr_ket: &Vec<Comp> , target_ket: &Vec<Comp>) ->Comp{
    let norm1 = cv_spaces::norm(&curr_ket) ;
    let norm2 = cv_spaces::norm(&target_ket) ;
    let norm_prod = norm1 * norm2 ;
    
// note c_vec_dot actually conjugates first to produce the hermittian product 
    let trans_amp = cv_spaces::c_vec_dot(&target_ket, &curr_ket) ;
    let new_comp = complex::new( cv_spaces::normalize(trans_amp.real /norm_prod) ,cv_spaces::normalize( trans_amp.imag/norm_prod ) ); 

    new_comp
  // target_conj
}
// only works for hermitian matrix 2 by 2  
pub fn eigen_values0(omega: &Vec<Vec<Comp>>) ->Vec<f64>{
    assert!(cv_spaces::herm_test(&omega) , "Must be Hermitian") ;
    let a = omega[0][0].clone() ;
    let b = omega[0][1].clone() ;
    let c = omega[1][0].clone() ;
    let d = omega[1][1].clone();
    let d_inv = complex::inv(&d) ;
    let a_min_d = complex::add(&a , &d_inv) ;
    let four = complex::new( 4f64 , 0f64) ;
    let fourb = complex::mul(&four , &b) ;
    let fourbc = complex::mul(&fourb , &c) ;
    let a_min_dsq = complex::mul(&a_min_d , &a_min_d) ;
    let det = complex::add(&a_min_dsq , &fourbc).real.sqrt() ;
    // for a hermittian matrix we are guaranteed that the eigen values will be real and we know the
    // diagonal elemnts will be zero so we can now focu on the real aspects 
    
    let lambda_1 = ((a.real + d.real) + det )/2f64  ;
     let lambda_2 = ((a.real + d.real) - det )/2f64  ;

     vec![lambda_1 , lambda_2 ]
     
}

pub fn eigen_vectors0(omega: &Vec<Vec<Comp>>) ->Vec<Vec<Comp>>{

    let eigen_vals = eigen_values0(&omega) ; 
    let lamb_comp1 = complex::new(eigen_vals[0].clone() , 0f64 ) ;
    let lamb_comp2 = complex::new(eigen_vals[1].clone() , 0f64 ) ;

    let a_inv = complex::inv(&omega[0][0]) ;
     let d_inv = complex::inv(&omega[1][1]) ;

 let mut final_: Vec<Vec<Comp>> = Vec::new() ;
 if omega[0][1].real != 0.0 {
  let eigen_vec1 = vec![omega[0][1].clone() ,complex::add(&lamb_comp1 , &a_inv)];
  let eigen_vec2 = vec![omega[0][1].clone() ,complex::add(&lamb_comp2 , &a_inv)];

  final_ = vec![eigen_vec1.clone() , eigen_vec2.clone() ];
 }
else if omega[1][0].real != 0.0 {
  let eigen_vec2 = vec![complex::add(&lamb_comp1 , &d_inv),omega[1][0].clone()];
  let eigen_vec1 = vec![complex::add(&lamb_comp2 , &d_inv),omega[1][0].clone()];


    final_ =  vec![eigen_vec1.clone() , eigen_vec2.clone() ];

}

else {
  let eigen_vec2 = vec![complex::new(1f64 , 0f64),complex::new(0f64 , 0f64)];
  let eigen_vec1 = vec![complex::new(0f64 , 0f64),complex::new(1f64 , 0f64)];



final_ =  vec![eigen_vec1.clone() , eigen_vec2.clone() ];

}

final_
}

pub fn observ_on_state(observable: &Vec<Vec<Comp>> , curr_state: &Vec<Comp>) ->States_N_Prob {

    let eigen_vecs: Vec<Vec<Comp>> = eigen_vectors0(&observable);
     let eigen_vals = eigen_values0(&observable) ; 
    let mut probs: Vec<Comp> = Vec::new() ;
    let mut mean: f64 = 0f64 ;
    for x in 0..eigen_vecs.len() {
        let trans_amp = trans_prob2(&curr_state , &eigen_vecs[x]) ;
        let modsq = complex::mod_(&trans_amp).powi(2) ;
        mean = mean + modsq*eigen_vals[x] ;
        probs.push(trans_amp)
}

States_N_Prob{
    mean: mean ,
    omega: observable.clone() ,
    fr_state: curr_state.clone() ,
    to_states: eigen_vecs.clone() ,
    probs: probs.clone() }
    }




