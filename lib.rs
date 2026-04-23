
use lagrangian_interpolation::lagrange_interpolate ;


pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

mod basic ;
mod poly2 ;
mod nthfib ;
mod merk ;

fn main () {
    println!("Seems we compiled well");
}






#[cfg(test)]
mod tests {
    use crate::basic::Omega ;
    use crate::basic ;
    use crate::basic::powers ;
    use crate::poly2 ;
    use crate::nthfib ;
    use crate::merk ;
    use num_bigint::BigUint ;
    use num_traits::cast::ToPrimitive ;
    
 
  
    #[test]
    fn it_works() {
     let fib_targ: usize = 9  ;// can go to 92 , but the parameters for now are for figures less
                              // than 100 k so we max out at fib 20 ;
     let fibs: Vec<i64> = nthfib::nthfib(fib_targ) ;
     

     // the function below takes a number n , usize , which is length of a collection and on its
     // basis decides what original subgroup order k to use in original lagrage interpolation 
     fn dec_k (n: u32) -> u32 {
         if n <= 1 { return 0 ; }
         // subtract 1 so that powers of 2( like 16 ) dont jump to  the next level 
         let adjusted_n  =  n - 1 ; 

        // calculate total bits minus leading zeros to get the position of the highest bit 
         32  - adjusted_n.leading_zeros()
     }

     let  k: u32 =(2_i32.pow( dec_k(fibs.len() as u32))) as u32 ;
     println!("k : {}" , k) ;
     let shift_par: i64 = 8 ;// later make sure there is a random generator here for the shift parameter 
     let fin_cont: Vec<Omega> = basic::prod_sub_gr(k , 0) ;
     let fin_proce: Vec<Omega>  =(fin_cont.iter()).map(move |x| x.process_bigger_than_kov2( k)).collect();
     let trace: Vec<i64> = fibs ;//vec![1,3,4,7,11,18 , 29 , 47];
      let indeces: Vec<usize> = (0..=(trace.len() - 1)).collect();
      let mut points: Vec<(i64 , i64)> = Vec::new();
      for x in indeces {points.push((basic::eval_Omega(&fin_proce[x] , k ) as i64 ,trace[x]  ))};
     
    let test: Vec<(i64 , i64)> = vec![(1 , 6), (2 , 15) , (3 , 28)]; 
     let polyn: Vec<i64> =  poly2::lagrange_interp(&points); 

     let factor: u32 = 4 ;
     let new_order: i64 = (factor * k ) as i64 ;
     let disg_trace_xs: Vec<Omega> = basic::prod_sub_gr((factor*k) , 0) ;
     //println!("Disguised Trace xs : {:?}" , disg_trace_xs ) ;
    
     let disg_proc: Vec<Omega> = (disg_trace_xs.iter()).map(move |x| x.process_bigger_than_kov2(factor * k)).collect();

     let mut disg_xs: Vec<i64> = Vec::new() ;

     //println!("Disguised processed Trace xs : {:?}" , disg_proc ) ;
     let indeces2: Vec<usize> = (0..=(disg_proc.len() - 1)).collect();
     let indeces3: Vec<usize> = indeces2.clone();
     let indeces4: Vec<usize> = indeces2.clone() ;
     let indeces5: Vec<usize> = indeces2.clone();
     let indeces6: Vec<usize> =  indeces2.clone();
     let mut values: Vec<i64> = Vec::new() ;
     for x in indeces2 {
         disg_xs.push(poly2::mod_mul(shift_par, poly2::evaluate_Omega(&disg_proc[x] ,new_order ).into() ));
         values.push(poly2::evaluate_poly(&polyn ,poly2::mod_mul(shift_par, poly2::evaluate_Omega(&disg_proc[x] ,new_order ).into() )))}

    // constraint polynomial , this is what has to change in adapting to a new trace 
   let disg_xs2: Vec<i64> = disg_xs.clone();
   let disg_xs3: Vec<i64> = disg_xs.clone();
   let disg_xs4: Vec<i64> = disg_xs.clone() ;
   let values_0_2 = values.clone() ;
   let values_0_3 = values.clone() ;
   let values_0_4 = values.clone();
   let root1: [u8 ; 32 ] =  merk::root_prod(&values) ;
   let alpha0_big  =  BigUint::from_bytes_le(&root1)     ;
   let p = BigUint::from(97u32);
   let alpha0 = &alpha0_big % &p ;//.to_i64().unwrap() ;
   
   
   println!("alpha_0 :{}" , alpha0) ;
   let mut values2: Vec<i64> = Vec::new() ;
   let mut values3: Vec<i64> = Vec::new() ;
   let mut values4: Vec<i64> = Vec::new() ;
   let mut values5: Vec<i64> = Vec::new() ;

   for x in indeces3 {
       values2.push( poly2::mod_mul(poly2::mod_sub(poly2::mod_mul(ToPrimitive::to_i64(&alpha0).unwrap() ,values_0_2[x]),1 ),poly2::mod_inv(poly2::mod_sub(disg_xs2[x], 1)))) }

 

   let root2: [u8 ; 32 ] = merk::root_prod(&values2);
   let alpha1_big = BigUint::from_bytes_le(&root2);
   let alpha1 = &alpha1_big % &p ;
 for x in indeces4 {
     let last_ys = values_0_3.last().unwrap().clone() ;
     let last_xs = disg_xs3.last().unwrap().clone() ;
      values3.push( poly2::mod_mul(poly2::mod_sub(poly2::mod_mul(ToPrimitive::to_i64(&alpha1).unwrap() ,values_0_3[x]),last_ys ),poly2::mod_inv(poly2::mod_sub(disg_xs3[x], last_xs)))) }

let root3: [u8 ; 32 ] = merk::root_prod(&values3);
   let alpha2_big = BigUint::from_bytes_le(&root3);
   let alpha2 = &alpha2_big % &p ;
 for x in indeces5 {
     let last_xs = disg_xs4.last().unwrap().clone();
      values4.push(poly2::mod_mul(poly2::mod_mul(poly2::mod_sub(poly2::mod_sub(values_0_4[((x + 2) as i64 % new_order) as usize] , values_0_4[((x + 1) as i64 % new_order) as usize] ) ,values_0_4[((x + 0) as i64 % new_order) as usize]),poly2::mod_sub(disg_xs4[x], last_xs)),poly2::mod_inv(poly2::mod_sub(poly2::mod_pow(disg_xs4[x], new_order), 1))))
          
  
 }

for x in indeces6 {
       
       values5.push(poly2::mod_add(poly2::mod_add(values2[x] , values3[x]) ,values4[x ])) }

 let root4: [u8 ; 32 ] = merk::root_prod(&values4);
   let index_big = BigUint::from_bytes_le(&root4);
   let for_np = new_order.clone();
   let np = BigUint::from(for_np as u32);

   let index = &alpha2_big % &np ;


fn con_poly_eval(xs: Vec<i64> , ys: Vec<i64>)-> (){
}
    println!("values are:{:?}" ,  values5.len());

     

    }
}
