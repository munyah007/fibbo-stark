
use lagrangian_interpolation::lagrange_interpolate ;
use rs_merkle::{MerkleProof , algorithms::Sha256} ;


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
    use rs_merkle::{MerkleProof , algorithms::Sha256} ;
    const MOD: i64 = 97 ;
    const MOD_2: i64 = 97 ;

    

 
  
    #[test]
    fn it_works() {

  fn dec_k (n: u32) -> u32 {
         if n <= 1 { return 0 ; }
         // subtract 1 so that powers of 2( like 16 ) dont jump to  the next level 
         let adjusted_n  =  n - 1 ; 

        // calculate total bits minus leading zeros to get the position of the highest bit 
         32  - adjusted_n.leading_zeros()
     }

     let fib_targ: usize = 9  ;// can go to 92 , but the parameters for now are for figures less
                              // than 100 k so we max out at fib 20 ;
                              //
      let  k: u32 =(2_i32.pow( dec_k(fib_targ as u32))) as u32 ;

     
     let fibs: Vec<i64> = nthfib::nthfib(k as usize) ;
     println!("Target {} , Length is{}" , fib_targ , fibs.len()) ;
     let out_fp: i64 = fibs.last().unwrap().clone();
     let out: i64 = fibs.last().unwrap().clone();
     let psn_out: usize = fibs.len() - 1 ;
     let root_rt: [u8 ; 32 ] =  merk::root_prod(&fibs) ;
     println!("Out position:{}", psn_out) ;

     // the function below takes a number n , usize , which is length of a collection and on its
     // basis decides what original subgroup order k to use in original lagrage interpolation 
     

     
     println!("k : {}" , k) ;
     let shift_par: i64 = 4 ;// later make sure there is a random generator here for the shift parameter 
     let fin_cont: Vec<Omega> = basic::prod_sub_gr(k , 0) ;
     let fin_proce: Vec<Omega>  =(fin_cont.iter()).map(move |x| x.process_bigger_than_kov2( k)).collect();
     let trace: Vec<i64> = fibs ;//vec![1,3,4,7,11,18 , 29 , 47];
    
      println!("fin_proce : {:?}" , fin_proce) ; 

     
     let x_at_out: i64 = poly2::evaluate_Omega(&fin_proce[psn_out] , k as i64) as i64 ;
     let x_pen_out: i64 = poly2::evaluate_Omega(&fin_proce[psn_out - 1], k as i64) as i64 ;
     println!("First , x at out {}", x_at_out) ;
      let indeces: Vec<usize> = (0..=(trace.len() - 1)).collect();
      let mut points: Vec<(i64 , i64)> = Vec::new();
      for x in indeces {points.push((poly2::evaluate_Omega(&fin_proce[x] , k as i64 ) as i64 ,trace[x]  ))};
     
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
      values3.push( poly2::mod_mul(poly2::mod_sub(poly2::mod_mul(ToPrimitive::to_i64(&alpha1).unwrap() ,values_0_3[x]),out_fp ),poly2::mod_inv(poly2::mod_sub(disg_xs3[x], x_at_out)))) }

let root3: [u8 ; 32 ] = merk::root_prod(&values3);
   let alpha2_big = BigUint::from_bytes_le(&root3);
   let alpha2 = (&alpha2_big % &p) ;
   let alpha2_ = ToPrimitive::to_i64(&alpha2).unwrap() ;
 for x in indeces5 {
    
     let fx_p2 = values_0_4[((x + 2) as i64 % new_order ) as usize ] ;
     let fx_p1 = values_0_4[((x + 1) as i64 % new_order ) as usize] ;
     let fx_ = values_0_4[(x as i64 % new_order ) as usize ] ;
     let num_fac1 = poly2::mod_sub_2(&disg_xs4[x] ,&x_at_out ) ;
     let num_fac2 = poly2::mod_sub_2(&disg_xs4[x] ,&x_pen_out);
     let x_G = poly2::mod_pow(disg_xs4[x], k as i64 ) ;
     let con: i64 = 1 ;
     let for_inv = poly2::mod_sub_2(&x_G ,&con) ;  
     let den = poly2::mod_inv(for_inv) ;
     let sub_t1 = poly2::mod_sub_2(&fx_p1 , &fx_);
     let num_fac_main = poly2::mod_sub_2(&fx_p2 , &sub_t1 );
     let num_fac_main2 = poly2::mod_mul_2(&num_fac1 , &num_fac2);
     let num = poly2::mod_mul_2(&num_fac_main ,&num_fac_main2) ;
     let fin = poly2::mod_mul_2(&num ,&den);
     let fin_ = poly2::mod_mul_2(&fin ,&alpha2_) ;
        values4.push(fin_) ;  
  
 }

for x in indeces6 {
       
       values5.push(poly2::mod_add(poly2::mod_add(values2[x] , values3[x]) ,values4[x ])) }

let root4: [u8 ; 32 ] = merk::root_prod(&values4);

 let root5: [u8 ; 32 ] = merk::root_prod(&values5);
   let index_big = BigUint::from_bytes_le(&root5);
   let for_np = new_order.clone();
   let np = BigUint::from(for_np as u32);

   let index = (&index_big % &np).to_usize().expect("Index too large")   ;
   let pindex_p1 = ((index + 1) as i64 % new_order ) as usize ;
   let pindex_p2 = ((index + 2) as i64 % new_order ) as usize ;
 

//#[derive(Debug)]
struct Proof {
    index: usize ,
    root_rt: [u8 ; 32 ],// root of real trace 
    root_fx: [u8 ; 32 ] , // root1
    root_p0x: [u8 ; 32 ],
    root_p1x: [u8 ; 32 ],
    root_p2x: [u8 ; 32 ],
    root_cx: [u8 ;32 ],
    proof_fx: MerkleProof<Sha256>,
    proof_fx_p1: MerkleProof<Sha256> ,
    proof_fx_p2: MerkleProof<Sha256> ,
    proof_p0x: MerkleProof<Sha256>,
    proof_p1x: MerkleProof<Sha256>,
    proof_p2x: MerkleProof<Sha256>,
    proof_cx: MerkleProof<Sha256>,
    field: i64 ,
    eval_order: i64 ,
    alpha0_p0x: i64 ,
    alpha1_p1x: i64 ,
    alpha2_p2x: i64 ,
    cx: i64 ,
    fx: i64 ,
    fx_p1: i64 ,
    fx_p2: i64 ,
    out: i64 ,
    shift: i64 ,
    fib_targ: usize ,
    orig_k: u32 ,

};
 let alpha_p0x: i64 = values2[index].clone() ;
 let alpha_p1x: i64 = values3[index].clone() ;
 let alpha_p2x: i64 = values4[index].clone() ;
 let cx: i64 = values5[index].clone();
 let fx: i64 = values[index].clone();
 let fx_p1: i64 = values[pindex_p1].clone() ;
 let fx_p2: i64 = values[pindex_p2].clone() ;
 let cx_cl: Vec<i64> = values5.clone() ;
let proof: Proof = Proof{
                   proof_fx: merk::proof_prod(&values , index),
                   proof_fx_p1: merk::proof_prod(&values , pindex_p1) ,
                   proof_fx_p2: merk::proof_prod(&values , pindex_p2),
                   proof_p0x: merk::proof_prod(&values2 ,index),
                   proof_p1x: merk::proof_prod(&values3 , index),
                   proof_p2x: merk::proof_prod(&values4 , index) ,
                   proof_cx: merk::proof_prod(&values5 , index) ,
                   index: index.clone() ,
                   root_fx: root1.clone() ,
                   root_p0x: root2.clone() ,
                   root_p1x: root3.clone() ,
                   root_p2x: root4.clone() ,
                   root_cx: root5.clone() ,
                   field: 97 ,
                   eval_order: new_order.clone() ,
                   alpha0_p0x: alpha_p0x,
                   alpha1_p1x: alpha_p1x,
                   alpha2_p2x: alpha_p2x ,
                   cx: cx ,
                   fx: fx, 
                   fx_p1: fx_p1 ,
                   fx_p2: fx_p2 ,
                   out: out ,
                   root_rt: root_rt,
                   shift: shift_par ,
                   fib_targ: fib_targ, 
                   orig_k: k };

fn reduce_layer(vals:  &Vec<i64>) -> Vec<i64> {
let order: usize = vals.len() ;
let half: usize = (order / 2)   ;
let half_ = half.saturating_sub(1);
let indeces: Vec<usize> = (0..=half_).collect();
let mut container: Vec<i64> = Vec::new() ;
 let xs: Vec<Omega> = basic::prod_sub_gr((order as u32) , 0) ;
  let xs_proc: Vec<Omega> = (xs.iter()).map(move |x| x.process_bigger_than_kov2(order as u32)).collect();
   let root: [u8 ; 32 ] = merk::root_prod(&vals);
   let alpha_big = BigUint::from_bytes_le(&root);
    

   let  p_div  = BigUint::from(MOD_2 as u32 );

   let alpha = (&alpha_big % &p_div) ;
   let alpha_ = ToPrimitive::to_i64(&alpha).unwrap(); 
for x in indeces {
    let x_at_ind = poly2::evaluate_Omega(&xs_proc[x] ,order as i64);
    let x_ = (x + half );
    let con = 2 ;
    let term1_num = poly2::mod_add_2(&vals[x] , &vals[x_]) ;
    let term1_den = poly2::mod_inv(2);
    let term2_num = poly2::mod_sub_2(&vals[x] ,&vals[x_]) ;
    let term2_den_pen = poly2::mod_mul_2(&con ,&x_at_ind) ;
    let term2_den = poly2::mod_inv(term2_den_pen);
    let term1 = poly2::mod_mul_2(&term1_num ,&term1_den);
    let term2 = poly2::mod_mul_2(&term2_num ,&term2_den);
    let term2_ = poly2::mod_mul_2(&term2 , &alpha_); 
    let term = poly2::mod_add_2(&term1 , &term2_ );
    container.push(term) ;
}
container 
}

struct LayerProof {
     index: usize ,
     root_vals: [u8 ; 32 ],// root of real trace 
     root_nl: [u8 ; 32 ] , // root1
     proof_fx: MerkleProof<Sha256>,
     proof_f_min_x: MerkleProof<Sha256>,
     proof_nl: MerkleProof<Sha256> ,
     fx: i64 ,
     f_min_x: i64 ,
     fx_nl: i64 ,
     length: usize ,

 }

 fn layer_a_to_b_proof (vals: &Vec<i64> ) -> LayerProof{
   let root: [u8 ; 32 ] = merk::root_prod(&vals);

   let index_big = BigUint::from_bytes_le(&root);
   let for_np = (vals.len() as i64 / 2 ) as usize ;// index should come from the lower half 
   let half = for_np.clone() ;
   let np = BigUint::from(for_np as u32);
   let index = (&index_big % &np).to_usize().expect("Index too large")   ;
   let next_layer: Vec<i64> = reduce_layer(&vals ) ;
   println!("New layer {:?}" ,next_layer);
   let fx = vals[index].clone() ;
   let f_min_x = vals[index + for_np] ;
   
   let red_val = next_layer[index] ;// the value of interest will be at the same index in its layer 
   let next_layer_root = merk::root_prod(&next_layer);
   
     LayerProof{
         length: vals.len() ,
         index: index , // has to be verified that it came from the root 
         root_vals: root ,
         root_nl: next_layer_root ,
         proof_fx:   merk::proof_prod(&vals , index)  ,
         proof_f_min_x:   merk::proof_prod(&vals , index + for_np)  ,
         proof_nl: merk::proof_prod(&next_layer , index)  , // used to prove that the commited value is indeed in the root 
         fx: fx ,
         f_min_x: f_min_x ,
         fx_nl: red_val ,

     }


 }

let layer2: Vec<i64> = reduce_layer(&cx_cl) ;
println!("Layer 2 {:?}" , layer2) ;

fn produce_comp_proofs(input: &Vec<i64>) -> Vec<LayerProof> {
let mut complete_proofs : Vec<LayerProof> = Vec::new() ;
let mut curr_input: Vec<i64> = reduce_layer(input) ;
complete_proofs.push(layer_a_to_b_proof(input));
let mut counter: usize = 2 ;
while curr_input.len() != 1  {
    complete_proofs.push(layer_a_to_b_proof(&curr_input)) ;
    curr_input = reduce_layer(&curr_input) ;
    println!(" Round {} done" , counter) ;
    counter += 1 ;
 }
println!("All proofs produced") ;
complete_proofs

}

let mut  comp_proof: Vec<LayerProof> = produce_comp_proofs(&cx_cl) ;
let max_index: usize  = produce_comp_proofs(&cx_cl).len() - 1 ;

fn verify_individual_proof(layer_proof: LayerProof ) -> bool {
//verify index came from the root 
 let root: [u8 ; 32 ] = layer_proof.root_vals;
 let root_1: [u8 ; 32 ] = root.clone() ;
   let index_big = BigUint::from_bytes_le(&root);
   let for_np = (layer_proof.length as i64 / 2 ) as usize ;// index should come from the lower half 
   let half = for_np.clone() ;
   let np = BigUint::from(for_np as u32);
   let index = (&index_big % &np).to_usize().expect("Index too large")   ;
   assert_eq!( index , layer_proof.index) ;
   println!("Index verified , that it came from root !") ;
  let fx_ver = merk::verify_proof(layer_proof.proof_fx ,layer_proof.root_vals, layer_proof.fx , layer_proof.index ,layer_proof.length) ;
    assert!(fx_ver);
    println!("fx verified to be at index in the root ") ;
 let f_minx_ver = merk::verify_proof(layer_proof.proof_f_min_x ,root_1 , layer_proof.f_min_x , layer_proof.index + half  ,layer_proof.length) ;
    assert!(f_minx_ver);
    println!("f_minx verified to be at index in the root ") ; 
 let f_x_nl_ver = merk::verify_proof(layer_proof.proof_nl ,layer_proof.root_nl, layer_proof.fx_nl , layer_proof.index  ,half) ;
    assert!(f_x_nl_ver);
    println!("f_x_nl_ver  verified to be at index in the root ") ;

    // calculating relevant x for verification in calculation
   let eval_xs: Vec<Omega> = basic::prod_sub_gr(layer_proof.length  as u32, 0) ;
  

  let orig_proc_out: Omega = eval_xs[index].process_bigger_than_kov2(layer_proof.length as u32 );
  
  let x_at_index = poly2::evaluate_Omega(&orig_proc_out ,layer_proof.length as i64) ;
  println!("x successfully calculated !!") ;

  // performing the calculation
  //
 let alpha_big = BigUint::from_bytes_le(&root_1);
    

   let  p_div  = BigUint::from(MOD_2 as u32 );

   let alpha = (&alpha_big % &p_div) ;
   let alpha_ = ToPrimitive::to_i64(&alpha).unwrap(); 
  let con = 2 ;
    let term1_num = poly2::mod_add_2(&layer_proof.fx , &layer_proof.f_min_x) ;
    let term1_den = poly2::mod_inv(2);
    let term2_num = poly2::mod_sub_2(&layer_proof.fx ,&layer_proof.f_min_x) ;
    let term2_den_pen = poly2::mod_mul_2(&con ,&x_at_index) ;
    let term2_den = poly2::mod_inv(term2_den_pen);
    let term1 = poly2::mod_mul_2(&term1_num ,&term1_den);
    let term2 = poly2::mod_mul_2(&term2_num ,&term2_den);
    let term2_ = poly2::mod_mul_2(&term2 , &alpha_); 
    let term = poly2::mod_add_2(&term1 , &term2_ );
    assert_eq!(term , layer_proof.fx_nl) ;
    println!("Reduction was valid!!") ;
   true 
}


fn verify_all_layer_proofs(mut lay_proofs:  Vec<LayerProof>, max_index: usize)-> bool {
let vals: Vec<i64> = vec![0,0] ; // dummy values for initialisation 
let mut  root_next_for_last: [u8 ; 32 ] = merk::root_prod(&vals);
for x in (0..= max_index) {
    let proof = lay_proofs.remove(0) ;
    if x == 0 {
        root_next_for_last = proof.root_nl.clone() 
         }
    else  {
        // asserting that what was produced as next layer and hence commited in the last layer is
        // what is being used at the vals for this layer thus ensuring linkage of the provided
        // layers .That the values used comes from those layers is what is accomplished by
        // verify_individual_proof .
        assert_eq!(root_next_for_last , proof.root_vals) ;
        root_next_for_last  = proof.root_nl.clone() ;
    }
    let verdict  = verify_individual_proof(proof) ;
    
    assert!(verdict) ;
    println!("We are at {}" , x) ;
}
println!("all proofs individually verified and linkage established .") ;

true 
}

let allproofs: bool = verify_all_layer_proofs(comp_proof, max_index) ;
println!("All layers proofs :{}", allproofs) ;

  fn verify (input: Proof) ->() {
  let eval_xs: Vec<Omega> = basic::prod_sub_gr((input.eval_order as u32) , 0) ;
  
  let orig_xs: Vec<Omega> = basic::prod_sub_gr(input.orig_k , 0 );

  let orig_proc_out: Omega = orig_xs[ (input.orig_k - 1) as usize].process_bigger_than_kov2(input.orig_k );
  let orig_proc_pen: Omega = orig_xs[(input.orig_k - 2) as usize].process_bigger_than_kov2(input.orig_k);
  let x_orig_at_out = poly2::evaluate_Omega(&orig_proc_out ,input.orig_k as i64) ;
  let x_orig_pen_out = poly2::evaluate_Omega(&orig_proc_pen , input.orig_k as i64);
  //println!("Disguised Trace xs : {:?}" , disg_trace_xs ) ;
    
     let eval_proc: Omega = eval_xs[input.index].process_bigger_than_kov2(input.eval_order as u32 );

    

      let rel_x = poly2::mod_mul(input.shift, poly2::evaluate_Omega(&eval_proc,input.eval_order ).into() );
      
        


    // stages in verification
    // clone roots 
    let root_fx = input.root_fx.clone() ;
    let root_fx_forp1 = input.root_fx.clone() ;
    let root_fx_forp2 = input.root_fx.clone() ;

    let root_p0x = input.root_p0x.clone() ;
    let root_p1x = input.root_p1x.clone() ;
    let fx = input.fx.clone() ;

    // Phase 1 : Correctness of proof
      //1. Test that the index came from the merkle root cx 

    let index_big = BigUint::from_bytes_le(&input.root_cx);
   let for_np = input.eval_order;
   let np = BigUint::from(for_np as u32);

   let index = (&index_big % &np).to_usize().expect("Index too large"); 

   assert_eq!(index ,input.index , "Indeces dont match") ;
   println!("Index passed") ;
   let index_vec: Vec<usize> = vec![index] ;
    
   //2.test that out is at the last position of root_rt (this may be included in the constraints
   //  so leave for now )
   //3.test that the value of fx given is at index in root_fx 
    
    let fx_ver = merk::verify_proof(input.proof_fx ,input.root_fx , input.fx , input.index ,input.eval_order as usize) ;
    assert!(fx_ver);
    println!("fx passed") ; 
   //4. calculate alph0 from root1
   let alpha0_big  =  BigUint::from_bytes_le(&root_fx)     ;
   let p = BigUint::from(97u32);
   let alpha0 = &alpha0_big % &p ;//.to_i64().unwrap() ;
  

   //5. test that the commited alpha0_p0x is in root_p0x
 let alpha0_ver = merk::verify_proof(input.proof_p0x ,input.root_p0x , input.alpha0_p0x , input.index ,input.eval_order as usize) ;
    assert!(alpha0_ver);
    println!("alpha0_p0x  passed") ; 


   //NB :no need to test equality for a value whose commitment is not first established 
   //6. calculate alpha0 *p0x which uses fx , and assert equality with tested value
  let con: i64 = 1 ;
  let den  = poly2::mod_inv(poly2::mod_sub_2(&rel_x, &con)) ;
  //by using locally produced alpha0 in calculation we are establishing correctness
  let mt_cant = ToPrimitive::to_i64(&alpha0).unwrap(); 
  let sub1  = poly2::mod_mul_2(&mt_cant ,&input.fx);
  let calc_alphap0x_pen =   poly2::mod_sub_2(&sub1,&con ); 
  let calc_alphap0x_fin = poly2::mod_mul_2(&calc_alphap0x_pen , &den ) ;


   assert_eq!(input.alpha0_p0x , calc_alphap0x_fin) ;
   println!("P0X passed ") ;
   //7. Repeat 5, and 6 for alpha1_p1x and the version wit 2
   //
     let alpha1_big  =  BigUint::from_bytes_le(&root_p0x)     ;

   let alpha1 = &alpha1_big % &p ;//.to_i64().unwrap() ;
  

   //5. test that the commited alpha1_p1x is in root_p1x
 let alpha1_ver = merk::verify_proof(input.proof_p1x ,input.root_p1x , input.alpha1_p1x, input.index ,input.eval_order as usize) ;
    assert!(alpha1_ver);
    println!("alpha1_p1x  passed") ; 
  ;
  
  let con: i64 = input.out ;

  let orig_xs: Vec<Omega> = basic::prod_sub_gr((input.orig_k) , 0) ;
     //println!("Disguised Trace xs : {:?}" , disg_trace_xs ) ;
    
     let orig_proc: Omega = orig_xs[(input.orig_k - 1) as usize ].process_bigger_than_kov2(input.orig_k );

    

      let x_at_out: i64  = poly2::evaluate_Omega(&orig_proc,input.orig_k as i64 ).into();
       println!("Second x at out: {}" , x_at_out);
  let den_p1  = poly2::mod_inv(poly2::mod_sub_2(&rel_x, &x_at_out)) ;
  // by using the locally produced alpha1 in calculation we establishing correctness 
  let mt_cant_p1 = ToPrimitive::to_i64(&alpha1).unwrap(); 
  let sub1_p1  = poly2::mod_mul_2(&mt_cant_p1 ,&input.fx);
  // below is the constraint that establishes that the provided out is actually in the trace at the
  // required index as per original trace
  let calc_alphap1x_pen =   poly2::mod_sub_2(&sub1_p1,&input.out); 
  let calc_alphap1x_fin = poly2::mod_mul_2(&calc_alphap1x_pen , &den_p1 ) ;
   assert_eq!(input.alpha1_p1x , calc_alphap1x_fin) ;
   println!("P1X passed ") ;
   // alpha2_p2x tests
     let pindex_p1 = ( (input.index + 1) as i64 % input.eval_order ) as usize ;
      let pindex_p2 = ( (input.index + 2) as i64 % input.eval_order ) as usize ;

 let alpha2_ver = merk::verify_proof(input.proof_p2x ,input.root_p2x , input.alpha2_p2x, index ,input.eval_order as usize) ;
    assert!(alpha2_ver);
    println!("alpha1_p2x commitment  passed") ; 

  

     let fx_p2 = input.fx_p2.clone() ;
     let fx_p1 = input.fx_p1.clone() ;
      // fx was destructured earlier
      let fx_p1_ver = merk::verify_proof(input.proof_fx_p1 ,input.root_fx , input.fx_p1 , pindex_p1 ,input.eval_order as usize) ;
    assert!(fx_p1_ver);
    println!("fx_p1  passed") ; 
  let fx_p2_ver = merk::verify_proof(input.proof_fx_p2 ,input.root_fx , input.fx_p2 , pindex_p2 ,input.eval_order as usize) ;
    assert!(fx_p1_ver);
    println!("fx_p2 passed") ; 

      //
      let alpha2_big  =  BigUint::from_bytes_le(&root_p1x)     ;

      let alpha2 = &alpha2_big % &p ;//.to_i64().unwrap() ;
       let alpha2_ = ToPrimitive::to_i64(&alpha2).unwrap() ;

       let num_fac1 = poly2::mod_sub_2(&rel_x ,&x_orig_at_out ) ;
     let num_fac2 = poly2::mod_sub_2(&rel_x ,&x_orig_pen_out);

     
     let x_G = poly2::mod_pow(rel_x, input.orig_k as i64 ) ;
     let con: i64 = 1 ;
     let den = poly2::mod_inv(poly2::mod_sub_2(&x_G ,&con)) ;
     let sub_t1 = poly2::mod_sub_2(&fx_p1 , &fx);
     let num_fac_main = poly2::mod_sub_2(&fx_p2 , &sub_t1 );
     let num_fac_main2 = poly2::mod_mul_2(&num_fac1 , &num_fac2);
     let num = poly2::mod_mul_2(&num_fac_main ,&num_fac_main2) ;
     let fin = poly2::mod_mul_2(&num ,&den);
     let fin_ = poly2::mod_mul_2(&fin ,&alpha2_) ;
      

       assert_eq!(input.alpha2_p2x , fin_) ;
       println!("P2X passed ") ;



   //8. test the commited value of cx - is it in the root at the index
  let cx_ver = merk::verify_proof(input.proof_cx,input.root_cx , input.cx ,index ,input.eval_order as usize) ;
    assert!(cx_ver);
    println!(" commited cx  passed") ; 

    let cx_eval_1 = poly2::mod_add_2(&fin_ ,&calc_alphap0x_fin) ;
    let cx_eval_fin = poly2::mod_add_2(&cx_eval_1 , &calc_alphap1x_fin);
     assert_eq!( cx_eval_fin , input.cx) ;
     println!("cx value accurate , accuracy phase successfully verified");
   //9. add the values verified in 5 and 6 to get cx and assert equality
   //10 .Verification phase 1 passes and we go to phase 2 
   //Phase 2 : Low degreeness of the commited cx polynomial 

  }
verify(proof) ;
 

    }
}
