const MOD: i64 = 97; //99991 ;
use crate::basic::Omega ;
pub fn mod_add (a: i64 , b: i64) -> i64{
    (a + b)%MOD } 
pub fn mod_add_2 (a: &i64 , b: &i64) -> i64{
    (a + b)%MOD } 

pub fn mod_sub (a: i64 , b: i64)-> i64 {
    (a - b).rem_euclid( MOD) }

pub fn mod_sub_2 (a: &i64 , b: &i64)-> i64 {
    (a - b).rem_euclid( MOD) }


pub fn mod_mul (a: i64 , b: i64) -> i64 {
    (a * b ) % MOD 
}

pub fn mod_mul_2 (a: &i64 , b: &i64) -> i64 {
    (a * b ) % MOD 
}


//fermat inverse 
pub fn mod_pow (mut base: i64 , mut exp: i64) -> i64 {
    let mut result = 1 ;
    base %= MOD ;

    while exp > 0 {
        if exp % 2 == 1 {
            result = mod_mul(result , base) ;
        }
        base = mod_mul(base , base ) ;

        exp /= 2 ;
    }
    result 
}

pub fn mod_inv (a: i64 ) -> i64{
    mod_pow( a, MOD - 2)

}



pub fn evaluate_Omega( oms: &Omega, order: i64) -> i64 {
    let mut  coef: i64 = oms.coef as i64  ;
    if coef < 0 {coef =  mod_add( coef.into() , MOD ) }
 mod_mul(coef as  i64 , mod_pow(order , oms.power as i64)) }
//polynomial helpers 
//
fn poly_mul (a: &Vec<i64> , b: &Vec<i64>) -> Vec<i64> {
    let mut res = vec![0 ; a.len() + b.len() - 1] ;
    for i in 0..a.len(){
        for j in 0..b.len() {
        res[i + j] = mod_add(res[i + j], mod_mul(a[i] , b[j]));}
    }
    res 
}

fn poly_add(a: &Vec<i64> , b: &Vec<i64>) -> Vec<i64> {
    let n = a.len().max(b.len());
    let mut res = vec![0;n] ;
    for i in 0..n {
        let av = if i < a.len(){ a[i] } else { 0};
        let bv = if i < b.len() { b[i] } else {0};
        res[i] = mod_add(av , bv );
    }

    res 
}

//scalar mult 
fn poly_scalar_mul(poly: &Vec<i64> , scalar: i64) -> Vec<i64>{
    poly.iter().map(|&c| mod_mul(c ,scalar)).collect()}

    // lagrange interpolation 
pub fn lagrange_interp (points: &Vec<(i64 ,i64)>) -> Vec<i64> {
    let mut result = vec![0] ;
    let tracks1: Vec<i64> = (0..=(points.len() - 1) as i64).collect();
     let tracks2: Vec<i64> = (0..= (points.len() - 1) as i64).collect();
     for i in tracks1 {
        let (xi , yi) = points[i as usize];
       //start with poly n = 1 
        let mut basis = vec![1] ;
        let mut denom = 1 ;
        for j in &tracks2 {
            if i == *j {
                continue ;  } 
            let (xj, _) = points[*j as usize] ;

            //multiply by (x-xj)
            let factor = vec![mod_sub(0,xj),1];
            basis = poly_mul(&basis , &factor) ;

            //multiply denom by xi - cj 
            denom = mod_mul(denom , mod_sub(xi, xj));
        
        
            }
        let inv_denom =  mod_inv(denom);
        
        let scale = mod_mul(yi , inv_denom);
         

        let term = poly_scalar_mul( &basis , scale);
        

        result = poly_add(&result , &term);
        

    }
    result
}

//polynomial evaluation 
pub fn evaluate_poly (poly: &Vec<i64>, val: i64) -> i64 {
    let mut res: i64 = 0 ;

    for x in 0..poly.len(){
             res = mod_add(res ,  mod_mul( poly[x] , mod_pow(val , x as i64) ))
                 
    }
    let resf: i64 = if res < 0 { MOD + res } else {res };
    resf }

    
    

// prity print 
fn print_poly(poly: &Vec<i64>) {
    for (i , coeff) in poly.iter().enumerate(){
        if *coeff != 0 {
            println!("{}x^{}", coeff,i)
        }}}

