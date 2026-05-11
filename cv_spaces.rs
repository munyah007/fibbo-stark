use crate::complex::Comp ;
use crate::complex ;

pub fn add (vec1: &Vec<Comp> , vec2: &Vec<Comp>) -> Vec<Comp> {
    let mut cont: Vec<Comp> = Vec::new() ;
    for x in 0..vec1.len() {
        cont.push(complex::add(&vec1[x] ,&vec2[x])) ;
            }
            cont }



pub fn inv (vec1: &Vec<Comp> ) -> Vec<Comp> {
    let mut cont: Vec<Comp> = Vec::new() ;
    for x in 0..vec1.len() {
        cont.push(complex::inv(&vec1[x])) ;
            }
            cont }


pub fn scal_mul (scal: &Comp , vec: &Vec<Comp>) -> Vec<Comp> {
    let mut cont: Vec<Comp> = Vec::new() ;
    for x in 0..vec.len() {
        cont.push(complex::mul(&scal ,&vec[x])) ;
            }
            cont }

pub fn trans (vec1: &Vec<Vec<Comp>> ) -> Vec<Vec<Comp>> {
    let mut cont: Vec<Vec<Comp>> = Vec::new() ;
    for x in 0..vec1[0].len() {// outer is over columns 
                               // inner over rows
        let mut cont2: Vec<Comp> = Vec::new() ; 
        for y in 0..vec1.len() {
        cont2.push(vec1[y][x]) ;
            }
          cont.push(cont2);
    }
            cont }
pub fn conj_mat (vec1: &Vec<Vec<Comp>> ) -> Vec<Vec<Comp>> {
    let mut cont: Vec<Vec<Comp>> = Vec::new() ;
    for x in 0..vec1.len() {
        let mut cont2: Vec<Comp> = Vec::new() ; 
        for y in 0..vec1[0].len() {
        cont2.push(complex::conj(&vec1[x][y])) ;
            }
          cont.push(cont2);
    }
            cont }

// if the original n by n equals its dagger its called Hermitian
pub fn dagger (vec1: &Vec<Vec<Comp>> ) -> Vec<Vec<Comp>> {
    let  conjugate = conj_mat(vec1) ;
    
            trans(&conjugate) }

pub fn scal_mul_mat (scal: &Comp , vec: &Vec<Vec<Comp>>) -> Vec<Vec<Comp>> {
    let mut cont: Vec<Vec<Comp>> = Vec::new() ;
    for x in 0..vec.len() {
    let mut cont2: Vec<Comp> = Vec::new() ;
        for y in 0..vec[0].len(){
        cont2.push(complex::mul(&scal ,&vec[x][y])) ;
            }
            cont.push(cont2) ; }
    cont
}


pub fn inv_mat ( vec: &Vec<Vec<Comp>>) -> Vec<Vec<Comp>> {
    let mut cont: Vec<Vec<Comp>> = Vec::new() ;
    for x in 0..vec.len() {
    let mut cont2: Vec<Comp> = Vec::new() ;
        for y in 0..vec[0].len(){
        cont2.push(complex::inv(&vec[x][y])) ;
            }
            cont.push(cont2) ; }
    cont
}


pub fn add_mat ( vec: &Vec<Vec<Comp>> , vec2: &Vec<Vec<Comp>>) -> Vec<Vec<Comp>> {
    let mut cont: Vec<Vec<Comp>> = Vec::new() ;
    for x in 0..vec.len() {
    let mut cont2: Vec<Comp> = Vec::new() ;
        for y in 0..vec[0].len(){
        cont2.push(complex::add(&vec[x][y], &vec2[x][y])) ;
            }
            cont.push(cont2) ; }
    cont
}

pub fn mul_rows ( vec: &Vec<Comp> , vec2: &Vec<Comp>) -> Comp {
    let mut cont: Vec<Comp> = Vec::new() ;
    
    for x in 0..vec.len() {
    
        
        cont.push(complex::mul(&vec[x], &vec2[x])) ;
            
             }

    let mut add_cont: Comp = cont[0].clone() ;
    for x in 1..cont.len() {
        add_cont = complex::add(&add_cont ,&cont[x]) ;
    }

    add_cont 
    
}

pub fn mul_mat ( vec: &Vec<Vec<Comp>> , vec2: &Vec<Vec<Comp>>) -> Vec<Vec<Comp>> {
    let transposed = trans(vec2) ;
    let mut cont: Vec<Vec<Comp>> = Vec::new() ;
    for x in 0..vec.len() {
    let mut cont2: Vec<Comp> = Vec::new() ;
        for y in 0..transposed.len(){
        cont2.push(mul_rows(&vec[x], &transposed[y])) ;
            }
            cont.push(cont2) ; }
    cont
}
pub fn action ( vec: &Vec<Vec<Comp>> , vec2: &Vec<Comp>) -> Vec<Comp> {
    let matricised = vec![vec2.clone()] ;
    let transposed = trans(&matricised) ;
    let mut cont: Vec<Vec<Comp>> = Vec::new() ;
    for x in 0..vec.len() {
    let mut cont2: Vec<Comp> = Vec::new() ;
        for y in 0..transposed.len(){
        cont2.push(mul_rows(&vec[x], &transposed[y])) ;
            }
            cont.push(cont2) ; }
    cont[0].clone() 
}

pub fn c_vec_dot ( vec: &Vec<Comp> , vec2: &Vec<Comp>) -> Comp {
    let mut cont: Vec<Comp> = Vec::new() ;
    
    
    for x in 0..vec.len() {
       let conjg = complex::conj(&vec[x]) ;
        
        cont.push(complex::mul(&conjg, &vec2[x])) ;
            
             }

    let mut add_cont: Comp = cont[0].clone() ;
    for x in 1..cont.len() {
        add_cont = complex::add(&add_cont ,&cont[x]) ;
    }

    add_cont 
    
}
pub fn norm ( vec: &Vec<Comp> ) -> f64 {
    let normsqd = c_vec_dot(&vec , &vec).real ;

   let realp =  normsqd.sqrt();

   realp
    
}   

pub fn dist ( vec1: &Vec<Comp> , vec2: &Vec<Comp>) -> f64 {
    let inversed = inv(&vec2) ;
    let  diff = add(&vec1 , &inversed );
    let norm = norm(&diff) ;

norm
    
}




pub fn herm_test(mat: &Vec<Vec<Comp>>) -> bool {
    let dag = dagger(&mat) ;
    let test = &dag == mat ;

    test }

pub fn gen_id (n: usize) ->Vec<Vec<Comp>>{
    let mut cont: Vec<Vec<Comp>> = Vec::new() ;
    for x in 0..n {
        let mut cont2: Vec<Comp> = Vec::new() ;
        for y in 0..n{
            if y==x{
                cont2.push(Comp {real: 1f64 ,imag: 0f64}) ;}
                else {
                    cont2.push(Comp{real: 0f64 , imag: 0f64})}
              }
        cont.push(cont2)
    }

    cont }
pub fn normalize (x: f64) -> f64 {
    let eps = 1e-12 ;
    if (x - 1.0).abs() < eps {
        1.0}
    else if x.abs() < eps{
        0.0 }
    else {
        x}}
pub fn normalize_mat(mat: Vec<Vec<Comp>>) -> Vec<Vec<Comp>>{
 let mut cont: Vec<Vec<Comp>> = Vec::new() ;
    for x in 0..mat.len() {
        let mut cont2: Vec<Comp> = Vec::new() ;
        for y in 0..mat.len(){
        let newComp: Comp =complex::new(normalize( mat[x][y].real) , normalize(mat[x][y].imag)) ;  
        cont2.push(newComp) ;
        }
        cont.push(cont2)
    }

    cont }


pub fn unitarian_test(mat: &Vec<Vec<Comp>>) -> bool {
    let dag = dagger(&mat) ;
    let lhs = normalize_mat( mul_mat(&mat , &dag));
    let rhs = normalize_mat(mul_mat(&dag , &mat )) ;
    let size: usize  = mat.len() ;
    let id_mat = gen_id(size);
    let test = (&lhs == &id_mat) & (&rhs == &id_mat) ;

   println!("lhs {:?}" , lhs ) ;
     println!("rhs {:?}" , rhs ) ;

    test }

pub fn tens_mul_mat ( vec: &Vec<Vec<Comp>> , vec2: &Vec<Vec<Comp>>) ->Vec< Vec<Vec<Vec<Comp>>>> {
    
    let mut cont: Vec<Vec<Vec<Vec<Comp>>>> = Vec::new() ;
    for x in 0..vec.len() {
        let mut cont2: Vec<Vec<Vec<Comp>>>  = Vec::new() ;
        for y in 0..vec[x].len() {
            let scal_prod:  Vec<Vec<Comp>> = scal_mul_mat(&vec[x][y] ,&vec2);
            cont2.push(scal_prod) ; }
        cont.push(cont2) }
    cont
}
fn id (comp: &Comp) -> Comp{
    comp.clone()    }
// this assumes square matrices 
pub fn combine_row (coll: &Vec<Vec<Vec<Comp>>>) -> Vec<Vec<Comp>>{
    let rows = coll[0].len() ;
    let result: Vec<Vec<_>> = (0..rows).map(|r| {
                                     coll.iter().flat_map(|m| m[r].iter().copied()).collect()}).collect();

    result
}
    
pub fn tens_prod( vec: &Vec<Vec<Comp>> , vec2: &Vec<Vec<Comp>>)   ->Vec<Vec<Comp>>{
    let block: Vec<Vec<Vec<Vec<Comp>>>> = tens_mul_mat(vec , vec2) ;
    let mut cont1: Vec<Vec<Vec<Comp>>> = Vec::new() ;

    for x in 0..block.len() {
        cont1.push(combine_row(&block[x]));
    }

    let mut cont2: Vec<Vec<Comp>> = Vec::new() ;
    for x in 0..cont1.len() {
        for y in 0..cont1[0].len() {
            cont2.push(cont1[x][y].clone());
        }}

    cont2

}
