// calculate , with a trace , th fibbonacci up to the 100th 
pub fn nthfib (n: usize ) -> Vec<i64> {
   let  mut  coll: Vec<i64> = vec![ 1 , 1 ] ;

    for x in 2..(n + 0) {
        coll.push(coll[x - 1] + coll[x - 2])
    }

    coll 
}

