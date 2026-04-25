use rs_merkle::{ MerkleTree , MerkleProof , Hasher ,algorithms::Sha256 };

pub fn root_prod (data: &Vec<i64>) -> [u8 ; 32 ]{
    let leaves: Vec<[u8; 32]> = data.iter().map(|x| Sha256::hash(&x.to_le_bytes())).collect();

    let merkle_tree =  MerkleTree::<Sha256>::from_leaves(&leaves);
    let root_hash = merkle_tree.root().expect("Tree should not be empty");
    println!("Merkle Root:{:?}" , root_hash);
    root_hash 
 }

pub fn proof_prod (data: &Vec<i64> , index: usize) -> MerkleProof<Sha256>{
    let leaves: Vec<[u8; 32]> = data.iter().map(|x| Sha256::hash(&x.to_le_bytes())).collect();

    let merkle_tree =  MerkleTree::<Sha256>::from_leaves(&leaves);
    let indeces_to_prove = vec![index] ;
    let merkle_proof = merkle_tree.proof(&indeces_to_prove) ;
    let proof_hashes = merkle_proof.to_bytes() ;
    println!("Proof type:{:?}" , std::any::type_name_of_val(&proof_hashes));
    merkle_proof 
}

pub fn verify_proof (proof: MerkleProof<Sha256> , root: [u8 ;32] , val: i64 , index: usize , eval_order: usize)->bool {
    let indeces = vec![index] ;
    let hashed = vec![Sha256::hash(&val.to_le_bytes())] ; 
    proof.verify( root ,&indeces ,&hashed , eval_order)
}
