use rs_merkle::{ MerkleTree , MerkleProof , Hasher ,algorithms::Sha256 };

pub fn root_prod (data: &Vec<i64>) -> [u8 ; 32 ]{
    let leaves: Vec<[u8; 32]> = data.iter().map(|x| Sha256::hash(&x.to_le_bytes())).collect();

    let merkle_tree =  MerkleTree::<Sha256>::from_leaves(&leaves);
    let root_hash = merkle_tree.root().expect("Tree should not be empty");
    println!("Merkle Root:{:?}" , root_hash);
    root_hash 
}
