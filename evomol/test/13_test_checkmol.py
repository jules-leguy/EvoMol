from evomol.evaluation_entropy import extract_shingles

# print(extract_checkmol(MolGraph(MolFromSmiles("CC(=O)CC(=O)C"))))

print(extract_shingles("c1ccc1", 1))

print(extract_shingles("c1ccc1", 2))
