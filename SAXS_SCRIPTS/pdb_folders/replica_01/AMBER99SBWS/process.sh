for pdbfile in *.pdb
do
  crysol $pdbfile  --explicit-hydrogens --alternative-names
done

