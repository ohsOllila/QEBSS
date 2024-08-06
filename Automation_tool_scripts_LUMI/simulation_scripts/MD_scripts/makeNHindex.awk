#!/bin/gawk -f

##########################################
# Set initial things
##########################################
BEGIN {

      #-- Variables
      ind = 0      
      frame = 0
      #-- For reading .gro-format
      FIELDWIDTHS = "5 5 5 5 8 8 8 8 8 8"
}

##########################################
# Read .gro file into memory line by line
##########################################
{
      #-- First line: comment
    if($1~"frame"){
	frame++;
	ind=0;
	next
    }
    if(NR == 1) {
	COMMENT=$0
	next
      }
 
      #-- Second line: number of atoms
      if(NR == 2) {
         NATOMS=$0
         next
      }

      #-- Last line: box dimensions 
      if(NR == NATOMS+3) {
         boxx = substr($0,  1, 10) + 0.0
         boxy = substr($0, 11, 20) + 0.0
         boxz = substr($0, 21, 30) + 0.0
         next
      }

      #-- Other lines: atom coordinates
      ind++
      MOLNUM[ind] = $1 
      MOLNAME[ind] = $2 
      ATOMNAME[ind] = $3 
      ATOMNUM[ind] = $4 
      X[ind] = $5 
      Y[ind] = $6 
      Z[ind] = $7 
     
}

##########################################
# Write .gro file
##########################################
END {

  

      #-- Print header
#      print COMMENT    
#      print NATOMS 
      #-- Print atoms 
      run_num = 0
      for(ind = 1; ind <= NATOMS; ind++) {
	  if(ATOMNAME[ind]=="    N" && ATOMNAME[ind+1]=="   HN"){
	      print "["MOLNUM[ind]" ]"
	      print ATOMNUM[ind] ATOMNUM[ind+1]
	  }
#            run_num++
#            printf("%5s", ATOMNAME[ind])
#            printf("%8.3f%8.3f%8.3f\n", X[ind], Y[ind], Z[ind])
      }

      #-- Print box dimensions
#      printf("%10.5f%10.5f%10.5f\n", boxx, boxy, boxz)      

}
