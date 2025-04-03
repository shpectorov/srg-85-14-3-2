# srg-85-14-3-2
Enumeration of SRG(85,14,3,2)

The code in this directory was used for the enumeration of strongly regular graphs with parameters (85,14,3,2).

- segments.g creating the lists of segments, segment pairs, and for each segment pair the list of compatible 
  third segments
- type.g an auxiliary function forming the type of segment
- trees.g generating small and big enumeration trees
- addone.g iterative LDLT algorithm
- enumeration13.g main code realising the four-step enumeration algorithm

One file is currently missing: filtering out segment triples without a favourite segment (Section 9 in the paper.)

The code was run each time on a certain interval of the entire range of cases [1..86333]. This interval can be set 
in line 930 of enumeration13.g.

If you want to run this code, it might be a good idea to consult us, as the code assumes a certain directory 
structure, which can, in principle, be gleaned from the code itself, but it might be easier to talk to us. 
