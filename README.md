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
