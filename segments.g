#
# Input type functions
#

Type:=function(type)
  return true;
end;

Read("type.g");

#
# Building the list of segments including a separate
# list of thirds
#

LoadPackage("grape");

goodgraphs:=[];
Read("data/goodgraphs.g");

# building the list of segments

Print("Building the list of segments...\n");

segments:=[];

for n in [1..Length(goodgraphs)] do
  G:=goodgraphs[n];
  A:=AutomorphismGroup(G);
  edges:=[];
  for a in [1..14] do
    for b in Adjacency(G,a) do
      if Intersection(Adjacency(G,a),Adjacency(G,b))=[] then 
        Add(edges,[a,b]);
      fi;
    od;
  od;
  orbs:=Orbits(A,edges,OnPairs);
  for o in orbs do
    e:=o[1];
    a:=e[1];
    b:=e[2];

    aa:=Set(Difference(Adjacency(G,a),[b]));
    bb:=Set(Difference(Adjacency(G,b),[a]));

    if aa[2] in Adjacency(G,aa[1]) then
      ta:=6;
    else
      ta:=4;
    fi;
    if bb[2] in Adjacency(G,bb[1]) then
      tb:=6;
    else
      tb:=4;
    fi;

    # we only use orbits of type other than [4,6]

    if ta=6 or tb=4 then

      mm:=Difference([1..14],Union(Set(e),aa,bb));
      ma:=Difference(mm,Union(Adjacency(G,aa[1]),
                              Adjacency(G,aa[2])));
      ua:=Difference(mm,ma);
      mb:=Difference(mm,Union(Adjacency(G,bb[1]),
                              Adjacency(G,bb[2])));
      ub:=Difference(mm,mb);

      none:=Intersection(ua,ub);
      right:=Intersection(ua,mb);
      left:=Intersection(ma,ub);
      both:=Intersection(ma,mb);

      type:=[Length(none),Length(right),
             Length(left),Length(both)];

      # we now reorder the segment

      vv:=Concatenation(aa,bb,none,right,left,both);
      H:=InducedSubgraph(G,vv);
      B:=AutomorphismGroup(H);
      C:=Stabilizer(B,[1,2],OnSets);
      if Length(Orbit(C,1))=2 then
        flip:=true;
      else
        flip:=false;
      fi;
      
      # create the Gram matrix for the segment
      
      gram:=List([1..12],i->List([1..12],j->-1/14));
      for i in [1..12] do
        gram[i][i]:=1;
        for j in Adjacency(H,i) do
          gram[i][j]:=2/7;
        od;
      od;

      # record the segment
                       
      Add(segments,rec(G:=n,
                       gram:=StructuralCopy(gram),
                       type:=ShallowCopy(type),flip:=flip));
    fi;
  od;
od;

StableSort(segments,function(x,y) if Type(x.type)>Type(y.type) 
                then return true; else return false; fi; end);

PrintTo("segment_data.g","segments:=",segments,";\n");

Print("...done\n");


#
# list of pairs
#

Print("Building the list of pairs...\n");

segment_pairs:=[]; 

# Each segment pair looks like [n,k,[c,d]], where n is 
# the first segment, k is the second segment and c and d 
# the images of the (first) pair of vertices under gluing 

for n in [1..Length(segments)] do
  s:=segments[n];
  for k in [n..Length(segments)] do
    t:=segments[k];
    if Type(s.type) in [10,12] and Type(t.type)=8 then
      continue;
    fi;
    Add(segment_pairs,rec(n:=n,k:=k,im:=[1,2]));
    if not s.flip and not t.flip then 
      Add(segment_pairs,rec(n:=n,k:=k,im:=[2,1]));
    fi;
  od;
od;

PrintTo("data/segment_pairs.g","segment_pairs:=",
                           segment_pairs,";\n");

Print("...done\n");

    
#
# List of triples
#

Print("Producing the list of triples files\n");

# For each pair p, segment_triples[p] looks like 
# a list of records [m,[c,d],[e,f]],
# where m is the third segment and 
# [c,d] and [e,f] are the images of pairs from 
# the segment pair p under gluing
#
# This list is saved in a separate file
# in the data/triples subdirectory

for p in [1..Length(segment_pairs)] do
  Print("                \rPair ",p,"\c");
  segment_triples:=[];
  n:=segment_pairs[p].n;
  k:=segment_pairs[p].k;
  tn:=Sum(segments[n].type{[2,4]});
  tk:=Sum(segments[k].type{[2,4]});
  for m in [k..Length(segments)] do
    if Sum(segments[m].type{[3,4]})<>tn or
       Sum(segments[m].type{[2,4]})<>tk then
      continue;
    fi;
    
    # aa=[1,2] and bb=[3,4] in every segment

    Add(segment_triples,rec(m:=m,im1:=[1,2],im2:=[3,4]));
    Add(segment_triples,rec(m:=m,im1:=[1,2],im2:=[4,3]));
    Add(segment_triples,rec(m:=m,im1:=[2,1],im2:=[3,4]));
    Add(segment_triples,rec(m:=m,im1:=[2,1],im2:=[4,3]));
  od;
  file_name:=Concatenation("data/triples/",String(p),".g");
  PrintTo(file_name,"segment_triples:=",segment_triples,";\n");
od;

Print("...done\n");

