#
# The main enumeration part
#

# randomize GAP
Reset(GlobalMersenneTwister,CurrentDateTimeString());

#
# Global variables
#

# exact: where we found exactly 38 vertex sets
exact:=[];

# begt, endt: time monitoring variables
begt:=0;
endt:=0;

# print_string: array of what we print for each harder case
# PrintCase: print function
print_string:=[];
PrintCase:=function()
  Print("\r                                                                            \c");
  Print("\r",print_string[1],print_string[2],print_string[3],"\c");
end;

#
# Case in hand
#

# p: segment pair number
p:=0;

# n: first segment number
# k: second segment number
# first: historic first segment
n:=0;
k:=0;
first:=0;

# im: how the pair is glued together ([1,2] or [2,1])
im:=[];

# s: first segment
# t: second segment
s:=rec();
t:=rec();

# match: size of matching between s and t
match:=0;

# tree: tree of matchings between s and t
# u: node in tree
tree:=[];
u:=0;

# segment_triples: list of triples
# z: pair within the triple
# w: third segment
segment_triples:=[];
z:=0;
w:=rec();

# big_tree: tree of matchings between w and s plus t
# v: node in big_tree
big_tree:=[];
v:=0;

# Gram matrix and its decomposition
# G,L,D (such that G=LDL^t) in the list form, R=L^-1
#
# P the projection matrix
G:=List([1..30],i->List([1..30],j->0));
L:=List([1..30],i->List([1..30],j->0));
D:=[];
R:=List([1..30],i->List([1..30],j->0));

P:=List([1..30],i->List([1..30],j->0));

# lists of vertices of the three components
seg1:=[1..12];
seg2:=[1,2,13,14,15,16,17,18,19,20,21,22];
seg3:=[3,4,13,14,23,24,25,26,27,28,29,30];

# downs: precomputed sets of possible neighbours for further 
#        vertices adjacent to 13
downs:=[];
for p1 in Combinations(seg1,2) do
  for p2 in Combinations(seg2,2) do
    if 13 in p2 then
      if Intersection(p1,[1,2])=Intersection(p2,[1,2]) then
        for p3 in Combinations(seg3,2) do
          if 13 in p3 then
            if Intersection(p1,[3,4])=Intersection(p3,[3,4]) and
               Intersection(p2,[14])=Intersection(p3,[14]) then
              Add(downs,rec(neibs:=Union(p1,p2,p3)));
            fi;
          fi;
        od;
      fi;
    fi;
  od;
od;

# further: list of further vertices adjacent to 13
further:=[];

# loading segment data

Print("Loading segment data...\n");

if not IsBound(segment_data) then 
  segments:=[];
  Read("data/segment_data.g");
fi;

Print("...done\n");

# loading doubles

Print("Loading segment pairs...\n");

if not IsBound(segment_pairs) then 
  segment_pairs:=[];
  Read("data/segment_pairs.g");
fi;

Print("...done\n");

# loading the LDLT function AddOne

Print("Loading the AddOne function...\n");

AddOne:=function(r)
  return 1;
end;

Read("addone.g");

Print("...done\n");

#
# Level 3+4 functions
#

# Vertex projection

ComputeProjMat:=function()
  local i;
  P:=0*IdentityMat(30);
  for i in [1..30] do
    if D[i]<>0 then
      P:=P+TransposedMat([R[i]])*[R[i]]/D[i];
    fi;
  od;
end;

# checking whether an element of downs can be a further vertex

IsGoodVert:=function(d)
  local i,proj;
  for i in [1..30] do
    if D[i]=0 then
      if d.r*R[i]<>0 then
        return false;
      fi;
    fi;
  od;
  proj:=d.r*P;
  if d.r*proj>1 then
    return false;
  fi;
  d.proj:=proj;
  return true;
end;

# additional check for a further vertex
#   - no triple removed vertices
#   - computing the halo

NoTripleRemoved:=function(d)
  local h,i,j,c;
  h:=[];
  for i in [1..30] do
    if i in d.neibs then
      continue;
    fi;
    c:=0;
    for j in d.neibs do
      if G[i][j]=2/7 then
        c:=c+1;
      fi;
    od;
    if c>=3 then
      return false;
    elif c=2 then
      Add(h,i);
    fi;
  od;
  d.h:=h;
  return true;
end;

# check further vertices for compatibility

AreGoodNeighbours:=function(d,e)
  if Length(Intersection(d.neibs,e.neibs))>3 then 
    return false;
  fi;
  if Length(Intersection(d.neibs,e.h))<>0 or
     Length(Intersection(d.h,e.neibs))<>0 then
    return false;
  fi;
  if (2/7-d.r*e.proj)^2>(1-d.r*d.proj)*(1-e.r*e.proj) then
    return false;
  fi;
  return true;
end;

AreGoodNonNeighbours:=function(d,e)
  if Length(Intersection(d.neibs,e.neibs))>2 then 
    return false;
  fi; 
  if (-1/14-d.r*e.proj)^2>(1-d.r*d.proj)*(1-e.r*e.proj) then
    return false;
  fi;
  return true;
end;

AreCompatibleVerts:=function(d,e)
  return AreGoodNeighbours(d,e) or AreGoodNonNeighbours(d,e);
end;

#
# Step 4
#
# checking the eight further vertices for compatibility
# 

# number of all neighbours of 13

tot:=14;

# edge options
yes:=1;
maybe:=0;
no:=-1;

# the iterator going through options

Step4Iterator:=function(dem,sup,edge)
  local ndem,nsup,nedge,change,i,j,e,k,M,r;

  ndem:=StructuralCopy(dem);
  nsup:=StructuralCopy(sup);
  nedge:=StructuralCopy(edge);

  # we iterate the checks
  change:=true;
  while change do
    change:=false;

    # vertex has all edges already
    for i in [1..tot] do
      if ndem[i]=0 and nsup[i]>0 then
        for j in [1..tot] do
          if nedge[i][j]=maybe then
            nedge[i][j]:=no;
            nedge[j][i]:=no;
            nsup[i]:=nsup[i]-1;
            nsup[j]:=nsup[j]-1;
            if nsup[j]<ndem[j] then
              return;
            fi;
            change:=true;
          fi;
        od;
      fi;
    od;

    # vertex needs all maybes
    for i in [1..tot] do
      if nsup[i]=ndem[i] and nsup[i]>0 then
        for j in [1..tot] do
          if nedge[i][j]=maybe then
            nedge[i][j]:=yes;
            nedge[j][i]:=yes;
            nsup[i]:=nsup[i]-1;
            nsup[j]:=nsup[j]-1;
            ndem[i]:=ndem[i]-1;
            ndem[j]:=ndem[j]-1;
            if ndem[j]<0 then
              return;
            fi;
            change:=true;
          fi;
        od;
      fi;
    od;

    # number of common neighbours
    for i in [1..tot-1] do
      for j in [i+1..tot] do
        e:=0;
        for k in [1..tot] do
          if not (k in [i,j]) then
            if nedge[i][k]=yes and nedge[j][k]=yes then
              e:=e+1;
            fi;
          fi; 
        od;
        # too many
        if e>2 then
          return;
        fi;
        # only a 4-clique  
        if e=2 then
          if nedge[i][j]=no then
            return;
          fi;
          if nedge[i][j]=maybe then
            change:=true;
            nedge[i][j]:=yes;
            nedge[j][i]:=yes;
            nsup[i]:=nsup[i]-1;
            nsup[j]:=nsup[j]-1;
            ndem[i]:=ndem[i]-1;
            ndem[j]:=ndem[j]-1;
            if ndem[i]<0 or ndem[j]<0 then
              return;
            fi;
          fi;
          # and no more common neibs
          if nedge[i][j]=yes then
            for k in [1..tot] do
              if nedge[i][k]=yes and nedge[j][k]=maybe then
                change:=true;
                nedge[j][k]:=no;
                nedge[k][j]:=no;
                nsup[j]:=nsup[j]-1;
                nsup[k]:=nsup[k]-1;
                if nsup[j]<ndem[j] or nsup[k]<ndem[k] then
                  return;
                fi;
              fi;         
              if nedge[i][k]=maybe and nedge[j][k]=yes then
                change:=true;
                nedge[i][k]:=no;
                nedge[k][i]:=no;
                nsup[i]:=nsup[i]-1;
                nsup[k]:=nsup[k]-1;
                if nsup[i]<ndem[i] or nsup[k]<ndem[k] then
                  return;
                fi;
              fi;         
            od;
          fi;
        fi;
        if e=1 then
          if nedge[i][j]=no then
            for k in [1..tot] do
              if nedge[i][k]=yes and nedge[j][k]=maybe then
                change:=true;
                nedge[j][k]:=no;
                nedge[k][j]:=no;
                nsup[j]:=nsup[j]-1;
                nsup[k]:=nsup[k]-1;
                if nsup[j]<ndem[j] or nsup[k]<ndem[k] then
                  return;
                fi;
              fi;         
              if nedge[i][k]=maybe and nedge[j][k]=yes then
                change:=true;
                nedge[i][k]:=no;
                nedge[k][i]:=no;
                nsup[i]:=nsup[i]-1;
                nsup[k]:=nsup[k]-1;
                if nsup[i]<ndem[i] or nsup[k]<ndem[k] then
                  return;
                fi;
              fi;         
            od; 
          fi;
        fi;        
      od;
    od;
  od;

  # if we find maybe then iterate
  for i in [1..tot-1] do
    for j in [i+1..tot] do
      if nedge[i][j]=maybe then
        nedge[i][j]:=no;
        nedge[j][i]:=no;
        nsup[i]:=nsup[i]-1;
        nsup[j]:=nsup[j]-1;
        Step4Iterator(ndem,nsup,nedge);
        nedge[i][j]:=yes;
        nedge[j][i]:=yes;
        ndem[i]:=ndem[i]-1;
        ndem[j]:=ndem[j]-1;
        Step4Iterator(ndem,nsup,nedge);
        return;        
      fi;
    od;
  od;

  # check positive definiteness and rank
  # compute Gram matrix in the perp of the triple
  M:=List([1..8],i->[]);
  for i in [1..8] do
    for j in [1..8] do
      if i=j then
        r:=1;
      elif nedge[6+i][6+j]=no then
        r:=-1/14;
      else
        r:=2/7;
      fi;
      M[i][j]:=r-further[i].proj*further[j].r;
    od;
  od;
  # an quick approximation of positive definiteness
  for i in [1..8] do
    if Determinant(M{[1..i]}{[1..i]})<0 then
      return;
    fi;
  od;
  # rank
  if Length(Filtered(D,d->d<>0))+Rank(M)>34 then
    return;
  fi;

  Print("\n==Passed Level 4 checks\n");
  Add(exact,StructuralCopy([G,further,edge]));
  return;

end;

# function preparing

Step4:=function()
  local b,edge,i,j,dem,sup;

  # neibours of 13 in the triple
  b:=Concatenation(Filtered(seg2,i->G[i][13]=2/7),Filtered(seg3,i->G[i][13]=2/7));
  
  # edges
  edge:=List([1..tot],i->List([1..tot],j->maybe));
  for i in [1..tot] do
    for j in [i..tot] do
      if i=j then
        edge[i][j]:=no;
      else
        if i=1 then
          if j in [2,3,4] then
            edge[i][j]:=yes;
            edge[j][i]:=yes;
          else
            edge[i][j]:=no;
            edge[j][i]:=no;
          fi;
        elif i=2 then
          if j in [5,6] then
            edge[i][j]:=yes;
            edge[j][i]:=yes;
          else
            edge[i][j]:=no;
            edge[j][i]:=no;
          fi;
        fi;
        if IsSubset([3..6],[i,j]) then
          if G[b[i-2]][b[j-2]]=2/7 then
            edge[i][j]:=yes;
            edge[j][i]:=yes;
          else
            edge[i][j]:=no;
            edge[j][i]:=no;
          fi;
        fi;
        if i in [3..6] and j>6 then
          if b[i-2] in further[j-6].neibs then
            edge[i][j]:=yes;
            edge[j][i]:=yes;
          else
            edge[i][j]:=no;
            edge[j][i]:=no;
          fi;
        fi;
        if i>6 then
          if not AreGoodNeighbours(further[i-6],further[j-6]) then
            edge[i][j]:=no;
            edge[j][i]:=no;
          fi;
          if not AreGoodNonNeighbours(further[i-6],further[j-6]) then
            edge[i][j]:=yes;
            edge[j][i]:=yes;
          fi;
        fi;
      fi;
    od;
  od;

  # demand in valency
  dem:=List([1..tot],i->3);
  for i in [1..tot] do
    for j in [1..tot] do
      if edge[i][j]=yes then
        dem[i]:=dem[i]-1;
        if dem[i]<0 then
          return;
        fi;
      fi;
    od;
  od;
  # supply of possible valencies
  sup:=List([1..tot],i->0);
  for i in [1..tot] do
    for j in [1..tot] do
      if edge[i][j]=maybe then
        sup[i]:=sup[i]+1;
      fi; 
    od;
    if sup[i]<dem[i] then
      return;
    fi;
  od;

  Step4Iterator(dem,sup,edge);
  return;
end;

#
# The iterative part of Step 3
#

# vertices that are not 13
others:=Difference([1..30],[13]);

Step3Iterator:=function(verts,all)
  local i,j,add,off,nall,nverts,d,e,done,p,q,good;

  # if insufficient further vertices
  if Length(further)+Length(verts)<tot then
    return;
  fi;  

  # can we find forced vertices to add?
  add:=[];
  if Length(further)+Length(verts)=tot then
    add:=StructuralCopy(verts);
  else
    # compute the offer (only involving 13)
    off:=[];
    for i in others do
      off[i]:=0;
    od;
    for d in verts do
      for i in d.neibs do
        if i<>13 then
          off[i]:=off[i]+1;
        fi;
      od;
    od;

    # compare offer and demand
    for i in others do
      p:=Set([i,13]);
      if all[p[1]][p[2]]<>0 then
        if off[i]<all[p[1]][p[2]] then
          return;
        elif off[i]=all[p[1]][p[2]] then
          add:=Union(add,
                 Filtered(verts,e->i in e.neibs));
        fi;
      fi;
    od;
  fi;

  if add<>[] then
    if Length(further)+Length(add)>tot then
      return;
    fi;

    # are additional further vertices compartible?
    for i in [1..Length(add)-1] do
      for j in [i+1..Length(add)] do
        if not AreCompatibleVerts(add[i],add[j]) then
          return;
        fi;
      od;
    od;

    # compute new all
    nall:=StructuralCopy(all);
    for d in verts do
      for p in Combinations(d.neibs,2) do
        i:=p[1];
        j:=p[2];
        nall[i][j]:=nall[i][j]-1;
        if nall[i][j]<0 then
          return;
        fi;
      od;
    od;

    # if we are in a potential exact situation
    if Length(further)+Length(add)=tot then
      further:=Concatenation(further,add);
      Step4();
      further:=further{[1..Length(further)-Length(add)]};
      return;
    fi;

    # compute new verts
    nverts:=[];
    for d in verts do
      if d in add then 
        continue; 
      fi;
      good:=true;
      for e in add do 
        if not AreCompatibleVerts(d,e) then 
          good:=false;
          break;
        fi;
      od;
      if not good then
        continue;
      fi;
      for p in Combinations(d.neibs,2) do
        if nall[p[1]][p[2]]=0 then 
          good:=false;
          continue; 
        fi;
      od;
      if good then
        Add(nverts,d);
      fi;
    od;
    
    # iterate
    further:=Concatenation(further,add);
    Step3Iterator(nverts,nall);
    further:=further{[1..Length(further)-Length(all)]};
    return;
  fi;

  # now the case where add=[]

  # select a good vertex
  j:=0;
  for i in others do
    if (i<13 and all[i][13]=0) or (i>13 and all[13][i]=0) then
      continue;
    fi;
    if j=0 then
      j:=i;
    else
      p:=Set([i,13]);
      q:=Set([j,13]);
      if all[p[1]][p[2]]<all[q[1]][q[2]] then
        j:=i;
      fi;
      if all[p[1]][p[2]]=1 then
        break;
      fi;
    fi;
  od;
  if j=0 then
    return;
  fi;
    
  # now we have a good vertex j
  done:=[];
  for e in verts do
    if not (j in e.neibs) then
      continue;
    fi;

    # force a vertex

    # compute new all
    nall:=StructuralCopy(all);
    for p in Combinations(e.neibs,2) do
      all[p[1]][p[2]]:=all[p[1]][p[2]]-1;
    od;
    
    # if we are in a potential exact situation
    if Length(further)+1=tot then
      Add(further,e);
      Step4();
      Add(done,e);
      Unbind(further[Length(further)]);
      continue;
    fi;

    # compute new verts
    nverts:=[];
    for d in verts do
      if d=e or d in done then
        continue;
      fi;
      if not AreCompatibleVerts(d,e) then
        continue;
      fi;
      good:=true;
      for p in Combinations(d.neibs,2) do
        if nall[p[1]][p[2]]=0 then
          good:=false;
          break;
        fi;
      od;
      if good then
        Add(nverts,d);
      fi;
    od;

    # update further and iterate
    Add(further,e);
    Step3Iterator(nverts,nall);
    Add(done,e);
    Unbind(further[Length(further)]);
  od;

  return;    
end;

#
# Step 3 preparation
#

Step3:=function()
  local all,i,j,ij,verts,d,e,good,p,r;

  further:=[];

  ComputeProjMat();

  # compute demand for all pairs
  all:=List([1..29],i->[]);
  for i in [1..29] do
    for j in [i+1..30] do
      if G[i][j]=2/7 then
        all[i][j]:=3;
      else
        all[i][j]:=2;
      fi;
      for ij in [1..30] do
        if G[i][ij]=2/7 and G[j][ij]=2/7 then
          all[i][j]:=all[i][j]-1;
        fi;
      od;
      if IsSubset(seg1,[i,j]) then
        all[i][j]:=all[i][j]-1;
      fi;
      if IsSubset(seg2,[i,j]) then
        all[i][j]:=all[i][j]-1;
      fi;
      if IsSubset(seg3,[i,j]) then
        all[i][j]:=all[i][j]-1;
      fi;
      if all[i][j]<0 then
        return;
      fi;
    od;
  od;

  # set up possible further vertices
  verts:=[];

  for d in downs do
    e:=StructuralCopy(d);

    # check pairs
    good:=true;
    for p in Combinations(e.neibs,2) do
      if all[p[1]][p[2]]=0 then
        good:=false;
        break;
      fi;
    od;
    if not good then
      continue;
    fi;

    # check the no-triple-removed condition
    if not NoTripleRemoved(e) then
      continue;
    fi;

    # set up right side
    e.r:=List([1..30],i->-1/14);
    for i in e.neibs do
      e.r[i]:=2/7;
    od;

    # check the projection
    if IsGoodVert(e) then
      Add(verts,e);
    fi;
  od;
  verts:=Set(verts);

# print if large case
  if Length(verts)>=150 then
    print_string[3]:=Concatenation(" (",String(Length(verts))," verts)");
    PrintCase();
  fi;
 
  Step3Iterator(verts,all);
    
  return;
end;

#
# Matchings enumeration to be used at Step 2
#

gram2:=[[]];
im1:=[];
im2:=[];
lset:=[];
rset:=[];

Step2_node:=function(v)
  local k,r,ll,rr,good,son;

  k:=4+big_tree[v].k;

  r:=List([1..22],j->-1/14);

  r[2+1]:=gram2[k][im1[1]];
  r[2+2]:=gram2[k][im1[2]];
  r[12+1]:=gram2[k][im2[1]];
  r[12+2]:=gram2[k][im2[2]];

  ll:=big_tree[v].l;
  if ll<>0 then
    r[lset[ll]]:=2/7;
  fi;

  rr:=big_tree[v].r;
  if rr<>0 then
    r[rset[rr]]:=2/7;
  fi;

  Append(r,gram2[k]{[5..k]});
  good:=AddOne(r);
  if good then
    son:=big_tree[v].son;
    if son=0 then
      print_string[2]:=Concatenation(", big_tree ",String(v),"/",String(Length(big_tree)));
      print_string[3]:="";
      PrintCase();
      Step3();
    else
      while son<>0 do
        son:=Step2_node(son);
      od;
    fi;
  fi;
  return big_tree[v].bro;

end;

#
# Preparing Step 2
#

Step2:=function()
  local c,type,none,right,left,both;

  c:=segment_triples[z];

  w:=segments[c.m];
  gram2:=w.gram;
  type:=w.type;
  none:=type[1];
  right:=type[2];
  left:=type[3];
  both:=type[4];

  Read(Concatenation("data/trees/big",String(none),String(right),String(left),String(both),".g"));

  im1:=c.im1;
  im2:=c.im2;

  type:=s.type;
  none:=type[1];
  right:=type[2];
  left:=type[3];
  both:=type[4];
  lset:=4+Concatenation([none+1..none+right],[9-both..8]);

  type:=t.type;
  none:=type[1];
  right:=type[2];
  left:=type[3];
  both:=type[4];
  rset:=14+Concatenation([none+1..none+right],[9-both..8]);
    
  v:=1;
  while v<>0 do
    v:=Step2_node(v);
  od;
end;

#
# Matchings enumeration to be used at Step 1
#

Step1_node:=function(v)
  local k,r,son,good;

  k:=12-match+tree[v].k;
  r:=List([1..12],j->-1/14);
  r[1]:=t.gram[k][im[1]];
  r[2]:=t.gram[k][im[2]];
  r[12-match+tree[v].m]:=2/7;  
  Append(r,t.gram[k]{[3..k]});
  good:=AddOne(r);
  if good then
    son:=tree[v].son;
    if son=0 then
      print_string[1]:=Concatenation("matching ",String(v),"/",String(Length(tree)));
      Step2();
    else
      while son<>0 do
        son:=Step1_node(son);
      od;
    fi;
  fi;
  return tree[v].bro;
end; 


#
# main enumeration: Step 1
#

for p in [9901..10000] do

  workname:=Concatenation("data/working/",String(p),".g");
  outname:=Concatenation("data/record/out",String(p),".g");

  if IsExistingFile(workname) or IsExistingFile(outname) then
    continue;
  fi;

  PrintTo(workname,"Working...\n");

  Print("\n\n\n\n======= Starting segment pair ",p,"\n");

  allexact:=[];

# reading in the triples

  Read(Concatenation("data/abrev/",String(p),".g"));

# selecting a triple

  for z in [1..Length(segment_triples)] do

    exact:=[];

    filename:=Concatenation("data/record/",String(p),"_",String(z),".g");  
    if IsExistingFile(filename) then
      Read(filename);
      Append(allexact,exact);
      RemoveFile(filename);
      continue;
    fi; 

    Print("\n\n\nStarting segment pair ",p,", triple ",z,"/",Length(segment_triples),"\n\n");

    begt:=Runtime();

    # setting up the run

    n:=segment_pairs[p].n;
    k:=segment_pairs[p].k;
    im:=segment_pairs[p].im;

    # update the first segment if necessary

    if n<>first then
      first:=n;
      s:=segments[n];
      for i in [1..12] do
        r:=s.gram[i]{[1..i]};
        zero:=AddOne(r);
      od;
    fi;

    # add the unmatched part of second segment

    t:=segments[k];
    match:=Sum(t.type{[3,4]});

    for i in [3..12-match] do
      r:=[];
      r[1]:=t.gram[i][im[1]];
      r[2]:=t.gram[i][im[2]];
      Append(r,List([3..12],j->-1/14));
      Append(r,t.gram[i]{[3..i]});
      zero:=AddOne(r);
    od;

    # reading in the small tree

    Read(Concatenation("data/trees/small",String(match),".g"));

    # the matched part of second segment is added iteratively
    # in each case

    u:=1;
    while u<>0 do
      u:=Step1_node(u);
    od;

    endt:=Runtime();

    Print("\n");
    Print("Time used: ",StringTime(endt-begt),"\n");
    Print("Total cases still alive: ",Length(exact), "\n");

    Append(allexact,exact);
    Print("done\n");

  od;
  
  PrintTo(outname,"allexact=",allexact,";\n");
  RemoveFile(workname);
  Print("\n\nDone with segment pair ",p,"\n");
  
od;
