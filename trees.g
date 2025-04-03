#
# Building the enumeration trees
#

# S_n

Make_Tree:=function(n)
  local tree,newlevel;
  tree:=[];

  # function inserting the descendants 

  newlevel:=function(sg)
    local k,s,c,i,m;
    k:=Length(sg)+1;
    s:=Length(tree);
    c:=0;
    for i in [1..n] do
      if not (i in sg) then 
        m:=Length(tree)+1;
        Add(tree,rec(k:=k,m:=i,bro:=0,son:=0));
        if c<>0 then
          tree[c].bro:=m;
        fi;
        if k<n then 
          tree[m].son:=m+1;
          newlevel(Concatenation(sg,[i]));
        fi;
        c:=m;
      fi;
    od;
  end;

  # creating the tree

  newlevel([]);
  return tree;
end;

# saving the trees in a file

PrintTo("data/trees/small4.g","tree:=",Make_Tree(4),";\n");
PrintTo("data/trees/small6.g","tree:=",Make_Tree(6),";\n");


# Big trees

Make_Big_Tree:=function(none,right,left,both)
  local n,L,R,tree,biglevel;
  n:=none+right+left+both;
  L:=left+both;
  R:=right+both;
  tree:=[];

  # biglevel

  biglevel:=function(k,sg,tau)
    local s,c,i,j,m;
    s:=Length(tree);
    c:=0;
    if k<=none then
      Add(tree,rec(k:=k,l:=0,r:=0,bro:=0,son:=s+2));
      biglevel(k+1,sg,tau);
    elif k<=none+right then
      for i in [1..R] do
        if not (i in tau) then 
          m:=Length(tree)+1;
          Add(tree,rec(k:=k,l:=0,r:=i,bro:=0,son:=0));
          if c<>0 then
            tree[c].bro:=m;
          fi;
          if k<n then 
            tree[m].son:=m+1;
            biglevel(k+1,sg,Concatenation(tau,[i]));
          fi;
          c:=m;
        fi;
      od;
    elif k<=none+right+left then
      for j in [1..L] do
        if not (j in sg) then 
          m:=Length(tree)+1;
          Add(tree,rec(k:=k,l:=j,r:=0,bro:=0,son:=0));
          if c<>0 then
            tree[c].bro:=m;
          fi;
          if k<n then 
            tree[m].son:=m+1;
            biglevel(k+1,Concatenation(sg,[j]),tau);
          fi;
          c:=m;
        fi;
      od;
    else
      for i in [1..L] do
        if not (i in sg) then 
          for j in [1..R] do
            if not (j in tau) then
              m:=Length(tree)+1;
              Add(tree,rec(k:=k,l:=i,r:=j,bro:=0,son:=0));
              if c<>0 then
                tree[c].bro:=m;
              fi;
              if k<n then 
                tree[m].son:=m+1;
                biglevel(k+1,Concatenation(sg,[i]),
                             Concatenation(tau,[j]));
              fi;
              c:=m;
            fi;
          od;
        fi;
      od;
    fi;
  end;

  # creating the tree

  biglevel(1,[],[]);
  return tree;
end;

# making the list of shapes

shapes:=Cartesian([0..8],[0..8],[0..8],[0..8]);
shapes:=Filtered(shapes,s->Sum(s)=8);
shapes:=Filtered(shapes,s->s[2]<=s[3]);
shapes:=Filtered(shapes,s->s[2]+s[4] in [4,6]);
shapes:=Filtered(shapes,s->s[3]+s[4] in [4,6]);

# saving the big trees in a file

for s in shapes do
  none:=s[1];
  right:=s[2];
  left:=s[3];
  both:=s[4];
  file_name:=Concatenation("data/trees/big",String(none),String(right),String(left),String(both),".g");
  PrintTo(file_name,"big_tree:=",Make_Big_Tree(none,right,left,both),";");
od;
