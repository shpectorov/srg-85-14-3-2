#
# extending the LDLT decomposition by one dimension
# if semi positive definite
#
# based on Madeleine Whybrow's code
#

G:=List([1..30],i->List([1..30],j->0));
L:=List([1..40],i->List([1..40],j->0));
D:=[];

R:=List([1..30],i->List([1..30],j->0));
Id:=IdentityMat(30);

AddOne:=function(r)
  local n,i,sum,j;
  n:=Length(r);

  for i in [1..n] do
    sum:=0;
    for j in [1..i-1] do
      sum:=sum+L[n][j]*L[i][j]*D[j];
    od;
    if i<n then
      if D[i]=0 then
        if r[i]=sum then
          L[n][i]:=0;
        else
          return false;
        fi;
      else
        L[n][i]:=(r[i]-sum)/D[i];
      fi;
    else
      L[n][n]:=1;
      D[n]:=r[n]-sum;
      if D[n]<0 then
        return false;
      fi;
    fi;
  od;
  if n<=30 then
    for i in [1..n] do
      G[n][i]:=r[i];
      G[i][n]:=r[i];
    od;
    R[n]:=Id[n];
    for j in [1..n-1] do
      R[n]:=R[n]-L[n][j]*R[j];
    od;
  fi;

  return true;
end;
