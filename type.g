#
# Type functions
#

Left:=function(type)
  return type[3]+type[4];
end;

Right:=function(type)
  return type[2]+type[4];
end;

Type:=function(type)
  return type[2]+type[3]+2*type[4];
end;

