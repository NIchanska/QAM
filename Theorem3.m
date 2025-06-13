Iset := Subsets({1..n});
Iset0 := Iset;

IndSet := [];

while #Iset ne 0 do;
  el := Random(Iset);
  Append(~IndSet,el);
  Exclude(~Iset,el);
  for i in {1..n-1 by m} do;
    el1 := {(j+i) mod n : j in el};
    if 0 in el1 then;
      Exclude(~el1,0);
      Include(~el1,n);
    end if;
    Exclude(~Iset,el1);
  end for;
end while;

