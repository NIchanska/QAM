function SubsCleaner(Ss);
        Subs := {};
        for Elem in Ss do
                Include(~Subs,[{ar : ar in Elem[1]},{rr : rr in Elem[2]}]);
        end for;
        for el1 in Subs do;
                Eqset := {[{(ee + k) mod n: ee in el1[1]},{(ee + k) mod n: ee in el1[2]}]: k in [1..n-1]};
                for el2 in Subs diff {el1} do
                        if el2 in Eqset then;
                                Exclude(~Subs,el2);
                        end if;
                end for;
                Include(~Subs,el1);
        end for;
        Subs := SetToSequence(Subs);
        S := [];
        for i in [1..#Subs] do;
                Append(~S,[SetToSequence(Subs[i][1]), SetToSequence(Subs[i][2])]);
        end for;
        return S;
end function;
