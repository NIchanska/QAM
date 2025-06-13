function findLinear(lin_funcs)
    linearClasses := [];
    FieldF := {x : x in FF | x ne 0};
    while #FieldF ne 0 do;
	el := Random(FieldF);
	Sel := {Evaluate(l, el) : l in lin_funcs};
	Append(~linearClasses,Sel);
	FieldF := FieldF diff Sel;
    end while;
    return linearClasses;
end function;


function linearB(A, linfuncs)
    return {l : l in linfuncs | Evaluate(l,A) eq A};
end function;

function equivSet(set)
    out := [];
    for i in set do
        out:=Append(out, FF!Min(i));
    end for;
    out2 := [];
    for el in out do;
	mina := FF!Min(out);
	Append(~out2,mina);
	Exclude(~out,mina);
    end for;
    return out2;
end function;

function classORB(A, linfuncs)
    perms := linearB(A, linfuncs);
    //"b perms:", #perms;
    class := findLinear(perms);
    return equivSet(class),perms;
end function;

function estimateTreeSize(LP,level,number_of_variables,max_classes)
	/* Check if we are in a lef */
	if level ge number_of_variables then
		return 1;
	end if;

	total_size := 0;

	/* Partition according to LPs */
	classes := findLinear(LP);
	reps := equivSet(classes);

	/* Check termination condition */
	if #classes gt max_classes then
		/* n here is a global variable denoting the dimension of the field */
		return (2^n)^(number_of_variables-level-1) * #classes;
	
	end if;

	/* If the number of classes is small enough, handle this sub-tree recursively */
	for VAL in reps do
		LPprime := linearB(VAL,LP);
		total_size +:= $$(LPprime,level+1,number_of_variables,max_classes);
	end for;

	return total_size;
end function;
