/* classify if the system has a solution */
function Pows2l(el,n);
	return [el^(2^i): i in [0..n-1]];
end function;

function Pows2lm(el,n,m);
	return [el^(2^i): i in [0..n-1 by m]];
end function;

function Cat_Ind(n, m)
    F := GF(2^n);
    Vec,phi := VectorSpace(F,GF(2)); 
    Cat_m_I := {};

    for a in F do
        if a eq 0 then
            continue;
        end if;
        
        conjugates := Pows2lm(a,n,m);
        
        // Check if the conjugates are linearly independent
	conj_matr := Matrix(GF(2),#conjugates,n,[phi(conjugates[i]): i in [1..#conjugates]]);

        // If the conjugates are linearly independent, include a in Cat^m_I
        if Rank(conj_matr) eq #conjugates then; 
            Include(~Cat_m_I, a);
        end if;
    end for;

    // Return the number of elements in Cat^m_I
    return Cat_m_I;
end function;
