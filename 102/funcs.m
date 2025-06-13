FF<a> := GF(p^n);
R<y> := PolynomialRing(FF);
Q<x> := quo<R | y^#FF-y>;

function VM(basis)
	res := basis;
	for i in [1 .. n - 1] do
		res cat:= [b^(p^i) : b in basis];
	end for;
	return Matrix(FF, n, n, res);
end function;

//Creates a coefficient matrix of a function
function coefA(F)
	A := ZeroMatrix(FF, n, n); 
	for i in [1 .. n] do
        	for j in [1 .. n] do
			if i eq j then
                		A[i,j] := 2*Coefficient(F, p^(i-1)+p^(j-1));
            		else
                		A[i,j] := Coefficient(F, p^(i-1)+p^(j-1));
            		end if;
        	end for;
    	end for;
    	return A;
end function;

//Creates the derivative matrix given the Vandemonde matrix and coefficient matrix
function derB(B, A)
    return Transpose(B)*A*B;
end function;

//Creates the derivative matrix by taking every derivative
function derivativeMatrix(F, basis)
        M := ZeroMatrix(FF, n, n);
        for i in [1 .. n] do
        for j in [1 .. n] do
            M[i,j] := Evaluate(F, basis[i] + basis[j]) - Evaluate(F, basis[i]) - Evaluate(F, basis[j]);
                end for;
        end for;
        return M;
end function;	

//  Creates a univariate polynomial from a coefficient matrix:
function coefficientToPolynomial(A)
    poly := 0;
    dim := Ncols(A);
    for i in [1 .. dim] do
        for j in [i .. dim] do
            poly := poly + A[i,j]*x^(p^(i-1) + p^(j-1));
        end for;
    end for;
    return poly;
end function;

//Creates the coefficient matrix given the derivative matrix and Vandemonde matrix
function derivativeToCoefficient(M, B)
    return Transpose(B)^-1 * M * B^-1;
end function;

//Given a row (b1, b2, .. , bn) as a vector, returns the rank of said row.
function APNflat(vec)
    len := Ncols(vec);
    M := ZeroMatrix(FF, len, n);

    for i in [1..len] do
        tmp := Flat(vec[i]);
        for j in [1..n] do
            M[i, j] := tmp[j];
        end for;
    end for;
    rank := Rank(M);
    if rank ge (len - 1) then
        //"\n";
        //M;
        return true;
    else
        //"\n";
        //M;
        return false;
    end if;
end function;


//Testing if a derivative matrix is QAM
function isQAM(M)
    rowSize := Nrows(M);
    colSize := Ncols(M);
    V := VectorSpace(FF, colSize);

    for i in Subsets({1 .. rowSize}) do
        if #i eq 0 then
            continue;
        end if;
        //i;
        vec := Zero(V);
        for j in i do
            vec := vec + M[j];
        end for;

        if APNflat(vec) eq false then
            //"\n";
            //vec;
            return false;
        end if;
    end for;
    //"\n";
    //M;
    return true;
end function;

/*Creating derivative matrix by assigning unknown values A,B,C,D,E,F,G,H,I for dim10*/

function AssignA(M, A)
    M[2,1] := A; M[1,2] := M[2,1];
    M[4,3] := A^(p^2); M[3,4] := M[4,3];
    M[6,5] := A^(p^4); M[5,6] := M[6,5];
    M[8,7] := A^(p^6); M[7,8] := M[8,7];
    M[10,9] := A^(p^8); M[9,10] := M[10,9];
    return M;
end function;

function AssignB(M, B)
    M[3,1] := B; M[1,3] := M[3,1];
    M[5,3] := B^(p^2); M[3,5] := M[5,3];
    M[7,5] := B^(p^4); M[5,7] := M[7,5];
    M[9,7] := B^(p^6); M[7,9] := M[9,7];
    M[1,9] := B^(p^8); M[9,1] := M[1,9];
    return M;
end function;

function AssignC(M, C)
    M[4,1] := C; M[1,4] := M[4,1];
    M[6,3] := C^(p^2); M[3,6] := M[6,3];
    M[8,5] := C^(p^4); M[5,8] := M[8,5];
    M[10,7] := C^(p^6); M[7,10] := M[10,7];
    M[2,9] := C^(p^8); M[9,2] := M[2,9];
    return M;
end function;

function AssignD(M, D)
    M[5,1] := D; M[1,5] := M[5,1];
    M[7,3] := D^(p^2); M[3,7] := M[7,3];
    M[9,5] := D^(p^4); M[5,9] := M[9,5];
    M[1,7] := D^(p^6); M[7,1] := M[1,7];
    M[3,9] := D^(p^8); M[9,3] := M[3,9];
    return M;
end function;

function AssignE(M, E)
    M[6,1] := E; M[1,6] := M[6,1];
    M[3,8] := E^(p^2); M[8,3] := M[3,8];
    M[5,10] := E^(p^4); M[10,5] := M[5,10];
    M[7,2] := E^(p^6); M[2,7] := M[7,2];
    M[9,4] := E^(p^8); M[4,9] := M[9,4];
    return M; 
end function; 
    
function AssignF(M, F) 
    M[8,1] := F; M[1,8] := M[8,1];
    M[10,3] := F^(p^2); M[3,10] := M[10,3];
    M[2,5] := F^(p^4); M[5,2] := M[2,5];
    M[4,7] := F^(p^6); M[7,4] := M[4,7];
    M[6,9] := F^(p^8); M[9,6] := M[6,9];
    return M;
end function;

function AssignG(M, G)
    M[10,1] := G; M[1,10] := M[10,1];
    M[2,3] := G^(p^2); M[3,2] := M[2,3];
    M[4,5] := G^(p^4); M[5,4] := M[4,5];
    M[6,7] := G^(p^6); M[7,6] := M[6,7];
    M[8,9] := G^(p^8); M[9,8] := M[8,9];
    return M;
end function;

function AssignH(M, H)
    M[4,2] := H; M[2,4] := M[4,2];
    M[6,4] := H^(p^2); M[4,6] := M[6,4];
    M[8,6] := H^(p^4); M[6,8] := M[8,6];
    M[10,8] := H^(p^6); M[8,10] := M[10,8];
    M[2,10] := H^(p^8); M[10,2] := M[2,10];
    return M;
end function;

function AssignI(M, I)
    M[6,2] := I; M[2,6] := M[6,2];
    M[8,4] := I^(p^2); M[4,8] := M[8,4];
    M[10,6] := I^(p^4); M[6,10] := M[10,6];
    M[2,8] := I^(p^6); M[8,2] := M[2,8];
    M[4,10] := I^(p^8); M[10,4] := M[4,10];
    return M;
end function;


function mySplit(list)
	out1 := [list[i] : i in [1..#list] | IsEven(i)];
	out2 := [list[i] : i in [1..#list] | IsOdd(i)];
	return out1, out2;
end function;

function values2M(A,B,C,D,E,F,G,H,I);
	M := ZeroMatrix(FF,n,n);
	M := AssignA(M,A);
	M := AssignB(M,B);
	M := AssignC(M,C);
	M := AssignD(M,D);
	M := AssignE(M,E);
	M := AssignF(M,F);
	M := AssignG(M,G);
	M := AssignH(M,H);
	M := AssignI(M,I);
	return M;
end function;

function values(A,B,C,D,E,F,G,H,I);
	return A,B,C,D,E,F,G,H,I;
end function;
