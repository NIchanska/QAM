load "IndSet";

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

function VM(basis)
	res := basis;
	for i in [1 .. n - 1] do
		res cat:= [b^(p^i) : b in basis];
	end for;
	return Matrix(FF, n, n, res);
end function;

function equivSet(set)
    out := [];
    for i in set do
        out:=Append(out, FF!Min(i));
    end for;
    return out;
end function;

function mySplit(list)
        out1 := [list[i] : i in [1..#list] | IsEven(i)];
        out2 := [list[i] : i in [1..#list] | IsOdd(i)];
        return [out1, out2];
end function;

function mySplit4(list)
	out := mySplit(list);
        return mySplit(out[1]) cat mySplit(out[2]);
end function;

function myElt(el);
	res := [[a1,a2,a3,a4,a5,a6,a7,a8] : a1,a2,a3,a4,a5,a6,a7,a8 in GF(2)| &+[[a1,a2,a3,a4,a5,a6,a7,a8][i]*basis[i]: i in [1..n]] eq element where element is el]; 
	return res[1];
end function;

function myVector(el);
	res := [[a1,a2,a3,a4,a5,a6,a7,a8] : a1,a2,a3,a4,a5,a6,a7,a8 in GF(2)| &+[[a1,a2,a3,a4,a5,a6,a7,a8][i]*basis[i]: i in [1..n]] eq element where element is el]; 
	return Matrix(GF(2),n,1,res[1]);
end function;

function myProd(M,cand);
	res := InnerProduct(Vector(basis),Vector(M*myVector(cand)));
	return res;
end function;

function map3(M,x,y)
	return M*myVector(x) eq myVector(y) or M*myVector(y) eq myVector(x);
end function;

function AssignA(M, A)
	for i in [0..n div m - 1] do;
		ind1 := (2+i*m) mod n;
		ind2 := (1+i*m) mod n;

		if ind1 eq 0 then;
			ind1 := n;
		elif ind2 eq 0 then;
			ind2 := n;
		end if;

		M[ind1,ind2] := A^(p^(i*m)); M[ind2,ind1] := M[ind1,ind2];
	end for;
    return M;
end function;

function AssignB(M, B)
	for i in [0..n div m - 1] do;
		ind1 := (3+i*m) mod n;
		ind2 := (1+i*m) mod n;

		if ind1 eq 0 then;
			ind1 := n;
		elif ind2 eq 0 then;
			ind2 := n;
		end if;

		M[ind1,ind2] := B^(p^(i*m)); M[ind2,ind1] := M[ind1,ind2];
	end for;
    return M;
end function;

function AssignC(M, C)
	for i in [0..n div m - 1] do;
		ind1 := (4+i*m) mod n;
		ind2 := (1+i*m) mod n;

		if ind1 eq 0 then;
			ind1 := n;
		elif ind2 eq 0 then;
			ind2 := n;
		end if;

		M[ind1,ind2] := C^(p^(i*m)); M[ind2,ind1] := M[ind1,ind2];
	end for;
    return M;
end function;

function AssignD(M, D)
		ind1 := 5;
		ind2 := 1;
		M[ind1,ind2] := D; M[ind2,ind1] := M[ind1,ind2];
    return M;
end function;

function AssignE(M, E)
	for i in [0..n div m - 1] do;
		ind1 := (6+i*m) mod n;
		ind2 := (1+i*m) mod n;

		if ind1 eq 0 then;
			ind1 := n;
		elif ind2 eq 0 then;
			ind2 := n;
		end if;

		M[ind1,ind2] := E^(p^(i*m)); M[ind2,ind1] := M[ind1,ind2];
	end for;
    return M; 
end function; 

function AssignF(M, F)
	for i in [0..n div m - 1] do;
		ind1 := (7+i*m) mod n;
		ind2 := (1+i*m) mod n;

		if ind1 eq 0 then;
			ind1 := n;
		elif ind2 eq 0 then;
			ind2 := n;
		end if;

		M[ind1,ind2] := F^(p^(i*m)); M[ind2,ind1] := M[ind1,ind2];
	end for;
    return M;
end function;

function AssignG(M, G)
	for i in [0..n div m - 1] do;
		ind1 := (8+i*m) mod n;
		ind2 := (1+i*m) mod n;

		if ind1 eq 0 then;
			ind1 := n;
		elif ind2 eq 0 then;
			ind2 := n;
		end if;

		M[ind1,ind2] := G^(p^(i*m)); M[ind2,ind1] := M[ind1,ind2];
	end for;
    return M;
end function;

function AssignH(M, H)
	for i in [0..n div m - 1] do;
		ind1 := (3+i*m) mod n;
		ind2 := (2+i*m) mod n;

		if ind1 eq 0 then;
			ind1 := n;
		elif ind2 eq 0 then;
			ind2 := n;
		end if;

		M[ind1,ind2] := H^(p^(i*m)); M[ind2,ind1] := M[ind1,ind2];
	end for;
    return M;
end function;

function AssignI(M, I)
	for i in [0..n div m - 1] do;
		ind1 := (4+i*m) mod n;
		ind2 := (2+i*m) mod n;

		if ind1 eq 0 then;
			ind1 := n;
		elif ind2 eq 0 then;
			ind2 := n;
		end if;

		M[ind1,ind2] := I^(p^(i*m)); M[ind2,ind1] := M[ind1,ind2];
	end for;
    return M;
end function;

function AssignJ(M, J)
		ind1 := 6;
		ind2 := 2;
		M[ind1,ind2] := J; M[ind2,ind1] := M[ind1,ind2];
    return M;
end function;

function AssignK(M, K)
	for i in [0..n div m - 1] do;
		ind1 := (7+i*m) mod n;
		ind2 := (2+i*m) mod n;

		if ind1 eq 0 then;
			ind1 := n;
		elif ind2 eq 0 then;
			ind2 := n;
		end if;

		M[ind1,ind2] := K^(p^(i*m)); M[ind2,ind1] := M[ind1,ind2];
	end for;
    return M; 
end function; 

function AssignL(M, L)
	for i in [0..n div m - 1] do;
		ind1 := (8+i*m) mod n;
		ind2 := (2+i*m) mod n;

		if ind1 eq 0 then;
			ind1 := n;
		elif ind2 eq 0 then;
			ind2 := n;
		end if;

		M[ind1,ind2] := L^(p^(i*m)); M[ind2,ind1] := M[ind1,ind2];
	end for;
    return M; 
end function; 

function AssignM(M, MM)
	for i in [0..n div m - 1] do;
		ind1 := (4+i*m) mod n;
		ind2 := (3+i*m) mod n;

		if ind1 eq 0 then;
			ind1 := n;
		elif ind2 eq 0 then;
			ind2 := n;
		end if;

		M[ind1,ind2] := MM^(p^(i*m)); M[ind2,ind1] := M[ind1,ind2];
	end for;
    return M; 
end function; 

function AssignN(M, N)
		ind1 := 7;
		ind2 := 3;
		M[ind1,ind2] := N; M[ind2,ind1] := M[ind1,ind2];
    return M; 
end function; 

function AssignO(M, O)
	for i in [0..n div m - 1] do;
		ind1 := (8+i*m) mod n;
		ind2 := (3+i*m) mod n;

		if ind1 eq 0 then;
			ind1 := n;
		elif ind2 eq 0 then;
			ind2 := n;
		end if;

		M[ind1,ind2] := O^(p^(i*m)); M[ind2,ind1] := M[ind1,ind2];
	end for;
    return M; 
end function; 

function AssignP(M, P)
		ind1 := 4;
		ind2 := 8;
		M[ind1,ind2] := P; M[ind2,ind1] := M[ind1,ind2];
    return M; 
end function; 

function vandemonde(basis, n)
	res := basis;
	for i in [1 .. n - 1] do
		res cat:= [b^(p^i) : b in basis];
	end for;

	return Matrix(FF, n, n, res);
end function;


//Creates a coefficient matrix of a function
function coeffMatrix(F)
    M := ZeroMatrix(FF, n, n);
    for i in [1 .. n] do
        for j in [1 .. n] do
            if i eq j then
                M[i,j] := 2*Coefficient(F, p^(i-1)+p^(j-1));
            else
                M[i,j] := Coefficient(F, p^(i-1)+p^(j-1));
            end if;
        end for;
    end for;
    return M;
end function;

//Creates the derivative matrix given the Vandemonde matrix and coefficient matrix
function derivativeBAB(B, A)
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

//Creates the coefficient matrix given the derivative matrix and Vandemonde matrix
function derivativeToCoefficient(M, B)
    return Transpose(B)^-1 * M * B^-1;
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

function isQAM_F(M)
    rowSize := Nrows(M);
    colSize := Ncols(M);
    V := VectorSpace(FF, colSize);
    
    for i in IndSet do
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

function SubsCleaner(Ss);
	Subs := {};
	for Elem in Ss do
		Include(~Subs,[{ar : ar in Elem[1]},{rr : rr in Elem[2]}]); 
	end for;
	for el1 in Subs do;
		Eqset := {[{(ee + k) mod n: ee in el1[1]},{(ee + k) mod n: ee in el1[2]}]: k in [1..9]}; 
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

function check_subs(M, subs)
    for subInd in subs do
        S := Submatrix(M, subInd[1], subInd[2]);
        if isQAM(S) eq false then
            return false;
        end if;
    end for;
    return true;
end function;

function values2M(A,B,C,D,E,F,G,H,I,J,K,L,MM,N,O,P);
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
	M := AssignJ(M,J);
	M := AssignK(M,K);
	M := AssignL(M,L);
	M := AssignM(M,MM);
	M := AssignN(M,N);
	M := AssignO(M,O);
	M := AssignP(M,P);
	return M;
end function;

function values(tup);
	return tup[1],tup[2],tup[3],tup[4],tup[5],tup[6],tup[7],tup[8],tup[9],tup[10],tup[11],tup[12],tup[13],tup[14],tup[15],tup[16];
end function;

function M2values(M);
	A := M[1,2];
	B := M[1,3];
	C := M[1,4];
	D := M[1,5];
	E := M[1,6];
	F := M[1,7];
	G := M[1,8];
	H := M[2,3];
	I := M[2,4];
	J := M[2,6];
	K := M[2,7];
	L := M[2,8];
	MM := M[3,4];
	N := M[3,7];
	O := M[3,8];
	P := M[4,8];
	return A,B,C,D,E,F,G,H,I,J,K,L,MM,N,O,P;
end function;

function M2Tvalues(M);
	A,B,C,D,E,F,G,H,I,J,K,L,MM,N,O,P := M2values(M);
	return <A,B,C,D,E,F,G,H,I,J,K,L,MM,N,O,P>;
end function;

function is_nm(M);  // checks if diagonalisation property is preserved 
	n := #Rows(M);
	outp := true;
	for i in [1..n] do;
		for j in [1..n] do;
			for k in [0..n-1 by m] do;
				ind1 := (i+k) mod n;
				ind2 := (j+k) mod n;
				if ind1 eq 0 then;
					ind1 := n;
				elif ind2 eq 0 then;
					ind2 := n;
				end if;

				if i ne j and (M[i,j])^(2^k) ne M[ind1,ind2] then;
					outp := false;
					return outp;
				end if;
			end for;
		end for;
	end for;
	return outp;
end function;
