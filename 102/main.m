p := 2;
n := 10; // quadratic APN functions over GF(2^n)
m := 2;  // with coefficients in GF(2^m)

load "funcs.m";
//load "~/dermat/GenSearch/utils.m";
load "SUBS.m"; 

Fm := [x : x in FF | x in GF(2^m)];

BHK := x^288 + a^682*x^96 + a^341*x^9 + x^3;

base := NormalElement(FF);
basis := [base^(p^x) : x in [0 .. (n-1)]];

V := VM(basis);

/* prep part

Start := Realtime();
lin_perm := { &+[ [c0,c1,c2,c3,c4,c5,c6,c7,c8,c9][i+1]*x^(2^i) : i in [0..n-1]] :  c0,c1,c2,c3,c4,c5,c6,c7,c8,c9 in Fm};
"O(n) of finding all linear functions:", Realtime()-Start;

perm_classes := findLinear(lin_perm);
"O(n) of finding all equivalence classes generated from permutations:", Realtime()-Start;
Write("lin_perm", perm_classes);

*/

//As := baseClass(lin_perm);
/* all info about partitioning*/

//Bs := [];
//Cs := [];
//Ds := [];
//Es := [];
//Fs := [];
//Gs := [];
//Hs := [];
//Is := [];

/*
for A in As do
	Blist, L1 := classUni(A,lin_perm);
	Append(~Bs, <A,Blist,L1>);
	Append(~B_orb,<A,Blist>);
end for;

q
Write("orb_res","Bs :=");
Write("orb_res",Bs);
Write("orb_res",";");

"We're done with Bs";

for B in Bs do
	A := B[1];
	Blist := B[2];
	L1 := B[3];
	for bo in Blist do
		Clist,L2 := classC(bo,L1);
		Append(~Cs, <A,bo,Clist,L2>);
	end for;
end for;

Write("orb_res","Cs :=");
Write("orb_res",Cs);
Write("orb_res",";");

"We're done with Cs";

for C in Cs do
	A := C[1];
	B := C[2];
	Clist := C[3];
	L2 := C[4]; 
        for co in Clist do
	        Dlist,L3 := classD(A,B,co,L2);
              	Append(~Ds, <A,B,co,Dlist,L3>);
	end for;
end for;

Write("orb_res","Ds :=");
Write("orb_res",Ds);
Write("orb_res",";");

"We're done with Ds";

for D in Ds do
        A := D[1];
        B := D[2];
	C := D[3];
        Dlist := D[4];
        L3 := D[5];
        for d_o in Dlist do
                Elist,L4 := classE(d_o,L3);
                Append(~Es, <A,B,C,d_o,Elist,L4>);
        end for;
end for;

Write("orb_res","Es :=");
Write("orb_res",Es);
Write("orb_res",";");

"We're done with Es";

for E in Es do
        A := E[1];
        B := E[2];
        C := E[3];
	D := E[4];
        Elist := E[5];
        L4 := E[6];
        for eo in Elist do
		Flist,L5 := classF(eo,L4);
		Append(~Fs, <A,B,C,D,eo,Flist,L5>);
        end for;
end for;

Write("orb_res","Fs :=");
Write("orb_res",Fs);
Write("orb_res",";");

"We're done with Fs";
		

for F in Fs do
        A := F[1];
        B := F[2];
        C := F[3];
        D := F[4];
	E := F[5];
        Flist := F[6];
        L5 := F[7];
        for fo in Flist do
		Glist,L6 := classG(fo,L5);
		Append(~Gs, <A,B,C,D,E,fo,Glist,L6>);
        end for;
end for;

Write("orb_res","Gs :=");
Write("orb_res",Gs);
Write("orb_res",";");

"We're done with Gs";

for G in Gs do
        A := G[1];
        B := G[2];
        C := G[3];
        D := G[4];
        E := G[5];
	F := G[6];
        Glist := G[7];
        L6 := G[8];
        for go in Glist do
		Hlist,L7 := classH(go,L6);
		Append(~Hs, <A,B,C,D,E,F,go,Hlist,L7>);
        end for;
end for;

Write("orb_res","Hs :=");
Write("orb_res",Hs);
Write("orb_res",";");

"We're done with Hs";

for H in Hs do
        A := H[1];
        B := H[2];
        C := H[3];
        D := H[4];
        E := H[5];
        F := H[6];
	G := H[7];
        Hlist := H[8];
        L7 := H[9];
        for ho in Hlist do
		Ilist,L8 := classI(ho,L7);
		Append(~Is, <A,B,C,D,E,F,G,ho,Ilist,L8>);
        end for;
end for;

Write("orb_res","Is :=");
Write("orb_res",Is);
Write("orb_res",";");

"We're done with Is";
*/
