p := 2;
n := 10; // quadratic APN functions over GF(2^n)
m := 2;  // with coefficients in GF(2^m)

load "funcs.m";
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

