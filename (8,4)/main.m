n := 8;
m := 4;
p := 2;

FF<a> := GF(2^n);
P<y> := PolynomialRing(FF);
Q<x> := quo<P|y^(2^n)+y>;

Fm := [x : x in FF | x in GF(2^m)];
base := NormalElement(FF);
basis := [base^(p^x) : x in [0 .. (n-1)]];

load "funcs.m";
load "SUBS.m";
V := VM(basis);
