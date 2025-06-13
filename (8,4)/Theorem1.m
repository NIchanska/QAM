/* classify if the system has a solution */
/* for cases where m = n div 2 */

function Pows2l(el);
	return [el^(2^i): i in [0..n-1]];
end function;

function Pows2lm(el);
	return [el^(2^i): i in [0..n-1 by m]];
end function;

function equivSet(set)
    out := [];
    for i in set do
        out:=Append(out, FF!Min(i));
    end for;
    return out;
end function;

function Equiv(x,y);
	Xl := Pows2l(x);
	M := Matrix(FF,[Rotate(Xl,-i) : i in [0..n-1 by m]]);
	V := Vector(FF,Pows2lm(y));

	bol, sol := IsConsistent(Transpose(M),V);
	return bol;
end function;

function EquivB(x,y,A);
	Xl := Pows2l(x);
	Al := Pows2l(A); 
	
	Mx := [Rotate(Xl,-i) : i in [0..n-1 by m]];
	MA := [Rotate(Al,-i) : i in [0..n-1 by m]];

	M := Matrix(FF,Mx cat MA);
	V := Vector(FF, Pows2lm(y) cat Pows2lm(A));
	bol, sol := IsConsistent(Transpose(M),V);
	return bol;
end function;

function EquivC(x,y,A,B);
	Xl := Pows2l(x);
	Al := Pows2l(A); 
	Bl := Pows2l(B);

	Mx := [Rotate(Xl,-i) : i in [0..n-1 by m]];
	MA := [Rotate(Al,-i) : i in [0..n-1 by m]];
	MB := [Rotate(Bl,-i) : i in [0..n-1 by m]];

	M := Matrix(FF,Mx cat MA cat MB);
	V := Vector(FF, Pows2lm(y) cat Pows2lm(A) cat Pows2lm(B));
	bol, sol := IsConsistent(Transpose(M),V);
	return bol;
end function;

function EquivD(x,y,A,B,C);
	Xl := Pows2l(x);
	Al := Pows2l(A); 
	Bl := Pows2l(B);
	Cl := Pows2l(C);

	Mx := [Rotate(Xl,-i) : i in [0..n-1 by m]];
	MA := [Rotate(Al,-i) : i in [0..n-1 by m]];
	MB := [Rotate(Bl,-i) : i in [0..n-1 by m]];
	MC := [Rotate(Cl,-i) : i in [0..n-1 by m]];

	M := Matrix(FF,Mx cat MA cat MB cat MC);
	V := Vector(FF, Pows2lm(y) cat Pows2lm(A) cat Pows2lm(B) cat Pows2lm(C));
	bol, sol := IsConsistent(Transpose(M),V);
	return bol;
end function;

function EquivE(x,y,A,B,C,D);
	Xl := Pows2l(x);
	Al := Pows2l(A); 
	Bl := Pows2l(B);
	Cl := Pows2l(C);
	Dl := Pows2l(D);

	Mx := [Rotate(Xl,-i) : i in [0..n-1 by m]];
	MA := [Rotate(Al,-i) : i in [0..n-1 by m]];
	MB := [Rotate(Bl,-i) : i in [0..n-1 by m]];
	MC := [Rotate(Cl,-i) : i in [0..n-1 by m]];
	MD := [Rotate(Dl,-i) : i in [0..n-1 by m]];

	M := Matrix(FF,Mx cat MA cat MB cat MC cat MD);
	V := Vector(FF, Pows2lm(y) cat Pows2lm(A) cat Pows2lm(B) cat Pows2lm(C) cat Pows2lm(D));
	bol, sol := IsConsistent(Transpose(M),V);
	return bol;
end function;

function EquivF(x,y,A,B,C,D,E);
	Xl := Pows2l(x);
	Al := Pows2l(A); 
	Bl := Pows2l(B);
	Cl := Pows2l(C);
	Dl := Pows2l(D);
	El := Pows2l(E);

	Mx := [Rotate(Xl,-i) : i in [0..n-1 by m]];
	MA := [Rotate(Al,-i) : i in [0..n-1 by m]];
	MB := [Rotate(Bl,-i) : i in [0..n-1 by m]];
	MC := [Rotate(Cl,-i) : i in [0..n-1 by m]];
	MD := [Rotate(Dl,-i) : i in [0..n-1 by m]];
	ME := [Rotate(El,-i) : i in [0..n-1 by m]];

	M := Matrix(FF,Mx cat MA cat MB cat MC cat MD cat ME);
	V := Vector(FF, Pows2lm(y) cat Pows2lm(A) cat Pows2lm(B) cat Pows2lm(C) cat Pows2lm(D) cat Pows2lm(E));
	bol, sol := IsConsistent(Transpose(M),V);
	return bol;
end function;

function EquivG(x,y,A,B,C,D,E,F);
	Xl := Pows2l(x);
	Al := Pows2l(A); 
	Bl := Pows2l(B);
	Cl := Pows2l(C);
	Dl := Pows2l(D);
	El := Pows2l(E);
	Fl := Pows2l(F);

	Mx := [Rotate(Xl,-i) : i in [0..n-1 by m]];
	MA := [Rotate(Al,-i) : i in [0..n-1 by m]];
	MB := [Rotate(Bl,-i) : i in [0..n-1 by m]];
	MC := [Rotate(Cl,-i) : i in [0..n-1 by m]];
	MD := [Rotate(Dl,-i) : i in [0..n-1 by m]];
	ME := [Rotate(El,-i) : i in [0..n-1 by m]];
	MF := [Rotate(Fl,-i) : i in [0..n-1 by m]];

	M := Matrix(FF,Mx cat MA cat MB cat MC cat MD cat ME cat MF);
	V := Vector(FF, Pows2lm(y) cat Pows2lm(A) cat Pows2lm(B) cat Pows2lm(C) cat Pows2lm(D) cat Pows2lm(E) cat Pows2lm(F));
	bol, sol := IsConsistent(Transpose(M),V);
	return bol;
end function;

function EquivH(x,y,A,B,C,D,E,F,G);
	Xl := Pows2l(x);
	Al := Pows2l(A); 
	Bl := Pows2l(B);
	Cl := Pows2l(C);
	Dl := Pows2l(D);
	El := Pows2l(E);
	Fl := Pows2l(F);
	Gl := Pows2l(G);

	Mx := [Rotate(Xl,-i) : i in [0..n-1 by m]];
	MA := [Rotate(Al,-i) : i in [0..n-1 by m]];
	MB := [Rotate(Bl,-i) : i in [0..n-1 by m]];
	MC := [Rotate(Cl,-i) : i in [0..n-1 by m]];
	MD := [Rotate(Dl,-i) : i in [0..n-1 by m]];
	ME := [Rotate(El,-i) : i in [0..n-1 by m]];
	MF := [Rotate(Fl,-i) : i in [0..n-1 by m]];
	MG := [Rotate(Gl,-i) : i in [0..n-1 by m]];

	M := Matrix(FF,Mx cat MA cat MB cat MC cat MD cat ME cat MF cat MG);
	V := Vector(FF, Pows2lm(y) cat Pows2lm(A) cat Pows2lm(B) cat Pows2lm(C) cat Pows2lm(D) cat Pows2lm(E) cat Pows2lm(F) cat Pows2lm(G));
	bol, sol := IsConsistent(Transpose(M),V);
	return bol;
end function;

function EquivI(x,y,A,B,C,D,E,F,G,H);
	Xl := Pows2l(x);
	Al := Pows2l(A); 
	Bl := Pows2l(B);
	Cl := Pows2l(C);
	Dl := Pows2l(D);
	El := Pows2l(E);
	Fl := Pows2l(F);
	Gl := Pows2l(G);
	Hl := Pows2l(H);

	Mx := [Rotate(Xl,-i) : i in [0..n-1 by m]];
	MA := [Rotate(Al,-i) : i in [0..n-1 by m]];
	MB := [Rotate(Bl,-i) : i in [0..n-1 by m]];
	MC := [Rotate(Cl,-i) : i in [0..n-1 by m]];
	MD := [Rotate(Dl,-i) : i in [0..n-1 by m]];
	ME := [Rotate(El,-i) : i in [0..n-1 by m]];
	MF := [Rotate(Fl,-i) : i in [0..n-1 by m]];
	MG := [Rotate(Gl,-i) : i in [0..n-1 by m]];
	MH := [Rotate(Hl,-i) : i in [0..n-1 by m]];

	M := Matrix(FF,Mx cat MA cat MB cat MC cat MD cat ME cat MF cat MG cat MH);
	V := Vector(FF, Pows2lm(y) cat Pows2lm(A) cat Pows2lm(B) cat Pows2lm(C) cat Pows2lm(D) cat Pows2lm(E) cat Pows2lm(F) cat Pows2lm(G) cat Pows2lm(H));
	bol, sol := IsConsistent(Transpose(M),V);
	return bol;
end function;

function EquivJ(x,y,A,B,C,D,E,F,G,H,I);
	Xl := Pows2l(x);
	Al := Pows2l(A); 
	Bl := Pows2l(B);
	Cl := Pows2l(C);
	Dl := Pows2l(D);
	El := Pows2l(E);
	Fl := Pows2l(F);
	Gl := Pows2l(G);
	Hl := Pows2l(H);
	Il := Pows2l(I);

	Mx := [Rotate(Xl,-i) : i in [0..n-1 by m]];
	MA := [Rotate(Al,-i) : i in [0..n-1 by m]];
	MB := [Rotate(Bl,-i) : i in [0..n-1 by m]];
	MC := [Rotate(Cl,-i) : i in [0..n-1 by m]];
	MD := [Rotate(Dl,-i) : i in [0..n-1 by m]];
	ME := [Rotate(El,-i) : i in [0..n-1 by m]];
	MF := [Rotate(Fl,-i) : i in [0..n-1 by m]];
	MG := [Rotate(Gl,-i) : i in [0..n-1 by m]];
	MH := [Rotate(Hl,-i) : i in [0..n-1 by m]];
	MI := [Rotate(Il,-i) : i in [0..n-1 by m]];

	M := Matrix(FF,Mx cat MA cat MB cat MC cat MD cat ME cat MF cat MG cat MH cat MI);
	V := Vector(FF, Pows2lm(y) cat Pows2lm(A) cat Pows2lm(B) cat Pows2lm(C) cat Pows2lm(D) cat Pows2lm(E) cat Pows2lm(F) cat Pows2lm(G) cat Pows2lm(H) cat Pows2lm(I));
	bol, sol := IsConsistent(Transpose(M),V);
	return bol;
end function;

function EquivK(x,y,A,B,C,D,E,F,G,H,I,J);
	Xl := Pows2l(x);
	Al := Pows2l(A); 
	Bl := Pows2l(B);
	Cl := Pows2l(C);
	Dl := Pows2l(D);
	El := Pows2l(E);
	Fl := Pows2l(F);
	Gl := Pows2l(G);
	Hl := Pows2l(H);
	Il := Pows2l(I);
	Jl := Pows2l(J);

	Mx := [Rotate(Xl,-i) : i in [0..n-1 by m]];
	MA := [Rotate(Al,-i) : i in [0..n-1 by m]];
	MB := [Rotate(Bl,-i) : i in [0..n-1 by m]];
	MC := [Rotate(Cl,-i) : i in [0..n-1 by m]];
	MD := [Rotate(Dl,-i) : i in [0..n-1 by m]];
	ME := [Rotate(El,-i) : i in [0..n-1 by m]];
	MF := [Rotate(Fl,-i) : i in [0..n-1 by m]];
	MG := [Rotate(Gl,-i) : i in [0..n-1 by m]];
	MH := [Rotate(Hl,-i) : i in [0..n-1 by m]];
	MI := [Rotate(Il,-i) : i in [0..n-1 by m]];
	MJ := [Rotate(Jl,-i) : i in [0..n-1 by m]];

	M := Matrix(FF,Mx cat MA cat MB cat MC cat MD cat ME cat MF cat MG cat MH cat MI cat MJ);
	V := Vector(FF, Pows2lm(y) cat Pows2lm(A) cat Pows2lm(B) cat Pows2lm(C) cat Pows2lm(D) cat Pows2lm(E) cat Pows2lm(F) cat Pows2lm(G) cat Pows2lm(H) cat Pows2lm(I) cat Pows2lm(J));
	bol, sol := IsConsistent(Transpose(M),V);
	return bol;
end function;

function EquivL(x,y,A,B,C,D,E,F,G,H,I,J,K);
	Xl := Pows2l(x);
	Al := Pows2l(A); 
	Bl := Pows2l(B);
	Cl := Pows2l(C);
	Dl := Pows2l(D);
	El := Pows2l(E);
	Fl := Pows2l(F);
	Gl := Pows2l(G);
	Hl := Pows2l(H);
	Il := Pows2l(I);
	Jl := Pows2l(J);
	Kl := Pows2l(K);

	Mx := [Rotate(Xl,-i) : i in [0..n-1 by m]];
	MA := [Rotate(Al,-i) : i in [0..n-1 by m]];
	MB := [Rotate(Bl,-i) : i in [0..n-1 by m]];
	MC := [Rotate(Cl,-i) : i in [0..n-1 by m]];
	MD := [Rotate(Dl,-i) : i in [0..n-1 by m]];
	ME := [Rotate(El,-i) : i in [0..n-1 by m]];
	MF := [Rotate(Fl,-i) : i in [0..n-1 by m]];
	MG := [Rotate(Gl,-i) : i in [0..n-1 by m]];
	MH := [Rotate(Hl,-i) : i in [0..n-1 by m]];
	MI := [Rotate(Il,-i) : i in [0..n-1 by m]];
	MJ := [Rotate(Jl,-i) : i in [0..n-1 by m]];
	MK := [Rotate(Kl,-i) : i in [0..n-1 by m]];

	M := Matrix(FF,Mx cat MA cat MB cat MC cat MD cat ME cat MF cat MG cat MH cat MI cat MJ cat MK);
	V := Vector(FF, Pows2lm(y) cat Pows2lm(A) cat Pows2lm(B) cat Pows2lm(C) cat Pows2lm(D) cat Pows2lm(E) cat Pows2lm(F) cat Pows2lm(G) cat Pows2lm(H) cat Pows2lm(I) cat Pows2lm(J) cat Pows2lm(K));
	bol, sol := IsConsistent(Transpose(M),V);
	return bol;
end function;

function Parti(REM);
	if #REM eq 0 then;
		return {};
	end if;
	x := Random(REM);
	Orb := {x};
	for y in REM diff {x} do
		if Equiv(x,y) then; // here is important Equiv to choose, depending on the depth of the value at the graph 
			Include(~Orb,y);
		else;
			continue;
		end if;
	end for;
	return Orb;
end function;

function LinClassesSol(Os,REM1,REM2);
	//"starts";
	Orb1 := Parti(REM1);
	Orb2 := Parti(REM2);
	Include(~Os,Orb1);
	Include(~Os,Orb2);
	REM1 := REM1 diff Orb1;
	REM2 := REM2 diff Orb2;
	return Os,REM1,REM2;
end function;

function ClassifySol();
	Os := {};
	REM := {x : x in FF | x ne 0};
	REM1 := {x : x in REM | x + x^(2^m) eq 0};
	REM2 := REM diff REM1; 
	while #(REM1 join REM2 ) ne 0 do;
		Os,REM1,REM2:= LinClassesSol(Os,REM1,REM2);
	end while;
	Os := Exclude(Os,{});
	Bs := equivSet(Os);
	return Bs;
end function;

function PartiB_(REM,A);
	if #REM eq 0 then;
		return {};
	end if;
	x := Random(REM);
	Orb := {x};
	for y in REM diff {x} do
		if EquivB(x,y,A) then; // here is important Equiv to choose, depending on the depth of the value at the graph 
			Include(~Orb,y);
		else;
			continue;
		end if;
	end for;
	return Orb;
end function;

function PartiB(REM,A);
	if #REM eq 0 then;
		return {};
	end if;
	x := Random(REM);
	Orb := {x};
	for y in REM diff {x} do
		if EquivB(x,y,A) and EquivB(y,x,A) then; // here is important Equiv to choose, depending on the depth of the value at the graph 
			Include(~Orb,y);
		else;
			continue;
		end if;
	end for;
	return Orb;
end function;

function LinClassesSolB(Os,REM1,REM2,A);
	//"starts";
	Orb1 := PartiB(REM1,A);
	Orb2 := PartiB(REM2,A);
	Include(~Os,Orb1);
	Include(~Os,Orb2);
	REM1 := REM1 diff Orb1;
	REM2 := REM2 diff Orb2;
	return Os,REM1,REM2;
end function;

function ClassifySolB(A);
	Os := {};
	REM := {x : x in FF | x ne 0};
	REM1 := {x : x in REM | x + x^(2^m) eq 0};
	REM2 := REM diff REM1; 
	while #(REM1 join REM2 ) ne 0 do;
		Os,REM1,REM2:= LinClassesSolB(Os,REM1,REM2,A);
	end while;
	Os := Exclude(Os,{});
	Bs := equivSet(Os);
	return Bs;
end function;

function PartiC(REM,A,B);
	if #REM eq 0 then;
		return {};
	end if;
	x := Random(REM);
	Orb := {x};
	for y in REM diff {x} do
		if EquivC(x,y,A,B) and EquivC(y,x,A,B) then; // here is important Equiv to choose, depending on the depth of the value at the graph 
			Include(~Orb,y);
		else;
			continue;
		end if;
	end for;
	return Orb;
end function;

function LinClassesSolC(Os,REM1,REM2,A,B);
	//"starts";
	Orb1 := PartiC(REM1,A,B);
	Orb2 := PartiC(REM2,A,B);
	Include(~Os,Orb1);
	Include(~Os,Orb2);
	REM1 := REM1 diff Orb1;
	REM2 := REM2 diff Orb2;
	return Os,REM1,REM2;
end function;

function ClassifySolC(A,B);
	Os := {};
	REM := {x : x in FF | x ne 0};
	REM1 := {x : x in REM | x + x^(2^m) eq 0};
	REM2 := REM diff REM1;
	while #(REM1 join REM2) ne 0 do;
		Os,REM1,REM2:= LinClassesSolC(Os,REM1,REM2,A,B);
	end while;
	Os := Exclude(Os,{});
	Cs := equivSet(Os);
	return Cs;
end function;

function PartiD(REM,A,B,C);
	if #REM eq 0 then;
		return {};
	end if;
	x := Random(REM);
	Orb := {x};
	for y in REM diff {x} do
		if EquivD(x,y,A,B,C) and EquivD(y,x,A,B,C) then; // here is important Equiv to choose, depending on the depth of the value at the graph 
			Include(~Orb,y);
		else;
			continue;
		end if;
	end for;
	return Orb;
end function;

function LinClassesSolD(Os,REM1,REM2,A,B,C);
	//"starts";
	Orb1 := PartiD(REM1,A,B,C);
	Orb2 := PartiD(REM2,A,B,C);
	Include(~Os,Orb1);
	Include(~Os,Orb2);
	REM1 := REM1 diff Orb1;
	REM2 := REM2 diff Orb2;
	return Os,REM1,REM2;
end function;

function ClassifySolD(A,B,C);
	Os := {};
	REM := {x : x in FF | x ne 0};
	REM1 := {x : x in REM | x + x^(2^m) eq 0};
	REM2 := REM diff REM1; 
	while #(REM1 join REM2) ne 0 do;
		Os,REM1,REM2:= LinClassesSolD(Os,REM1,REM2,A,B,C);
	end while;
	Os := Exclude(Os,{});
	Ds := equivSet(Os);
	return Ds;
end function;

function PartiE(REM,A,B,C,D);
	if #REM eq 0 then;
		return {};
	end if;
	x := Random(REM);
	Orb := {x};
	for y in REM diff {x} do
		if EquivE(x,y,A,B,C,D) and EquivE(y,x,A,B,C,D) then;
			Include(~Orb,y);
		else;
			continue;
		end if;
	end for;
	return Orb;
end function;

function LinClassesSolE(Os,REM1,REM2,A,B,C,D);
	//"starts";
	Orb1 := PartiE(REM1,A,B,C,D);
	Orb2 := PartiE(REM2,A,B,C,D);
	Include(~Os,Orb1);
	Include(~Os,Orb2);
	REM1 := REM1 diff Orb1;
	REM2 := REM2 diff Orb2;
	return Os,REM1,REM2;
end function;

function ClassifySolE(A,B,C,D);
	Os := {};
	REM := {x : x in FF | x ne 0};
	REM1 := {x : x in REM | x + x^(2^m) eq 0};
	REM2 := REM diff REM1; 
	while #(REM1 join REM2) ne 0 do;
		Os,REM1,REM2:= LinClassesSolE(Os,REM1,REM2,A,B,C,D);
	end while;
	Os := Exclude(Os,{});
	Es := equivSet(Os);
	return Es;
end function;

function PartiF(REM,A,B,C,D,E);
	if #REM eq 0 then;
		return {};
	end if;
	x := Random(REM);
	Orb := {x};
	for y in REM diff {x} do
		if EquivF(x,y,A,B,C,D,E) and EquivF(y,x,A,B,C,D,E) then;
			Include(~Orb,y);
		else;
			continue;
		end if;
	end for;
	return Orb;
end function;

function LinClassesSolF(Os,REM1,REM2,A,B,C,D,E);
	//"starts";
	Orb1 := PartiF(REM1,A,B,C,D,E);
	Orb2 := PartiF(REM2,A,B,C,D,E);
	Include(~Os,Orb1);
	Include(~Os,Orb2);
	REM1 := REM1 diff Orb1;
	REM2 := REM2 diff Orb2;
	return Os,REM1,REM2;
end function;

function ClassifySolF(A,B,C,D,E);
	Os := {};
	REM := {x : x in FF | x ne 0};
	REM1 := {x : x in REM | x + x^(2^m) eq 0};
	REM2 := REM diff REM1; 
	while #(REM1 join REM2) ne 0 do;
		Os,REM1,REM2:= LinClassesSolF(Os,REM1,REM2,A,B,C,D,E);
	end while;
	Os := Exclude(Os,{});
	Fs := equivSet(Os);
	return Fs;
end function;

function PartiG(REM,A,B,C,D,E,F);
	if #REM eq 0 then;
		return {};
	end if;
	x := Random(REM);
	Orb := {x};
	for y in REM diff {x} do
		if EquivG(x,y,A,B,C,D,E,F) and EquivG(y,x,A,B,C,D,E,F) then;
			Include(~Orb,y);
		else;
			continue;
		end if;
	end for;
	return Orb;
end function;

function LinClassesSolG(Os,REM1,REM2,A,B,C,D,E,F);
	//"starts";
	Orb1 := PartiG(REM1,A,B,C,D,E,F);
	Orb2 := PartiG(REM2,A,B,C,D,E,F);
	Include(~Os,Orb1);
	Include(~Os,Orb2);
	REM1 := REM1 diff Orb1;
	REM2 := REM2 diff Orb2;
	return Os,REM1,REM2;
end function;

function ClassifySolG(A,B,C,D,E,F);
	Os := {};
	REM := {x : x in FF | x ne 0};
	REM1 := {x : x in REM | x + x^(2^m) eq 0};
	REM2 := REM diff REM1; 
	while #(REM1 join REM2) ne 0 do;
		Os,REM1,REM2:= LinClassesSolG(Os,REM1,REM2,A,B,C,D,E,F);
	end while;
	Os := Exclude(Os,{});
	Gs := equivSet(Os);
	return Gs;
end function;

function PartiH(REM,A,B,C,D,E,F,G);
	if #REM eq 0 then;
		return {};
	end if;
	x := Random(REM);
	Orb := {x};
	for y in REM diff {x} do
		if EquivH(x,y,A,B,C,D,E,F,G) and EquivH(y,x,A,B,C,D,E,F,G) then;
			Include(~Orb,y);
		else;
			continue;
		end if;
	end for;
	return Orb;
end function;

function LinClassesSolH(Os,REM1,REM2,A,B,C,D,E,F,G);
	//"starts";
	Orb1 := PartiH(REM1,A,B,C,D,E,F,G);
	Orb2 := PartiH(REM2,A,B,C,D,E,F,G);
	Include(~Os,Orb1);
	Include(~Os,Orb2);
	REM1 := REM1 diff Orb1;
	REM2 := REM2 diff Orb2;
	return Os,REM1,REM2;
end function;

function ClassifySolH(A,B,C,D,E,F,G);
	Os := {};
	REM := {x : x in FF | x ne 0};
	REM1 := {x : x in REM | x + x^(2^m) eq 0};
	REM2 := REM diff REM1;
	while #(REM1 join REM2) ne 0 do;
		Os,REM1,REM2:= LinClassesSolH(Os,REM1,REM2,A,B,C,D,E,F,G);
	end while;
	Os := Exclude(Os,{});
	Hs := equivSet(Os);
	return Hs;
end function;

function PartiI(REM,A,B,C,D,E,F,G,H);
	if #REM eq 0 then;
		return {};
	end if;
	x := Random(REM);
	Orb := {x};
	for y in REM diff {x} do
		if EquivI(x,y,A,B,C,D,E,F,G,H) and EquivI(y,x,A,B,C,D,E,F,G,H) then;
			Include(~Orb,y);
		else;
			continue;
		end if;
	end for;
	return Orb;
end function;

function LinClassesSolI(Os,REM1,REM2,A,B,C,D,E,F,G,H);
	//"starts";
	Orb1 := PartiI(REM1,A,B,C,D,E,F,G,H);
	Orb2 := PartiI(REM2,A,B,C,D,E,F,G,H);
	Include(~Os,Orb1);
	Include(~Os,Orb2);
	REM1 := REM1 diff Orb1;
	REM2 := REM2 diff Orb2;
	return Os,REM1,REM2;
end function;

function ClassifySolI(A,B,C,D,E,F,G,H);
	Os := {};
	REM := {x : x in FF | x ne 0};
	REM1 := {x : x in REM | x + x^(2^m) eq 0};
	REM2 := REM diff REM1;
	while #(REM1 join REM2) ne 0 do;
		Os,REM1,REM2:= LinClassesSolI(Os,REM1,REM2,A,B,C,D,E,F,G,H);
	end while;
	Os := Exclude(Os,{});
	Is := equivSet(Os);
	return Is;
end function;

function PartiJ(REM,A,B,C,D,E,F,G,H,I);
	if #REM eq 0 then;
		return {};
	end if;
	x := Random(REM);
	Orb := {x};
	for y in REM diff {x} do
		if EquivJ(x,y,A,B,C,D,E,F,G,H,I) and EquivJ(y,x,A,B,C,D,E,F,G,H,I) then;
			Include(~Orb,y);
		else;
			continue;
		end if;
	end for;
	return Orb;
end function;

function LinClassesSolJ(Os,REM1,REM2,A,B,C,D,E,F,G,H,I);
	//"starts";
	Orb1 := PartiJ(REM1,A,B,C,D,E,F,G,H,I);
	Orb2 := PartiJ(REM2,A,B,C,D,E,F,G,H,I);
	Include(~Os,Orb1);
	Include(~Os,Orb2);
	REM1 := REM1 diff Orb1;
	REM2 := REM2 diff Orb2;
	return Os,REM1,REM2;
end function;

function ClassifySolJ(A,B,C,D,E,F,G,H,I);
	Os := {};
	REM := {x : x in FF | x ne 0};
	REM1 := {x : x in REM | x + x^(2^m) eq 0};
	REM2 := REM diff REM1;
	while #(REM1 join REM2) ne 0 do;
		Os,REM1,REM2:= LinClassesSolJ(Os,REM1,REM2,A,B,C,D,E,F,G,H,I);
	end while;
	Os := Exclude(Os,{});
	Js := equivSet(Os);
	return Js;
end function;

function PartiK(REM,A,B,C,D,E,F,G,H,I,J);
	if #REM eq 0 then;
		return {};
	end if;
	x := Random(REM);
	Orb := {x};
	for y in REM diff {x} do
		if EquivK(x,y,A,B,C,D,E,F,G,H,I,J) and EquivK(y,x,A,B,C,D,E,F,G,H,I,J) then;
			Include(~Orb,y);
		else;
			continue;
		end if;
	end for;
	return Orb;
end function;

function LinClassesSolK(Os,REM1,REM2,A,B,C,D,E,F,G,H,I,J);
	//"starts";
	Orb1 := PartiK(REM1,A,B,C,D,E,F,G,H,I,J);
	Orb2 := PartiK(REM2,A,B,C,D,E,F,G,H,I,J);
	Include(~Os,Orb1);
	Include(~Os,Orb2);
	REM1 := REM1 diff Orb1;
	REM2 := REM2 diff Orb2;
	return Os,REM1,REM2;
end function;

function ClassifySolK(A,B,C,D,E,F,G,H,I,J);
	Os := {};
	REM := {x : x in FF | x ne 0};
	REM1 := {x : x in REM | x + x^(2^m) eq 0};
	REM2 := REM diff REM1;
	while #(REM1 join REM2) ne 0 do;
		Os,REM1,REM2:= LinClassesSolK(Os,REM1,REM2,A,B,C,D,E,F,G,H,I,J);
	end while;
	Os := Exclude(Os,{});
	Ks := equivSet(Os);
	return Ks;
end function;

function PartiL(REM,A,B,C,D,E,F,G,H,I,J,K);
	if #REM eq 0 then;
		return {};
	end if;
	x := Random(REM);
	Orb := {x};
	for y in REM diff {x} do
		if EquivL(x,y,A,B,C,D,E,F,G,H,I,J,K) and EquivL(y,x,A,B,C,D,E,F,G,H,I,J,K) then;
			Include(~Orb,y);
		else;
			continue;
		end if;
	end for;
	return Orb;
end function;

function LinClassesSolL(Os,REM1,REM2,A,B,C,D,E,F,G,H,I,J,K);
	//"starts";
	Orb1 := PartiL(REM1,A,B,C,D,E,F,G,H,I,J,K);
	Orb2 := PartiL(REM2,A,B,C,D,E,F,G,H,I,J,K);
	Include(~Os,Orb1);
	Include(~Os,Orb2);
	REM1 := REM1 diff Orb1;
	REM2 := REM2 diff Orb2;
	return Os,REM1,REM2;
end function;

function ClassifySolL(A,B,C,D,E,F,G,H,I,J,K);
	Os := {};
	REM := {x : x in FF | x ne 0};
	REM1 := {x : x in REM | x + x^(2^m) eq 0};
	REM2 := REM diff REM1;
	while #(REM1 join REM2) ne 0 do;
		Os,REM1,REM2:= LinClassesSolL(Os,REM1,REM2,A,B,C,D,E,F,G,H,I,J,K);
	end while;
	Os := Exclude(Os,{});
	Ls := equivSet(Os);
	return Ls;
end function;


function Category1(n, m)
    // Define the field F_2^n and subfield F_2^m
    F := GF(2^n);
    F_sub := GF(2^m);

    // Initialize two empty sets for independent (Cat^m_I) and dependent (Cat^m_D) categories
    Cat_Ind := {};

    // Loop over all non-zero elements in F
    for a in F do
        // Skip 0 as it has no conjugates
        if a eq 0 then
            continue;
        end if;

        // Collect the conjugates of a in F with respect to F_sub
        conjugates := [a^(2^(m*k)) : k in [0..n div m - 1]];

        // Check if the conjugates are linearly independent
        V := VectorSpace(F, #conjugates);
        conj_vectors := [V!Eltseq(conjugates[i]) : i in [1..#conjugates]];

        // If the conjugates are linearly independent, include a in Cat^m_I
        if Rank(Matrix(conj_vectors)) eq #conjugates then
            Include(~Cat_Ind, a);
	   end if;
   end for;

   return Cat_Ind
end function;
