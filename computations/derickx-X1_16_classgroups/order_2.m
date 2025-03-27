// Case 1.a;
_<u>:=PolynomialRing(Rationals());
C1:=HyperellipticCurve(u^4+1);
C2:=HyperellipticCurve(2*(u^4+1));
E1,phi1 := EllipticCurve(C1, Points(C1 : Bound:=100)[1]);
E2,phi2 := EllipticCurve(C2, Points(C2 : Bound:=100)[1]);
r1, flag1 := Rank(E1);
r2, flag2 := Rank(E2);

// Claim 1: E1 and E2 have rank 0
assert  r1 eq 0 and flag1;
assert  r2 eq 0 and flag2;
print "Claim 1 succesfully verified!";

T1, mu1 := TorsionSubgroup(E1);
T2, mu2 := TorsionSubgroup(E2);

// Claim 2: E1 and E2 have (Z/2Z)^2 as torsion subgroups
assert Invariants(T1) eq [ 2, 2 ];
assert Invariants(T2) eq [ 2, 2 ];
print "Claim 2 succesfully verified!";

print "Claim 3: (u:y:v)-coordinates of the torsion points";
[Inverse(phi1)(mu1(t)) : t in T1];
[Inverse(phi2)(mu2(t)) : t in T2];


// Case 1.b
C3:=HyperellipticCurve(u^4+2*u^2-1);
C4:=HyperellipticCurve(2*(u^4+2*u^2-1));
C5:=HyperellipticCurve(u^4-2*u^2-1);
C6:=HyperellipticCurve(2*(u^4-2*u^2-1));

// Claim 4: C6 has no points over Z/16Z
assert {y^2 mod 16 : y in [0..15]} eq { 0, 1, 4, 9 };
assert {2*(u^4-2*u^2*v^2-v^4) mod 16 : u,v in [0..15] | [u mod 2,v mod 2] ne [0,0]} eq { 2, 12, 14 };
print "Claim 4 succesfully verified!";

E3,phi3 := EllipticCurve(C3, Points(C3 : Bound:=100)[1]);
E4,phi4 := EllipticCurve(C4, Points(C4 : Bound:=100)[1]);
E5,phi5 := EllipticCurve(C5, Points(C5 : Bound:=100)[1]);


r3, flag3 := Rank(E3);
r4, flag4 := Rank(E4);
r5, flag5 := Rank(E5);


// Claim 5: E3 and E4 and E5 have ranks 0, 1 and 0
assert  r3 eq 0 and flag3;
assert  r4 eq 1 and flag4;
assert  r5 eq 0 and flag5;
print "Claim 5 succesfully verified!";

T3, mu3 := TorsionSubgroup(E3);
T4, mu4 := TorsionSubgroup(E4);
T5, mu5 := TorsionSubgroup(E5);

// Claim 6: E3 and E4 and E5 have (Z/2Z)^2 as torsion subgroups
assert Invariants(T3) eq [ 2 ];
assert Invariants(T4) eq [ 2 ];
assert Invariants(T5) eq [ 2 ];
print "Claim 6 succesfully verified!";

print "Claim 7: (u:y:v)-coordinates of the torsion points";
[Inverse(phi3)(mu3(t)) : t in T3];
[Inverse(phi5)(mu5(t)) : t in T5];



// Case 2
_<x>:=PolynomialRing(Rationals());
C7:=HyperellipticCurve((x^2+1)*(x^2+2*x-1));
C8:=HyperellipticCurve(-(x^2+1)*(x^2+2*x-1));
E7, phi := EllipticCurve(C7, Points(C7 : Bound:=100)[1]);
E8, phi := EllipticCurve(C8, Points(C8 : Bound:=100)[1]);

rank7, flag7 := Rank(E7);
rank8, flag8 := Rank(E8);

// Claim 8: E7 and E8  have rank 1
assert flag7 and flag8;
assert rank7 eq 1 and rank8 eq 1;
print "Claim 8 succesfully verified!";



// Case 2.a.i
P<v> := PolynomialRing(Rationals());
m := v^2+1;
n := 2*v;
u := v^2-1;
z := m^2+n^2;
r := m^2-n^2;
s := 2*m*n;

print "Claim 9: Expressing r^2+2*r*s-s^2 in terms of v";
print r^2+2*r*s-s^2;

C9 := HyperellipticCurve(r^2+2*r*s-s^2);
C10 := HyperellipticCurve(-(r^2+2*r*s-s^2));

auts9 := Automorphisms(C9);
auts10 := Automorphisms(C10);

D9, phi9 := CurveQuotient(AutomorphismGroup(C9, [auts9[4]]));
D10, phi10 := CurveQuotient(AutomorphismGroup(C10, [auts10[4]]));

print "Claim 10: equation for genus 2 quotient curve";
assert IsIsomorphic(D9, D10);
print D10;

J10 := Jacobian(D10);
// Claim 11: The Mordell-Weil group of J10 is Z/2Z x Z
MW10, mu10, flag1, flag2 := MordellWeilGroup(J10);
assert flag1 and flag2;
assert Invariants(MW10) eq [ 2, 0 ];
print "Claim 11 successfully verified!";

print "Claim 12: the points on D10";
ptsD10 :=  Chabauty(mu10(MW10.2));
print ptsD10;

ptsD9 := Points(D9:  Bound := 10);
assert #ptsD9 eq #ptsD10;
// since D9 and D10 are isomorphic this means we also have all points on D9

print "Claim 13: the points on C9";
for p in ptsD9 do
inverse_image := Support(Pullback(phi9, Place(p)));
print <p,[<q,Degree(q)> : q in inverse_image]>;
end for;

print "Claim 14: the points on C10";
for p in ptsD10 do
inverse_image := Support(Pullback(phi10, Place(p)));
print <p,[<q,Degree(q)> : q in inverse_image]>;
end for;


P<v> := PolynomialRing(Rationals());

// Case 2.a.ii
// equations in subcase A)
w := 1;
n := v^2;
m := 2*w^2;
z := m^2+n^2;
r := m^2-n^2;
s := 2*m*n;

print "Claim 15 A) expanding in terms of v:";
print r^2 + 2*r*s - s^2;

C11 := HyperellipticCurve(r^2 + 2*r*s - s^2);
C12 := HyperellipticCurve(-(r^2 + 2*r*s - s^2));

// equations in subcase A)
n := 2*v^2;
m := w^2;
z := m^2+n^2;
r := m^2-n^2;
s := 2*m*n;

print "Claim 15 B) expanding in terms of v:";
print r^2 + 2*r*s - s^2;

C13 := HyperellipticCurve(r^2 + 2*r*s - s^2);
C14 := HyperellipticCurve(-(r^2 + 2*r*s - s^2));

// Sanity check, not an actual claim in the paper as this can be seen by hand.
assert IsIsomorphic(C11, C13);
assert IsIsomorphic(C12, C14);


auts11 := Automorphisms(C11);
auts12 := Automorphisms(C12);

D11,phi11 := CurveQuotient(AutomorphismGroup(C11, [auts11[4]]));
D12,phi12 := CurveQuotient(AutomorphismGroup(C12, [auts12[4]]));

print "Claim 16: equation for genus 2 quotient curve";
assert IsIsomorphic(D11, D12);
print D12;


J12 := Jacobian(D12);
lowerbound, upperbound := RankBounds(J12);
// Claim 17: rank is 0
assert upperbound eq 0;
print "Claim 17 successfully verified!";
print "Claim 18: the points on D12";
ptsD12 := Chabauty0(J12);
print ptsD12;

ptsD11 := Points(D11:  Bound := 10);
assert #ptsD11 eq 2;
assert #ptsD12 eq 2;
// since D11 and D12 are isomorphic this means we also have all points on D11


print "Claim 19: the points on C11";
for p in ptsD11 do
inverse_image := Support(Pullback(phi11, Place(p)));
print <p,[<q,Degree(q)> : q in inverse_image]>;
end for;

// Somehow the points returned by Chabauty0(J12) are points on a copy of D12.
// This causes problems so we recompute the points in a different way.
ptsD12 := Points(D12:  Bound := 10);
assert #ptsD12 eq 2;


print "Claim 20: the points on C12";
for p in ptsD12 do
inverse_image := Support(Pullback(phi12, Place(p)));
print <p,[<q,Degree(q)> : q in inverse_image]>;
end for;


// Case 2.b.i
P<v> := PolynomialRing(Rationals());
w := 1;
k := 1;
u := k*(v^2 - 2*v*w - w^2);
n := k*(2*v*w - 2*w^2);
m := k*(v^2+w^2);

r := m^2 - 2*m*n - n^2;
assert r eq u^2; // sanity check
s := m^2 + 2*m*n - n^2;
z := m^2 + n^2;


print "Claim 21 expanding in terms of v:";
print (r^2 + 2*r*s - s^2)/2;
C15 := HyperellipticCurve((r^2 + 2*r*s - s^2)/2);
C16 := HyperellipticCurve(-(r^2 + 2*r*s - s^2)/2);

auts15 := Automorphisms(C15);
auts16 := Automorphisms(C16);

D15,phi15 := CurveQuotient(AutomorphismGroup(C15, [auts15[4]]));
D16,phi16 := CurveQuotient(AutomorphismGroup(C16, [auts16[4]]));


print "Claim 22: equation for genus 2 quotient curve";
assert IsIsomorphic(D15, D16);
print D16;

J16 := Jacobian(D16);
// Claim 23: The Mordell-Weil group of J16 is Z/2Z x Z
MW16, mu16, flag1, flag2 := MordellWeilGroup(J16);
assert flag1 and flag2;
assert Invariants(MW16) eq [ 2, 0 ];
print "Claim 23 successfully verified!";

print "Claim 24: the points on D16";
ptsD16 :=  Chabauty(mu16(MW16.2));
print ptsD16;

ptsD15 := Points(D15:  Bound := 10);
assert #ptsD15 eq #ptsD16;
// since D15 and D16 are isomorphic this means we also have all points on D15

print "Claim 25: the points on C15";
for p in ptsD15 do
inverse_image := Support(Pullback(phi15, Place(p)));
print <p,[<q,Degree(q)> : q in inverse_image]>;
end for;

print "Claim 26: the points on C16";
for p in ptsD16 do
inverse_image := Support(Pullback(phi16, Place(p)));
print <p,[<q,Degree(q)> : q in inverse_image]>;
end for;

// Case 2.b.ii
P<v> := PolynomialRing(Rationals());
w := 1;
k := 1;
u := k*(v^2 - 2*v*w - w^2);
n := k*(v^2 + 2*v*w - w^2);
m := k*(2*(v^2+w^2)+v^2 + 2*v*w - w^2);
assert (m-n)/2 eq v^2+w^2; // sanity check

r := (m^2 - 2*m*n - n^2)/2;
assert r eq u^2; // sanity check
s := (m^2 + 2*m*n - n^2)/2;
z := (m^2 + n^2)/2;

print "Claim 27 expanding in terms of v:";
print (r^2 + 2*r*s - s^2)/2;
C17 := HyperellipticCurve((r^2 + 2*r*s - s^2)/2);
C18 := HyperellipticCurve(-(r^2 + 2*r*s - s^2)/2);

auts15 := Automorphisms(C17);
auts16 := Automorphisms(C18);

D17,phi17 := CurveQuotient(AutomorphismGroup(C17, [auts15[4]]));
D18,phi18 := CurveQuotient(AutomorphismGroup(C18, [auts16[4]]));


print "Claim 28: equation for genus 2 quotient curve";
assert IsIsomorphic(D17, D18);
print D18;

J18 := Jacobian(D18);
// Claim 29: The Mordell-Weil group of J16 is Z/2Z x Z
MW18, mu18, flag1, flag2 := MordellWeilGroup(J18);
assert flag1 and flag2;
assert Invariants(MW18) eq [ 2, 0 ];
print "Claim 29 successfully verified!";

print "Claim 30: the points on D18";
ptsD18 :=  Chabauty(mu18(MW18.2));
print ptsD18;

ptsD17 := Points(D17:  Bound := 10);
assert #ptsD17 eq #ptsD18;
// since D17 and D18 are isomorphic this means we also have all points on D17

print "Claim 31: the points on C17";
for p in ptsD17 do
inverse_image := Support(Pullback(phi17, Place(p)));
print <p,[<q,Degree(q)> : q in inverse_image]>;
end for;

print "Claim 32: the points on C18";
for p in ptsD18 do
inverse_image := Support(Pullback(phi18, Place(p)));
print <p,[<q,Degree(q)> : q in inverse_image]>;
end for;

exit;
