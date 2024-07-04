load "elliptic_curve_chabauty.m";

_<x>:=PolynomialRing(Integers());
f := x*(x^2+1)*(x^2+2*x-1);
// a model of X_1(16)
C := HyperellipticCurve(f);

FF<x,y> := FunctionField(C);
p1,p2 := Explode(Zeros(x-1)); p1,p2;

CC := RegularModel(C,2);
IntersectionMatrix(CC);
ComponentGroup(CC);
P1_CC := PointOnRegularModel(CC, RepresentativePoint(p1));
P2_CC := PointOnRegularModel(CC, RepresentativePoint(p2));

//the following show that the two rational cusps p1,p2 lie on the same component
P1_CC`component_indices;
P2_CC`component_indices;



tf, h := IsPrincipal(5*(p1-p2));

// h is 0 mod 2 and 3 so we need to divide by 6

print h/6;
//prints
// (-6*x^2 - 4*x + 2)/(x^5 - 5*x^4 + 10*x^3 - 10*x^2 + 5*x - 1)*y + (x^5 + 13*x^4 -
//    2*x^3 + 10*x^2 - 7*x + 1)/(x^5 - 5*x^4 + 10*x^3 - 10*x^2 + 5*x - 1)
// exceptions give quadratic points on the degree 5 cover X -> X_1(16) given by z^5 = (h/6);
// note that t(x,y,z) = (x,-y,z^(-1)) is an automorphism of X, and a quadratic point on
// X with x in Q gives a rational point on X/<t>


//first compute the curve whose rational points we need to find

A3<x,y,z>:=AffineSpace(Rationals(),3);

C := Curve(A3,[y^2-(x*(x^2+1)*(x^2+2*x-1)), z^5*(x^5 - 5*x^4 + 10*x^3 - 10*x^2 + 5*x - 1)-((-6*x^2 - 4*x + 2)*y + x^5 + 13*x^4 - 2*x^3 + 10*x^2 - 7*x + 1)]);

C1 := ProjectiveClosure(C);
G := AutomorphismGroup(C1);
auts := [g : g in Automorphisms(C1) | Order(G ! g) eq 2];
assert #auts eq 1;

// t = auts[1]
// C2 = X/<t>, we need to find the points on this curve
C2 := CurveQuotient(AutomorphismGroup(C1, [G ! auts[1]]));

print C2;
//  prints: 
//      Hyperelliptic Curve defined by y^2 = x^6 - x^5 + 5*x^3 - x + 1 over Rational Field

//  now for the determination of the points

_<x> :=PolynomialRing(Integers());
f := x^6 - x^5 + 5*x^3 - x + 1;
// Uncomment the following line if you want to skip computing C2.
C2 := HyperellipticCurve(f);
L1 :=  NumberField(x^3 - x^2 + 2*x + 2);
g := [pe[1] : pe in Factorization(ChangeRing(f,L1)) |  Degree(pe[1]) eq 4][1];


time Hk , AtoHk, expvecs, factorbase := TwoCoverDescent(C2: Raw:=true);
  
// Get the map that allows us to construct covers explicitly as elements in an algebra A
HktoA := createHktoA(AtoHk, expvecs, factorbase);

Hk := Setseq(Hk);
points := {@ @};
hk_info := [* *];
for i in [1..#Hk] do
  hk := Hk[i];
  print "doing hk", i, hk;
  assert AtoHk(HktoA(hk)) eq hk;
  
  success, gamma_g, new_points, message := EllipticChabauty(C2, g, hk, HktoA);
  print "chabauty result", success, gamma_g, new_points, message;

  if success then
    points := points join new_points;
    continue;
  end if;

  Append(~hk_info,<gamma_g,message>);
  if not success then
    print "=============================== failed ======================";
    print hk, i, #Hk;
    print hk_info;
    break;
  end if;
end for;

print "====================== found points: ====================";
print points;

