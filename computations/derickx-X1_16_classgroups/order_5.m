load "elliptic_curve_chabauty.m";

_<x>:=PolynomialRing(Integers());
f16 := x*(x^2+1)*(x^2+2*x-1);
// a model of X_1(16)
X1_16 := HyperellipticCurve(f16);

FF<x,y> := FunctionField(X1_16);
p1 := Place(X1_16 ! [1,2]);
p2 := Place(X1_16 ! [1,-2]);


// Claim 1: p1 - w(p1) is or order 5, note p2 = w(p1);
is_of_order_5, g := IsPrincipal(5*(p1-p2));
assert is_of_order_5;
print "Claim 1 successfully verified";

// g is 0 mod 2 and 3 so we need to divide by 6
g := g/6;

// g evaluated at the point [0,0] is a unit so there are no vertical
// fibers in the support of div(g)
assert (g)(X1_16 ! [0,0]) eq -1;


print "Claim 2: equation for g";
print g;
// prints
// (-6*x^2 - 4*x + 2)/(x^5 - 5*x^4 + 10*x^3 - 10*x^2 + 5*x - 1)*y + (x^5 + 13*x^4 -
//    2*x^3 + 10*x^2 - 7*x + 1)/(x^5 - 5*x^4 + 10*x^3 - 10*x^2 + 5*x - 1)

XX1_16 := RegularModel(X1_16,2);
// Claim 3: there are irreducible components over F_2, both occurring with multiplicity 1
assert Multiplicities(XX1_16) eq [1,1];
print "Claim 3 successfully verified";

P1_XX := PointOnRegularModel(XX1_16, RepresentativePoint(p1));
P2_XX := PointOnRegularModel(XX1_16, RepresentativePoint(p2));

print P1_XX`component_indices;

// Claim 4: two rational points p1,p2 lie on the same component
assert P1_XX`component_indices eq P2_XX`component_indices;
print "Claim 4 successfully verified";




// exceptions give quadratic points on the degree 5 cover Y_{5,g,1} -> X_1(16) given by z^5 = g;
// note that tau(x,y,z) = (x,-y,z^(-1)) is an automorphism of Y_{5,g/6,1}, and a quadratic point on
// Y_{5,g,1} with x in Q gives a rational point on X/<tau>


//first compute the curve whose rational points we need to find

A3<x,y,z>:=AffineSpace(Rationals(),3);

// The curve Y_{5,g,1}
Y_5g1 := Curve(A3,[y^2-(x*(x^2+1)*(x^2+2*x-1)), z^5*(x^5 - 5*x^4 + 10*x^3 - 10*x^2 + 5*x - 1)-((-6*x^2 - 4*x + 2)*y + x^5 + 13*x^4 - 2*x^3 + 10*x^2 - 7*x + 1)]);


FF := FunctionField(Y_5g1);

// z+z^-1 is invariant under tau
// so f below defines a curve to which Y_{5,g,1}/tau has a map
f :=  MinimalPolynomial(FF.3+FF.3^-1);
_<y> := Parent(f);
f;
// prints y^5 - 5*y^3 + 5*y + (-2*x^5 - 26*x^4 + 4*x^3 - 20*x^2 + 14*x - 2)/(x^5 - 5*x^4 + 10*x^3 - 10*x^2 + 5*x - 1)


PY_5g1 := ProjectiveClosure(Y_5g1);
G := AutomorphismGroup(PY_5g1);
auts := [g : g in Automorphisms(PY_5g1) | Order(G ! g) eq 2];
assert #auts eq 1;


// Y16 = Y_{5,g,1}/tau, we need to find the points on this curve
Y16, h := CurveQuotient(AutomorphismGroup(PY_5g1, [G ! auts[1]]));

print Y16;
// prints: 
//     Hyperelliptic Curve defined by y^2 = x^6 - x^5 + 5*x^3 - x + 1 over Rational Field

// We look at the curve defined by f above
A2<x,y>:=AffineSpace(Rationals(),2);
f :=  (y^5 - 5*y^3 + 5*y)*(x^5 - 5*x^4 + 10*x^3 - 10*x^2 + 5*x - 1) + (-2*x^5 - 26*x^4 + 4*x^3 - 20*x^2 + 14*x - 2);
Y16_1 := Curve(A2,f);
Y16_2 := ProjectiveClosure(Y16_1);
_,Y16_3, phi := IsHyperelliptic(Y16_2);

// Y16 has a map to Y16_1 since FF.3+FF.3^-1 is invariant under tau. The following shows that the
// normalisation of Y16_1 is binational to Y16.
assert Genus(Y16) eq Genus(Y16_1);
// The curves Y16_2 and Y16_3 are binational to Y16_1 by construction. In particular, all the four
// curves Y16, Y16_1, Y_16_2 and Y_16_3 are birational models for the same curve.


print "Claim 5: The equation for Y_{16} in the paper";
print Y16_3;
// prints: 
//     Hyperelliptic Curve defined by y^2 = x^6 + x^5 - 5*x^3 + x + 1 over Rational Field


// its easy to see that (x,y) -> (-x,y) is an isomorphism between Y16 and Y16_3.

FFY16_2 := FunctionField(Y16_2);


// Inverse(phi) doesn't work
// IsIsomorphic(Y16_2,Y16_3) and Pusforward(phi,Divisor(FFY16_2.1)) seems to take forever
// The block of code below is just to get the x coordinate on Y16_2 as a function
// on Y16_3; Suggestions for a better/faster solution are welcome.

K := QuadraticField(5);
Y16_3K := ChangeRing(Y16_3,K);
pts3 := Points(Y16_3K : Bound := 100);
plcs3 := [];
for P in pts3 do
  P1 := Place(Y16_3(K) ! P);
  if &or[P1 eq P2 : P2 in plcs3] then continue; end if;
  Append(~plcs3, P1);
end for;
plcs2 := [Support(Pullback(phi,P))[1] : P in plcs3];
D := &+[Valuation(FFY16_2.1,plcs2[i])*plcs3[i] : i in [1..#plcs2]];
_,h16 := IsPrincipal(D);
print x/Pullback(phi,h16);  
// the above prints 3/2, if for some reason this is no longer the case replace 3/2 in the
// two lines below by the fraction that was printed instead.       
assert x/Pullback(phi,h16) eq 3/2;
h16 := 3/2*h16;
// h16 is now the x-coordinate on Y16_2 expressed as a function on Y16_3


R<z,y> := Parent(h16);
print "Claim 6: equation for h16+1";
// note the +1 is merely to get a smaller looking equation for the paper.
print h16+1;
// prints:           
// (2*z^2 - 6*z + 2)/(z^5 - 5*z^4 + 5*z^3 + 5*z^2 - 5*z - 3)*y + (10*z^2 - 10*z - 
//    2)/(z^5 - 5*z^4 + 5*z^3 + 5*z^2 - 5*z - 3)


J16_3 := Jacobian(Y16_3);
// Claim 7: The Jacobian has rank 3. This means that both abelian and
// quadratic chabauty are not possible
lower_bound, upper_bound := RankBounds(J16_3); 
assert lower_bound eq 3;
assert upper_bound eq 3;
print "Claim 7 successfully verified";

print "Claim 8: The automorphism group of Y_16:";
GroupName(AutomorphismGroup(Y16_3));

E := CurveQuotient(AutomorphismGroup(Y16_3,[Automorphisms(Y16_3)[3]]));
print "Claim 9: Equation for the elliptic curve E";
print E;
// Claim 10: rank of the elliptic curve
lower_bound, upper_bound := RankBounds(E);
assert lower_bound eq 2;
assert upper_bound eq 2;
print "Claim 10 successfully verified";

// now for the determination of the points on Y16_3

_<x> :=PolynomialRing(Integers());
f := x^6 + x^5 - 5*x^3 + x + 1;
K<alpha> :=  NumberField(x^3 - x^2 + 2*x + 2);
g1,g2 := Explode([pe[1] : pe in Factorization(ChangeRing(f,K))]);

R<z> := Parent(g1);
print "Claim 11: factorisation of z^6 + z^5 - 5z^3 + z + 1";
print g1,g2;


print "Starting two cover descent + elliptic curve chabauty";

Hk , AtoHk, expvecs, factorbase := TwoCoverDescent(Y16_3: Raw:=true);
  
// Get the map that allows us to construct covers explicitly as elements in an algebra A
HktoA := createHktoA(AtoHk, expvecs, factorbase);

Hk := Setseq(Hk);
points := {@ @};
hk_info := [* *];
for i in [1..#Hk] do
  hk := Hk[i];
  print "doing hk", i, hk;
  assert AtoHk(HktoA(hk)) eq hk;
  
  success, gamma_g, new_points, message := EllipticChabauty(Y16_3, g2, hk, HktoA);
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


 
points := [Y16_3 ! P : P in points];

print "Claim 12: the rational points on Y_{16}";
print points;

print "Claim 13: table of points where divisibility by 5 could potentially fail";
for P in points do
  xP := h16(P);
  if xP eq Infinity() then
     yP := xP;
     d := 1;
   else
     yP := Evaluate(f16,xP);
     if yP eq 0 then
       d := 1;
     else
       d := SquarefreeFactorization(Numerator(yP)*Denominator(yP));
     end if;
  end if;
  print P, xP, yP, d, ClassNumber(QuadraticField(d));
end for;




// Finally let's study the ideal classes in Q(\sqrt{-2030}) coming from the element of order 5 
// in J_1(16)
_<t> :=PolynomialRing(Rationals());     
x := 29/242;
K<y> := NumberField(t^2-x*(x^2+1)*(x^2+2*x-1));
g := (-6*x^2 - 4*x + 2)/(x^5 - 5*x^4 + 10*x^3 - 10*x^2 + 5*x - 1)*y + (x^5 + 13*x^4 -
    2*x^3 + 10*x^2 - 7*x + 1)/(x^5 - 5*x^4 + 10*x^3 - 10*x^2 + 5*x - 1);

// the ideal class associated with 29/242
I1 := &*[pe[1]^Sign(pe[2]) : pe in Factorization(MaximalOrder(K)*g)];

x := -242/29;
K<y> := NumberField(t^2-x*(x^2+1)*(x^2+2*x-1));
g := (-6*x^2 - 4*x + 2)/(x^5 - 5*x^4 + 10*x^3 - 10*x^2 + 5*x - 1)*y + (x^5 + 13*x^4 -
    2*x^3 + 10*x^2 - 7*x + 1)/(x^5 - 5*x^4 + 10*x^3 - 10*x^2 + 5*x - 1);

// the ideal class associated with -242/29
I2 := &*[pe[1]^Sign(pe[2]) : pe in Factorization(MaximalOrder(K)*g)];

// Claim 14: the constructed ideal classes are trivial in Q(\sqrt{-2030})
assert IsPrincipal(I1);
assert IsPrincipal(I2);
print "Claim 14 successfully verified";


