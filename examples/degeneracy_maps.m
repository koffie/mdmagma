AttachSpec("mdmagma.spec");

X := MDX1(22,GF(3));
Y := MDX1(11,GF(3));

// we compute the two degeneracy maps X -> Y explicitly and verify
// that the resulting quadratic form should be 3y^2+3x^2-4xy

E := EllipticCurve(X`_E);
P := (E ! [0,0]);
xy := TateNormalForm_xy(E,2*P);
f := map< Curve(X) -> Curve(Y) | xy cat [1] >;
Degree(f);
assert Degree(f1) eq 3;


KX := BaseRing(E);
R<x> := PolynomialRing(KX);
E1,phi := IsogenyFromKernel(E,x-(11*P)[1]);
P1 := phi(P);
xy1 := TateNormalForm_xy(E1,2*P1);
f1 := map< Curve(X) -> Curve(Y) | xy1 cat [1] >;
assert Degree(f1) eq 3;

YE, psi := EllipticCurve(Curve(Y));

g := f*psi;
g1 := f1*psi;

FFX := FunctionField(Curve(X));

Pg := [Evaluate(DefiningEquations(g)[i],[FFX.1,FFX.2,1]) : i in [1..3]];
Pg1 := [Evaluate(DefiningEquations(g1)[i],[FFX.1,FFX.2,1]) : i in [1..3]];

Q := (YE(FFX) ! Pg) + (YE(FFX) ! Pg1);

h := map< Curve(X) -> YE | Eltseq(Q) >;

assert Degree(h) eq 2;