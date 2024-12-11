X := MDX1(21,GF(2));

procedure test_print()
  TSTAssertEQ(
    Sprint(X),
    "The modular curve X_1(21) over Finite field of size 2"
  );
end procedure;

procedure test_Level()
  TSTAssertEQ(Level(X), 21);
end procedure;

procedure test_Genus()
  TSTAssertEQ(Genus(X), 5);
end procedure;

procedure test_LevelStructure()
   x := Places(Curve(X), 4)[1];
   L := LevelStructure(X, x);
   TSTAssertEQ(Order(L`P), 21);
end procedure;

procedure test_ModuliPoint()
   x := Places(Curve(X), 4)[1];
   E := EllipticCurve(X,x);
   P := LevelStructure(X, x);
   TSTAssertEQ(ModuliPoint(X, E, P), x);
end procedure;

procedure test_ApplyIsogeny(X,p)
   x := NoncuspidalPlaces(X, 4)[1];
   E := EllipticCurve(X, x);
   P := LevelStructure(X, x);
   phi := MDIsogenies(E,p)[4];
   phi_t:= DualIsogeny(phi);
   P1 := ApplyIsogeny(X, phi, P);
   P2 := ApplyIsogeny(X, phi_t, P1);
   E1 := Codomain(phi_t);
   // applying phi and then phi_t should
   // be multiplication by deg phi = p.
   TSTAssertEQ(E1 ! (5*P`P), P2`P);
end procedure;

procedure test_DiamondOperatorPlcCrvElt()
   x := Places(Curve(X), 5)[1];
   dx := DiamondOperator(X, 2, x);
   ddx := DiamondOperator(X, 10, dx);
   // 2*10 = -1 mod 21 so this should be the original point again
   TSTAssertEQ(x, ddx);
end procedure;

procedure test_CongruenceSubgroup()
   TSTAssertEQ(CongruenceSubgroup(X), Gamma1(21));
end procedure;

procedure test_Cusps()
   cusp_degrees := Sort([Degree(c) : c in Cusps(X)]);
   TSTAssertEQ(cusp_degrees, [ 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 6 ]);
end procedure;

procedure test_CuspSignature()
   cusp_signatures := Sort([CuspSignature(X, c) : c in Cusps(X)]);
   expected_signatures := [
       [ -17, 1 ],
       [ -13, 2 ],
       [ -5, 4 ],
       [ -3, 1 ],
       [ -1, 0 ],
       [ -1, 5 ],
       [ 1, 1 ],
       [ 1, 1 ],
       [ 1, 2 ],
       [ 5, -3 ],
       [ 11, -1 ],
       [ 19, -17 ]
   ];
   TSTAssertEQ(cusp_signatures, expected_signatures);
end procedure;

procedure test_CuspsWithSignature()
   cusps := CuspsWithSignature(X, [1,1]);
   TSTAssertEQ(#cusps, 2);
end procedure;

procedure test_IsCuspSignatureUnique()
   assert IsCuspSignatureUnique(MDX1(13, Rationals()));
   assert not IsCuspSignatureUnique(MDX1(13, GF(5)));
end procedure;


test_print();
test_Level();
test_Genus();
test_LevelStructure();
test_ModuliPoint();
// this seems to fail because of a magma bug
// X := MDX1(21,GF(2));
// p := 5;
// So we test using the following values instead
X := MDX1(17,GF(3));
p := 5;
test_ApplyIsogeny(X,p);
test_DiamondOperatorPlcCrvElt();
/*
the following are already tested for X11 and the code is generic
test_HeckeOperatorPlcCrvElt();
test_HeckeOperatorDivCrvElt();
test_DiamondOperatorDivCrvElt();
*/
test_CongruenceSubgroup();
test_Cusps();
test_CuspSignature();
test_CuspsWithSignature();
test_IsCuspSignatureUnique();
