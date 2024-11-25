X := MDX11(2, 20, GF(3));

procedure test_Print()
  TSTAssertEQ(
    Sprint(X),
    "The modular curve X_1(2, 20) over Finite field of size 3"
  );
end procedure;

procedure test_Level()
  TSTAssertEQ(Level(X), 20);
end procedure;

procedure test_Genus()
  TSTAssertEQ(Genus(X), 9);
end procedure;

procedure test_LevelStructure()
   x := Places(Curve(X), 4)[2];
   PQ := LevelStructure(X, x);
   TSTAssertEQ(Order(PQ`P), 2);
   TSTAssertEQ(Order(PQ`Q), 20);
   TSTAssertNE(PQ`P,10*PQ`Q);
end procedure;

procedure test_ModuliPoint()
   x := Places(Curve(X), 4)[2];
   E := EllipticCurve(X,x);
   PQ := LevelStructure(X, x);
   TSTAssertEQ(ModuliPoint(X, E, PQ),x);
end procedure;

test_Print();
test_Level();
test_Genus();
test_LevelStructure();
test_ModuliPoint();