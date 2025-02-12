procedure test_MDRestriction()
  QQ := Rationals();
  V := VectorSpace(QQ,3);
  W := sub< V | [V.1,V.3]>;
  A := Matrix(QQ,[[1,0,0],[0,2,0],[0,0,3]]);
  B := Matrix(QQ,[[1,0],[0,3]]);
  B1 := MDRestriction(A, W);
  TSTAssertEQ(B1, B);
end procedure;

test_MDRestriction();