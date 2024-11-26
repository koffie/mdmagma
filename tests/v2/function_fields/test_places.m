procedure test_TwoGenerators(x, y, n)
    // n should be the degree of the field extension F(x,y)/F
    f,g := MDTwoGenerators(x,y);
    // f and g should vanish at x and y
    TSTAssertEQ(Evaluate(f,[x,y]), 0);
    TSTAssertEQ(Evaluate(g,[x,y]), 0);
    R := Parent(f);
    I := ideal<R | [f,g]>;
    h := HilbertPolynomial(Homogenization(I));
    // the following test ensures that f and g are actually generators
    TSTAssertEQ(h,n);

end procedure;

procedure test_Place(coordinates,values,P)
  P1 := MDPlace(coordinates,values);
  TSTAssertEQ(P1,P);
  for i in [1,2] do
    TSTAssertEQ(
      MinimalPolynomial(Evaluate(coordinates[i],P1)),
      MinimalPolynomial(values[i])
    );
  end for;
end procedure;


p := 7;
n := 12;
F := GF(p,n);
test_TwoGenerators(F.1^(7^6+1),F.1^(7^8+7^4+1),n);
test_TwoGenerators(F ! 1,F ! 2, 1);
test_TwoGenerators(F.1^(7^8+7^4+1)^5,F.1^(7^8+7^4+1)^3, 4);

F0 := GF(p);
P2<x,y,z> := ProjectivePlane(F0);
C := Curve(P2,x+y+z);
FF := FunctionField(C);
values := [F.1,-F.1-1];
test_Place([FF.1,FF.2], values, Places(C(F) ! values)[1]);
test_Place([FF.1^p,FF.2], [values[1]^(p),values[2]], Places(C(F) ! values)[1]);
test_Place(
  [Evaluate(MinimalPolynomial(F.1),FF.1),Evaluate(MinimalPolynomial(-F.1-1),FF.2)],
  [F! 0,F ! 0],
  Places(C(F) ! values)[1]
);
