intrinsic TwoGenerators(x::FldFinElt,y::FldFinElt) -> RngMPolElt, RngMPolElt
{ Assumes x,y are elements of GF(p,n) and return generators of the ideal
  in PolynomialRing(GF(p),2) of bivariate polynomial vanishing at x,y
}
  Fmax := Parent(x);
  assert Fmax eq Parent(y);
  F := BaseField(Fmax);
  f := MinimalPolynomial(x);
  FF<t> := sub<Fmax|x>;
  g := MinimalPolynomial(y,FF);
  // at this point f, and g are generators, however we need to make them elements
  // of the polynomial ring R[u,v];
  R<u,v> := PolynomialRing(F,2);
  f1:=Evaluate(f,u);
  g1:=R!0;
  for i in [0..Degree(g)] do
    ai:=Eltseq(Coefficient(g,i));
    ai := &+ [ai[j]*u^(j-1) : j in [1..#ai]];
	g1:=g1+ai*v^i;
  end for;
  return f1,g1;
end intrinsic;

intrinsic Place(coordinates::[FldFunFracSchElt], values::[FldElt]) -> PlcCrvElt
{ Return the place on the curve corresponding coordinates[i]-value[i]=0, the values
  are allowed to lie in a field extension
}
  assert #coordinates eq 2;
  assert #values eq 2;
  a,b := Explode(values);
  f,g := TwoGenerators(a,b);
  // the following doesn't work:
  //   return Place([Evaluate(fi, coordinates) : fi in [f,g]]);
  // because of the error:
  //   Runtime error in 'Place': Elements of argument must be contained in the maximal
  //   finite or infinite equation order.
  f1 := Evaluate(f, coordinates);
  g1 := Evaluate(g, coordinates);
  places := [x : x in Zeros(f1) | Evaluate(g1,x) eq 0];
  assert #places eq 1;
  return places[1];

end intrinsic