intrinsic HeckeFilter(X::MDCrvMod, q::RngIntElt, S::SeqEnum) -> SeqEnum
{ Return all the divisors D in S in such that (T_q-q<q>-1)D=0 }
  filtered := [];
  for D in S do
    A_qD := HeckeOperator(X,q,D)-q*DiamondOperator(X,q,D)-D;
    if IsPrincipal(A_qD) then
     Append(~filtered, D);
    end if;
  end for;
  return filtered;
end intrinsic;

intrinsic HeckeFilter(X::MDCrvMod, Y::MDCrvMod, q::RngIntElt, S::SeqEnum) -> SeqEnum
{ Return all the divisors D in S in such that (T_q-q<q>-1)D'=0 where D' is the image
 of D in Y}
  filtered := [];
  for D in S do
    D1 := DegeneracyMap(X, Y, D);
    time A_qD := HeckeOperator(Y,q,D1)-q*DiamondOperator(Y,q,D1)-D1;
    if IsPrincipal(A_qD) then
     print Index(S,D);
     Append(~filtered, D);
    end if;
  end for;
  return filtered;
end intrinsic;