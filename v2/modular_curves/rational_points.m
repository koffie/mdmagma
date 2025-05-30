intrinsic HeckeSieve(X::MDCrvMod, q::RngIntElt, S::SeqEnum) -> SeqEnum
{ Return all the divisors D in S in such that (T_q-q<q>-1)D=0 }
  sieved := [];
  for D in S do
    A_qD := HeckeOperator(X,q,D)-q*DiamondOperator(X,q,D)-D;
    if IsPrincipal(A_qD) then
     Append(~sieved, D);
    end if;
  end for;
  return sieved;
end intrinsic;

intrinsic HeckeSieve(X::MDCrvMod, Y::MDCrvMod, q::RngIntElt, S::SeqEnum : i:=0) -> SeqEnum
{ Return all the divisors D in S in such that (T_q-q<q>-1)D'=0 where D' is the image
 of D in Y, i is an optional parameter passed to DegeneracyMap }
  sieved := [];
  for D in S do
    if i eq 0 then
      D1 := DegeneracyMap(X, Y, D);
    else
      D1 := DegeneracyMap(X, Y, D : i:=i);
    end if;
    A_qD := HeckeOperator(Y,q,D1)-q*DiamondOperator(Y,q,D1)-D1;
    if IsPrincipal(A_qD) then
     Append(~sieved, D);
    end if;
  end for;
  return sieved;
end intrinsic;