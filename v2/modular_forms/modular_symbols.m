QQ := Rationals();

intrinsic MDCongruenceSubgroup(M :: ModSym) -> GrpPSL2
{
  Return the congruence subgroup that was used to define the ambient space of this modular
  symbols space.
} 
  M := AmbientSpace(M);
  if IsTrivial(DirichletCharacter(M)) then
    return Gamma0(Level(M));
  elif IsMultiChar(M) then
    if assigned M`isgamma1 and M`isgamma1 then
      return Gamma1(Level(M));
    elif assigned M`isgamma and M`isgamma then
      return CongruenceSubgroup(Integers()!Sqrt(Level(M)));
    end if;
  end if;
  require false : "Only implemented for Gamma0(N), Gamma1(N) and Gamma(N)";
end intrinsic;




intrinsic MDDiagonalALBasis(decomposition :: SeqEnum[ModSym]) -> SeqEnum[SeqEnum[ModFrmElt]], SeqEnum[SeqEnum[SeqEnum[FldRatElt]]]
{
  Given a decomposition of the modular symbols of weight k for some group Gamma0(N) it gives a corresponding basis
  of the modular in S_k(Gamma0(N)) respecting this decomposition such that further more the Atkin-Lehner operators 
  act diagonally. The second return value are the Atkin-Lehner eigenvalues on this basis
}
  if #decomposition eq 0 then
    return [],[];
  end if;
  G := MDCongruenceSubgroup(decomposition[1]);
  N := Level(G);
  require G eq Gamma0(N) : "Only implemented for Gamma0(N)";
  g :=  Genus(G);
  prec := 2*g;
  k := Weight(decomposition[1]);

  M := ModularForms(G,k);
  S := CuspidalSubspace(M);
  V, VtoS, StoV := RSpace(S);

  al_divisors := [pe[1]^pe[2] : pe in Factorisation(N)];
  al := [AtkinLehnerOperator(S,d) : d in al_divisors];
  q_basis := [qIntegralBasis(Si, prec) : Si in decomposition];

  V_basis := [[StoV(S ! f) : f in B] : B in q_basis];
  V_dec := [ChangeRing(sub< V | B>, QQ) : B in V_basis];

  V_basis_diag := [];
  al_diag := [[] : w_p in al];

  for Vi in V_dec do
    al_res :=  [MDRestriction(w_p, Vi) : w_p in al];
    diag, t :=  Diagonalisation(al_res);
    basis := [ClearDenominator(v) : v in Rows(t*BasisMatrix(Vi))];
    Append(~V_basis_diag, basis);
    for i in [1..#diag] do
      Append(~al_diag[i],Diagonal(diag[i]));
    end for;
  end for;


  q_basis_diag := [[VtoS(f) : f in B] : B in V_basis_diag];

  return q_basis_diag,al_diag;

end intrinsic;
