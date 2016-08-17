function TwistOfSimpleModularSymbolsSpace(S,chi);
//On input of a irreducible new modular symbols space S and a character chi of prime power 
//modulus output the modular symbols space corresponding to the twist of M by the primitive character associated to chi
//i.e. if f is the modular form associated to S then outputs the modular symbols space
//corresponding to the modular form f_chi
//the sign of M should be 1 or -1, and the sign of the output will be the same
  if Conductor(chi) eq 1 then;
    return S;
  end if;
  chi := AssociatedPrimitiveCharacter(chi);
  chi_S := DirichletCharacter(S);
  m := Modulus(chi_S);
  n := Modulus(chi);
  assert IsPrimePower(n);
  p := PrimeDivisors(n)[1];
  sign := Sign(S);
  k := Weight(S);
  chi_t := Extend(chi_S*chi^2,m*n^2);
  Mt := ModularSymbols(chi_t,k,sign);
  St := CuspidalSubspace(Mt);
  for Si in NewformDecomposition(St) do;
    Snew := AssociatedNewSpace(Si);
    tf,chi_i := IsTwist(S,Snew,p);
    if tf and AssociatedPrimitiveCharacter(chi_i) eq chi then;
      return Snew;
    end if;
  end for;
  print "Did not find a twist while we should have!!!!";
  assert false;
end function;




function PostiveRankNewFactors(m,n);
  //Let m,n be two integer and let 
  //G be the congruence subgroup given by the matrices of the form
  //
  //   [a b]
  //   [c d]
  //
  //with a,d congruent to 1 modulo mn
  //and c congruent to 0 modulo m^2n
  //Then this function returns one modular symbols space
  //for every isogeny class of simple abelian varieties that occurs as an isogeny factor
  //of J(G) and obtains positive rank over Q(zeta_m)
  pr_new_factors := [];
  D := FullDirichletGroup(m*n);
  Chi := FullDirichletGroup(m);
  for chi in Elements(Chi) do;
    for d in GaloisConjugacyRepresentatives(D) do;
      d1 := Extend(d,m^2*n);
      //if IsOdd(d) then;
      //  continue;
      //end if;

      M := ModularSymbols(d1,2,1);
      S := CuspidalSubspace(M);
      for Si in NewformDecomposition(S) do;
        Snew := AssociatedNewSpace(Si);
        St := TwistOfSimpleModularSymbolsSpace(Snew,chi);
        if Dimension(St) ne Dimension(WindingSubmodule(St)) then;
          Append(~pr_new_factors,Snew);
        end if;
      end for;
    end for;
  end for;
  return pr_new_factors;
end function;





