intrinsic MDClassGroup(C::Crv, S::SeqEnum[DivCrvElt] : classgroup:=<>) -> GrpAb
{ The subgroup of the class group generated by the divisors in S }
    if #classgroup eq 0 then
        G,m1,m2 := ClassGroup(C);
    else
        G,m1,m2 := Explode(classgroup);
    end if;
    M := FreeAbelianGroup( #S );
    f := hom< M -> G | [m2(i) : i in S] >;
    // N := Kernel(f);
    // N := Lattice( Matrix( [ Eltseq(M ! i) : i in Generators(N)] ) );
    H := Image(f);
    return H,m1,m2;
end intrinsic;