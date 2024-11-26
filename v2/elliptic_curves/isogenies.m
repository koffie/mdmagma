intrinsic MDIsogenies(E::CrvEll, p::RngIntElt) -> SeqEnum[MapSch]
{ Returns the p+1 elliptic curves that are p-isogenous to an initial E
  curve defined over a finite field
}
    assert IsPrime(p);
    //assert p ne 2;
    Fq := BaseRing(E);
    fp := DivisionPolynomial(E,p);
    Fqn := SplittingField(fp);
    Fq2n := RandomExtension(Fqn, 2);
    xi := Roots(ChangeRing(fp,Fq2n));
    Eq2n := BaseChange(E,Fq2n);
    fE := DefiningEquation(Eq2n);
    yi := [Roots(UnivariatePolynomial(Evaluate(Evaluate(fE,3,1),1,x[1])))[1][1] : x in xi];
    Pi := [E(Fq2n) ! [xi[i][1],yi[i]] : i in [1..#xi]];
    //now Pi is a list with all torsion points

    R<X> := PolynomialRing(Fq2n);
    kernel_polynomials := {&*[X-(j*P)[1] : j in [1..Ceiling((p-1)/2)]] : P in Pi};
    isogenies := [PowerStructure(MapSch) | ];
    for f in kernel_polynomials do;
        Ef,phi := IsogenyFromKernel(Eq2n,f);

        Append(~isogenies,phi);
    end for;
    return isogenies;
end intrinsic;