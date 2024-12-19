procedure test_Isogenies()
    q := 5;
    F := GF(11);
    E := EllipticCurve([F ! 1, F ! 3]);
    j := jInvariant(E);
    isogenies := MDIsogenies(E, q);
    TSTAssertEQ(#isogenies, q+1);
    for phi in isogenies do
        TSTAssertEQ(jInvariant(Domain(phi)), j);
        TSTAssertEQ(Degree(phi), q);
    end for;
    jset := {jInvariant(Codomain(isogenies[i])) : i in [1..#isogenies]};
    TSTAssertEQ(#jset,q+1); //this might fail some times if one changes E or p
    // this is expected, this means there is "some CM".
end procedure;

procedure test_IsogeniesBugFix()
    // the following fails in magma V2.28-14, V2.28-15
    // and possibly some older versions as well.
    F := GF(2^4);
    // print F, MinimalPolynomial(F.1);
    R<x> := PolynomialRing(F);
    fp := x^12 + F.1^2*x^9 + F.1^3*x^8 + F.1^4*x^6 + F.1^6*x^4 + F.1^6*x^3 + F.1^7*x^2 + F.1^8*x + F.1^12;
    Fqn := SplittingField(fp);
    Fq2n := RandomExtension(Fqn, 2);
    // uncommenting the line below will make the error go away
    xi := Roots(ChangeRing(fp,Fqn));
    xi := Roots(ChangeRing(fp,Fq2n));

    E := EllipticCurve([0, F.1, F.1, 0, 0]);
    p := 5;

    Eq2n := BaseChange(E,Fq2n);
    fE := DefiningEquation(Eq2n);
    yi := [Roots(UnivariatePolynomial(Evaluate(Evaluate(fE,3,1),1,x[1])))[1][1] : x in xi];
    Pi := [E(Fq2n) ! [xi[i][1],yi[i]] : i in [1..#xi]];
    //now Pi is a list with all torsion points

    R<X> := PolynomialRing(Fq2n);
    kernel_polynomials := {&*[X-(j*P)[1] : j in [1..Ceiling((p-1)/2)]] : P in Pi};
    isogenies := [PowerStructure(MapSch) | ];
    for f in kernel_polynomials do;
        G := SubgroupScheme(Eq2n, f);
        /*
        print G, Order(G);
        pts := Points(G);
        print pts;
        print [n*P in pts : n in [2..5], P in pts];
        */
        Ef,phi := IsogenyFromKernel(Eq2n,f);

        Append(~isogenies,phi);

        DualIsogeny(phi);
    end for;
end procedure;

test_Isogenies();
// test_IsogeniesBugFix(); // todo run the test if magma is fixed