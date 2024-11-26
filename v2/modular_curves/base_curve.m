// Base type for a all MDMagma modular curves
declare type MDCrvMod;
declare attributes MDCrvMod: level, curve, _E;

intrinsic Print(X::MDCrvMod)
{ Print X }
    printf "The MDCrvMod class shouldn't be used directly";
end intrinsic;

intrinsic _InitMDCrvMod(X::MDCrvMod, level::RngIntElt, curve::CrvPln, _E)
{Initialize an MDCrvMod object}
   assert GCD(level,Characteristic(BaseRing(curve))) eq 1;
   X`level := level;
   X`curve := curve;
   X`_E := _E;
end intrinsic;

intrinsic Level(X::MDCrvMod) -> RngIntElt
{ Return the level }
    return X`level;
end intrinsic;

intrinsic Curve(X::MDCrvMod) -> CrvPln
{ Return the underlying CrvPln object used for calculations with this modular curve }
    return X`curve;
end intrinsic;

intrinsic Genus(X::MDCrvMod) -> CrvPln
{ Return the genus of this modular curve }
    return Genus(X`curve);
end intrinsic;

intrinsic BaseRing(X::MDCrvMod) -> Rng
{ Return the base ring over which this modular curve is defined }
    return (assigned X`curve ) select BaseRing(X`curve)  else "<unassigned>";
end intrinsic;

intrinsic EllipticCurve(X::MDCrvMod, x::PlcCrvElt) -> CrvEll
{ Return the elliptic curve corresponding to the place x on X }
    E := EllipticCurve([Evaluate(f, x) : f in X`_E]);
    return E;
end intrinsic;

intrinsic HeckeOperator(X::MDCrvMod, p::RngIntElt, x::PlcCrvElt) -> DivCrvElt
{ Return the result of applying the hecke operator T_p on x as a divisor on X }
    ZZ := Integers();
    d := Degree(x);
    E := EllipticCurve(X, x);
    L := LevelStructure(X, x);
    isogenies := MDIsogenies(E,p);
    tp := [<Codomain(phi),ApplyIsogeny(X, phi, L)> : phi in isogenies];
    tp := Multiset([ModuliPoint(X, EL[1], EL[2]): EL in tp]);
    Tpx := &+[ (ZZ ! (Multiplicity(tp,x)*d/Degree(x)))*x : x in MultisetToSet(tp)];
    return Tpx;
end intrinsic;

intrinsic HeckeOperator(X::MDCrvMod, p::RngIntElt, D::DivCrvElt) -> DivCrvElt
{ Return the result of applying the hecke operator T_p on D as a divisor on X }
    a,b := Support(D);
    TpD := [b[i]*HeckeOperator(X, p, a[i]) : i in [1..#a]];
    return &+TpD;
end intrinsic;

intrinsic DiamondOperator(X::MDCrvMod, d::RngIntElt, x::PlcCrvElt) -> PlcCrvElt
{ Return the result of applying the diamond operator <d> on x as a place on X }
    E := EllipticCurve(X, x);
    L := LevelStructure(X, x);
    return ModuliPoint(X, E, DiamondOperator(X, d, L));
end intrinsic;

intrinsic DiamondOperator(X::MDCrvMod, d::RngIntElt, D::DivCrvElt) -> DivCrvElt
{ Return the result of applying the diamond operator <d> on D as a divisor on X }
    a,b := Support(D);
    dD := [b[i]*DiamondOperator(X, d, a[i]) : i in [1..#a]];
    return &+dD;
end intrinsic;


