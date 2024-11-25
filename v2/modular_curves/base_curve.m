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
    Ex := EllipticCurve([Evaluate(f, x) : f in X`_E]);
    return Ex;
end intrinsic;




