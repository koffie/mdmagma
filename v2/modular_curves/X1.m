
declare type MDCrvMod1: MDCrvMod;
declare attributes MDCrvMod1: level, curve, base_ring, N;
//declare attributes MDCrvMod1: N;


intrinsic Print(X::MDCrvMod1)
{Print X}
   printf "The modular curve X_1(%o) over %o", X`N, "<basering not implemented>"; // todo BaseRing(X`curve);
end intrinsic;

intrinsic MDX1(N::RngIntElt, base_ring::Rng) -> MDCrvMod1
{ Create the modular curve X_1(N) }
  X := New(MDCrvMod1);
  _InitMDCrvMod1(X, N, base_ring);
  return X;
end intrinsic;

intrinsic _InitMDCrvMod1(X::MDCrvMod1, N::RngIntElt, base_ring::Rng)
{Initialize an MDCrvMod1 object}
    X`N := N;
    // todo curve := _equation_X(.....);
    // _InitMDCrvMod(X, N,  curve);
end intrinsic;