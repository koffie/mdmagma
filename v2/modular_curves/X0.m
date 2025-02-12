declare type MDCrvMod0: MDCrvMod;
declare attributes MDCrvMod0:level, curve, base_ring, N;


intrinsic Print(X::MDCrvMod0)
{Print X}
   base_ring := BaseRing(X);
   printf "The modular curve X_0(%o) over %o", X`N,  base_ring;
end intrinsic;

intrinsic MDX0(N::RngIntElt, base_ring::Rng) -> MDCrvMod0
{ Create the modular curve X_0(N) }
  assert N ge 10; //X_0(N) is only implemented when it is of gonality > 3
  X := New(MDCrvMod0);
  _InitMDCrvMod0(X, N, base_ring);
  return X;
end intrinsic;

intrinsic _InitMDCrvMod0(X::MDCrvMod0, N::RngIntElt,
    base_ring::Rng)
{Initialize an MDCrvMod0 object}
    X`N := N;
    X`base_ring := base_ring;
    _InitMDCrvMod(X, N);
end intrinsic;

intrinsic BaseRing(X::MDCrvMod0) -> Rng
{ Return the base ring over which this modular curve is defined }
    return (assigned X`base_ring) select X`base_ring  else "<unassigned>";
end intrinsic;
