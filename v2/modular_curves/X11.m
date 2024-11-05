//Code for the modular curve X_1(M,N)
declare type MDCrvMod11: MDCrvMod;
declare attributes MDCrvMod11: level, curve, base_ring, N, M;


intrinsic Print(X::MDCrvMod11)
{Print X}
   printf "The modular curve X_1(%o, %o) over %o", X`N, X`M, X`base_ring;
end intrinsic;

intrinsic MDX11(N::RngIntElt, M::RngIntElt, base_ring::Rng) -> MDCrvMod11
{ Create the modular curve X_1(N, M) }
  assert M mod N eq 0;
  X := New(MDCrvMod11);
  _InitMDCrvMod11(X, N, M, base_ring);
  return X;
end intrinsic;

intrinsic _InitMDCrvMod11(X::MDCrvMod11, N::RngIntElt, M::RngIntElt, base_ring::Rng)
{Initialize an MDCrvMod1 object}
    X`N := N;
    X`M := M;
    _InitMDCrvMod(X, M,  base_ring);
end intrinsic;