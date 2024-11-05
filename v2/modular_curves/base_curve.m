// Base type for a all MDMagma modular curves
declare type MDCrvMod;
declare attributes MDCrvMod: level, curve, base_ring;

//intrinsic Level(X: MDCrvMod)

intrinsic Print(X::MDCrvMod)
{Print X}
    printf "The MDCrvMod class shouldn't be used directly";
end intrinsic;

intrinsic _InitMDCrvMod(X::MDCrvMod, level::RngIntElt, base_ring::Rng)
{Initialize an MDCrvMod object}
   X`level := level;
   X`base_ring := base_ring;
end intrinsic;




