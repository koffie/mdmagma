//Code for the modular curve X_1(M,N)
declare type MDCrvMod11: MDCrvMod;
declare attributes MDCrvMod11: level, curve, base_ring, N, M;


intrinsic Print(X::MDCrvMod11)
{Print X}
   M := (assigned X`M ) select X`M else "<unassigned>";
   N := (assigned X`N ) select X`N else "<unassigned>";
   base_ring := BaseRing(X);
   printf "The modular curve X_1(%o, %o) over %o", M, N, base_ring;
end intrinsic;

intrinsic MDX11(N::RngIntElt, M::RngIntElt, base_ring::Rng) -> MDCrvMod11
{ Create the modular curve X_1(N, M) }
  assert M mod N eq 0;
  X := New(MDCrvMod11);
  _InitMDCrvMod11(X, N, M, base_ring);
  return X;
end intrinsic;

intrinsic _InitMDCrvMod11(X::MDCrvMod11, M::RngIntElt, N::RngIntElt,
    base_ring::Rng: equation_directory:="../models_X1_m_n", zeta_M:=0)
{Initialize an MDCrvMod1 object}
    X`M := M;
    X`N := N;
    curve := _equation_X11(M, N,base_ring : equation_directory:=equation_directory, zeta_m:=zeta_M);
    _InitMDCrvMod(X, N,  curve);
end intrinsic;

intrinsic _equation_X11(m,n,base_ring : equation_directory:="../models_X1_m_n", zeta_m:=0) -> CrvPln
{  Input: m,n - integers such that m divides n
          base_ring - a ring
          equation_directory - directory with files X1_m_n.txt containing models
          zeta_m - a primitive mth root of unity in the base_ring (if unspecified one will be chosen)
    Output: C - a curve
    Returns an algebraic model C of the modular curve X_1(m,n) as a curve over base_ring
}
    assert IsDivisibleBy(n,m);
    if m gt 2 then
        if zeta_m ne 0 then
            assert zeta_m^m eq 1 and &and[zeta_m^e ne 1: e in Divisors(m)| e ne m];
        else
            try
                zeta_m := RootOfUnity(m,base_ring);
            catch e
                printf "Specified base ring %o does not contain a %oth root of unity", base_ring, m;
                assert false;
            end try;
        end if;
        z:=zeta_m; i:=zeta_m;
    end if;
    n_str := IntegerToString(n);
    m_str := IntegerToString(m);
    file_name := equation_directory cat "/X1_" cat m_str cat "_" cat n_str cat ".txt";
    print GetCurrentDirectory();
    data := Read(file_name);
    data := Split(data);
    //example contents of the file X1_2_10.txt
    //X := v^2 + (u^2 - 1)*v - 1;
    //q := 1/u;
    //t := -4*u/(u^2*v + u^2 - v + 3);
    //E:=[0,t^2-2*q*t-2,0,-(t^2-1)*(q*t+1)^2,0];
    //P:=[(t+1)*(q*t+1),t*(q*t+1)*(t+1)];
    //Q:=[0,0];
    A<u,v> := AffineSpace(base_ring,2);
    for line in data do
        val := Split(Split(line,"=")[2],";")[1];
        if line[1] eq "X" then X := eval(val); end if;
        if line[1] eq "q" then q := eval(val); end if;
        if line[1] eq "r" then r := eval(val); end if;
        if line[1] eq "s" then s := eval(val); end if;
        if line[1] eq "t" then t := eval(val); end if;
        if line[1] eq "E" then E := eval(val); end if;
        if line[1] eq "P" then P := eval(val); end if;
        if line[1] eq "Q" then Q := eval(val); end if;
    end for;
    C := Curve(A,X);
    return ProjectiveClosure(C),E,P,Q;
end intrinsic;