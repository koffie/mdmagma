function X_1_m_n(m,n,base_ring : equation_directory:="models_X1_m_n", zeta_m:=0)
    //Input: m,n - integers such that m divides n
    //       base_ring - a ring
    //       equation_directory - directory with files X1_m_n.txt containing models
    //       zeta_m - a primitive mth root of unity in the base_ring (if unspecified one will be chosen)
    //Output: C - a curve
    //Returns an algebraic model C of the modular curve X_1(m,n) as a curve over base_ring
    assert IsDivisibleBy(n,m);
    if m gt 2 then
        if zeta_m ne 0 then
            assert zeta_m^m eq 1 and &and[zeta_m^e ne 1: e in Divisors(m)| e ne m];
        else
            zeta_m := RootOfUnity(m,base_ring);
        end if;
        z:=zeta_m; i:=zeta_m;
    end if;
    n_str := IntegerToString(n);
    m_str := IntegerToString(m);
    file_name := equation_directory cat "/X1_" cat m_str cat "_" cat n_str cat ".txt";
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
        if line[1] eq "t" then t := eval(val); end if;
        if line[1] eq "E" then E := eval(val); end if;
        if line[1] eq "P" then P := eval(val); end if;
        if line[1] eq "Q" then Q := eval(val); end if;
    end for;
    C := Curve(A,X);
    return ProjectiveClosure(C),E,P,Q;
end function;
