function X_1_m_n(m,n,base_ring : equation_directory:="models_X1_m_n")
    //Input: m,n - integers such that m divides n
    //       base_ring - a ring
    //Output: C - a curve
    //Returns an algebraic model C of the modular curve X_1(m,n) as a curve over base_ring
    assert IsDivisibleBy(n,m);
    n_str := IntegerToString(n);
    m_str := IntegerToString(m);
    file_name := equation_directory cat "/X1_" cat m_str cat "_" cat n_str cat ".txt";
    data := Read(file_name);
    data := Split(data);
    if m eq 2 then;
        //example contents of the file X1_2_10.txt
        //N := 5;
        //X := v^2 + (u^2 - 1)*v - 1;
        //q := 1/u;
        //t := -4*u/(u^2*v + u^2 - v + 3);
        //E:=[0,t^2-2*q*t-2,0,-(t^2-1)*(q*t+1)^2,0];
        //P:=[(t+1)*(q*t+1),t*(q*t+1)*(t+1)];
        //Q:=[0,0];

        A<u,v> := AffineSpace(base_ring,2);
        N := eval Substring(data[1], 6 ,#data[1]-6);
        X := eval Substring(data[2], 6 ,#data[2]-6);
        q := eval Substring(data[3], 6 ,#data[3]-6);
        t := eval Substring(data[4], 6 ,#data[4]-6);
        E:=eval Substring(data[5], 4 ,#data[5]-4);
        P:=eval Substring(data[6], 4 ,#data[6]-4);
        Q:=eval Substring(data[7], 4 ,#data[7]-4);
        C := Curve(A,X);
    else;
        print "m > 2 not yet implemented";
    end if;

    return ProjectiveClosure(C),E,P,Q;
end function;
     

