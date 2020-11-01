function X_1_n(n,base_ring : equation_directory:="models_X1_n")
    //Input: n - integer
    //       base_ring - a ring
    //       equation_directory - directory with files FFFc<n>.txt containing models
    //Output: C - a curve
    //Returns an algebraic model C of the modular curve X_1(m,n) as a curve over base_ring
    
    n_str := IntegerToString(n);
    file_name := equation_directory cat "/FFFc" cat n_str cat ".txt";
    data := Read(file_name);

    A<x,y> := AffineSpace(base_ring,2);
    X := eval(data);
    C := Curve(A,X);
    return ProjectiveClosure(C);
end function;