
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

function Functions_xyrsbcF2F3(curve)
    //Input: curve - the modular curve X_1(N) as returned by the function X_1
    //Output: x,y,r,s,b,c,F2,F3 - The modular units as in http://arxiv.org/pdf/1307.5719v1.pdf
    FF := FunctionField(curve);
    x := FF.1;
    y := FF.2;
    r := (x^2*y-x*y+y-1)/x/(x*y-1);
    s := (x*y-y+1)/x/y;
    b := r*s*(r-1);
    c := s*(r-1);
    F3 := b;
    F2 := b/(16*b^2+(1-20*c-8*c^2)*b + c*(c-1)^3);
    return x,y,r,s,b,c,F2,F3;
end function;

function X_1_n_jInvariant(curve)
    x,y,r,s,b,c,F2,F3 := Functions_xyrsbcF2F3(curve);
    jNum := 4096*b^6 - 6144*b^5*c^2 - 6144*b^5*c + 12288*b^5 + 3840*b^4*c^4 + 3072*b^4*c^3 - 4608*b^4*c^2 - 15360*b^4*c + 13056*b^4 
    - 1280*b^3*c^6 + 768*b^3*c^5 + 1536*b^3*c^4 - 2048*b^3*c^3 + 8448*b^3*c^2 - 13056*b^3*c + 5632*b^3 + 240*b^2*c^8 - 
    768*b^2*c^7 + 384*b^2*c^6 + 384*b^2*c^5 + 2400*b^2*c^4 - 7680*b^2*c^3 + 8448*b^2*c^2 - 4224*b^2*c + 816*b^2 - 
    24*b*c^10 + 168*b*c^9 - 432*b*c^8 + 288*b*c^7 + 1008*b*c^6 - 3024*b*c^5 + 4032*b*c^4 - 3168*b*c^3 + 1512*b*c^2 - 
    408*b*c + 48*b + c^12 - 12*c^11 + 66*c^10 - 220*c^9 + 495*c^8 - 792*c^7 + 924*c^6 - 792*c^5 + 495*c^4 - 220*c^3 + 
    66*c^2 - 12*c + 1;
    jDen := 16*b^5 - 8*b^4*c^2 - 20*b^4*c + b^4 + b^3*c^4 - 3*b^3*c^3 + 3*b^3*c^2 - b^3*c;
    return jNum/jDen;
end function;


function TateNormalForm_bc(E,P);
//Return the b,c of the tate normal form of (E,P) as in equation (2) of http://arxiv.org/pdf/1307.5719v1.pdf
    assert P[3] eq 1;
    x0:=P[1];
    y0:=P[2];

    a1,a2,a3,a4,a6:=Explode(aInvariants(E));
    aa1:=a1;
    aa3:=2*y0+a3+a1*x0;
    aa2:=3*x0+a2;
    aa4:=3*x0^2+2*x0*a2+a4-a1*y0;

    aaa1:=2*aa4/aa3+aa1;
    aaa3:=aa3;
    aaa2:=aa2-(aa4/aa3)^2-aa1*aa4/aa3;


    b:=-aaa2^3/aaa3^2;
    c:=-(aaa1*aaa2-aaa3)/aaa3;
    return [b,c];
end function;

function TateNormalForm_xy(E,P);
//return the x,y of the tate normal form of (E,P) as in section 2.1 http://arxiv.org/pdf/1307.5719v1.pdf
    b,c := Explode(TateNormalForm_bc(E,P));
    r := b/c;
    s := c^2/(b-c);
    t := (r*s-2*r+1);
    x := (s-r)/t;
    y := t/(s^2-s-r+1);
    return [x,y];
end function;

function EllipticCurveFromX1Place(P);
//Returns the associated elliptic curve corresponding to a place on X1N
//the elliptic curve is guaranteed to be in tate normal form, so that
//0,0 is the point of order N. The point 0,0 is returned as optional second element
    X1N := Curve(P);
    x,y,r,s,b,c,F2,F3:=Functions_xyrsbcF2F3(X1N);
    bP:=Evaluate(b,P);
    cP:=Evaluate(c,P);
    E:=EllipticCurve([1-cP,-bP,-bP,0,0]);
    return E, E ! [0,0];
end function;

function X1PlaceFromEllipticCurve(X1N, E, P)
//Returns a place on X_1(N) given an elliptic curve and a point of order N
    K := BaseRing(E);
    xy := TateNormalForm_xy(E,P);
    dP := Places(X1N(K) ! xy);
    assert #dP eq 1;
    return dP[1];
end function;

function ElementsUpToFrobenius(F)
  orbits := {{Frobenius(x,i) : i in [1..Degree(F)]}: x in F};
  return [Random(orbit) : orbit in orbits];
end function;

function EllipticCurvesOverField(F)
  return &cat[Twists(EllipticCurveFromjInvariant(j)) : j in F];
end function;

function EllipticCurvesOverFieldUpToFrobenius(F)
  return &cat[Twists(EllipticCurveFromjInvariant(j)) : j in ElementsUpToFrobenius(F)];
end function;

function EllipticCurvesWithPointOverFieldUpToFrobeniusAndDiamond(p,i,N)
  assert IsSquarefree(N);
  F := GF(p,i);
  ECs := [E: E in EllipticCurvesOverFieldUpToFrobenius(F) | (#E(F) mod N) eq 0];
  ECs_with_point := [* *];
  for E in ECs do
    gens := [P*(Order(P) div N) : P in Generators(E(F)) | (Order(P) mod N) eq 0];
    assert #gens eq 1;
    Append(~ECs_with_point,<E,gens[1]>);
  end for;
  return ECs_with_point;
end function;


function NonCuspidalPlacesUpToDiamond(C,i,N)
  assert IsSquarefree(N);
  p := Characteristic(BaseRing(C));
  F := GF(p,i);
  ECs := [E: E in EllipticCurvesOverFieldUpToFrobenius(F) | (#E(F) mod N) eq 0];
  places := [];
  for E in ECs do
    gens := [P*(Order(P) div N) : P in Generators(E(F)) | (Order(P) mod N) eq 0];
    assert #gens eq 1;
    time Append(~places,X1PlaceFromEllipticCurve(C,E,gens[1]));
  end for;
  return places;
end function;

