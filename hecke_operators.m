Attach("assorted_functions.m");
import "assorted_functions.m": Functions_xyrsbcF2F3;

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

function pIsogeniesFiniteField1(E,p);
//Returns the p+1 elliptic curves that are p-isogenous to an initial E
//curve defined over a finite field
//Raises an error in characteristic 2,3 or if E is supersingular 
    Fq := BaseRing(E);
    jE := jInvariant(E);
    fp := ClassicalModularPolynomial(p);
    fp := ChangeRing(Parent(fp),Fq) ! fp;
    fpjE := UnivariatePolynomial(Evaluate(fp,2,jE));
    Fqn := SplittingField(fpjE);
    
    A2 := AffineSpace(Fqn,2);
    X0p := ModularCurve(A2,"Canonical",p);
    mp := ModuliPoints(X0p,BaseChange(E,Fqn));
    assert #mp eq p+1;
    return [Isogeny(BaseChange(E,Fqn),x) : x in mp];
end function;

function pIsogeniesFiniteField(E,p);
//Returns the p+1 elliptic curves that are p-isogenous to an initial E
//curve defined over a finite field
    assert IsPrime(p);
    //assert p ne 2; 
    Fq := BaseRing(E);
    fp := DivisionPolynomial(E,p);
    Fqn := SplittingField(fp);
    Fq2n := RandomExtension(Fqn, 2);
    xi := Roots(ChangeRing(fp,Fq2n));
    Eq2n := BaseChange(E,Fq2n);
    fE := DefiningEquation(Eq2n);
    yi := [Roots(UnivariatePolynomial(Evaluate(Evaluate(fE,3,1),1,x[1])))[1][1] : x in xi];
    Pi := [E(Fq2n) ! [xi[i][1],yi[i]] : i in [1..#xi]];
    //now Pi is a list with all torsion points
    
    R<X> := PolynomialRing(Fq2n);
    kernel_polynomials := {&*[X-(j*P)[1] : j in [1..Ceiling((p-1)/2)]] : P in Pi};
    isogenies := [PowerStructure(MapSch) | ];
    for f in kernel_polynomials do;
        Ef,phi := IsogenyFromKernel(Eq2n,f);

        Append(~isogenies,phi);
    end for;
    return isogenies;
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


function MyMultiset(itterable);
    uniques := [];
    for i in itterable do;
        if i notin uniques then;
            Append(~uniques,i);
        end if;
    end for;
    
    return [<j,#[1 : i in itterable | i eq j]> : j in uniques];
end function;

function Tp_X1N_noncuspidal_place(P,p);
    assert IsPrime(p);
    ZZ := IntegerRing();
    E := EllipticCurveFromX1Place(P);
    X1N := Curve(P);
    assert Characteristic(BaseRing(X1N)) ne p;
    isogenies := pIsogeniesFiniteField(E,p);
    Eq2n := Domain(isogenies[1]);
    Fq2n := BaseRing(Eq2n);

    Pi := [<Codomain(phi),phi(Eq2n ! [0,0])> : phi in isogenies];
    xyi := [TateNormalForm_xy(P[1],P[2]) : P in Pi];
    places := [Places(X1N(Fq2n) ! xy) : xy in xyi];
    assert &and[#p eq 1 :  p in places];
    places := MyMultiset(&cat places);
    d := Degree(P);
    return &+[ (ZZ ! (place[2]*d/Degree(place[1])))*place[1] : place in places]; 
end function;

function Tp_X1N_noncuspidal(D,p);
  P,e := Support(D);
  return &+[e[i]*Tp_X1N_noncuspidal_place(P[i],p) : i in [1..#P]]; 
end function;

function diamond_operator_X1N_noncuspidal_place(P,d);
    E := EllipticCurveFromX1Place(P);
    X1N := Curve(P);
    Fq := BaseRing(E);
    xy := TateNormalForm_xy(E,d*(E ! [0,0]));
    dP := Places(X1N(Fq) ! xy);
    assert #dP eq 1;
    return dP[1];
end function;

function diamond_operator_X1N_noncuspidal(D,p);
  P,e := Support(D);
  return &+[e[i]*diamond_operator_X1N_noncuspidal_place(P[i],p) : i in [1..#P]]; 
end function;


function Tp_pdp_1_noncuspidal_place(P,p);
    return Tp_X1N_noncuspidal_place(P,p)-p*diamond_operator_X1N_noncuspidal_place(P,p)-P;
end function;

function Tp_pdp_1_noncuspidal(D,p);
  P,e := Support(D);
  return &+[e[i]*Tp_pdp_1_noncuspidal_place(P[i],p) : i in [1..#P]]; 
end function;


function PositiveRankHeckePolynomial(S,n);
//Returns the characteristic polynomial of the hecke operator n on the subspace
//of the cuspidal modular symbol spaces S corresponding to the part where the LRatio is 0
//Under BSD this is exactly the part corresponding to the part of S where the corresponding abelian variety has positive rank
 return &*[HeckePolynomial(Si,n) : Si in NewformDecomposition(S) | LRatio(AssociatedNewSpace(Si),1) eq 0];
end function;



function PositiveRankHeckePolynomialX1N(N,n,chars);
//The input space needs to be cuspidal of sign 0
//Returns the characteristic polynomial of the hecke operator n on the subspace
//of the cuspidal modular symbol spaces S corresponding to the part where the LRatio is 0
//with respect to at least one of the characters in chars
//Under BSD this is exactly the part corresponding to the part of S where the corresponding abelian variety when twisted by one of the characters has positive rank

  D := FullDirichletGroup(N);
  chars := [D ! chi : chi in chars];
  ann_pol := 1;
  for d in Elements(D) do;
    M := ModularSymbols(d,2,0);
    S := CuspidalSubspace(M);
    for Si in NewformDecomposition(S) do;
      Snew := AssociatedNewSpace(Si);
      rank_0 := &or[Dimension(Snew)/2 ne Dimension(TwistedWindingSubmodule(Snew,1,chi)) : chi in chars];
      if rank_0 then;
        ann_pol := ann_pol*Sqrt(HeckePolynomial(Si,n));
      end if;
    end for;
  end for;
  return ann_pol;
end function;




    