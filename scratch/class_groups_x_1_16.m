_<x>:=PolynomialRing(Integers());
f := x*(x^2+1)*(x^2+2*x-1);
C := HyperellipticCurve(f);

for t in [-100..100] do
  d := Evaluate(f,t);
  if d ge 0 then continue; end if;
  d_squarefree := - &* [pe[1] : pe in Factorization(d) | pe[2] mod 2 eq 1];
  d_square := &* [pe[1]^(pe[2] div 2) : pe in Factorization(d)];
  
  K<a0> := QuadraticField(Evaluate(f,t));
  a := a0*d_square;
  //print d,d_squarefree,K;

  OK := MaximalOrder(K);

  I1 := ideal< OK | [t-1,a-2]>;
  I2 := ideal< OK | [t-1,a+2]>;

  g := ((-6*t^2 - 4*t + 2)*a + (t^5 + 13*t^4 -  2*t^3 + 10*t^2 - 7*t + 1))/(t^5 - 5*t^4 + 10*t^3 - 10*t^2 + 5*t - 1);
  I := I1/I2;
  factors := [[Norm(pe[1]),pe[2]] : pe in Factorization(I^5*(1/g))];
  if #factors gt 0 then
    print [[Norm(pe[1]),pe[2]] : pe in Factorization(OK*g)];
    print [[Norm(pe[1]),pe[2]] : pe in Factorization(I)];
    print [[Norm(pe[1]),pe[2]] : pe in Factorization(I^5*(1/g))];
  end if;
  is_principal := IsPrincipal(I);
  print t,d_squarefree, is_principal, IsPrincipal(I^5);//,Invariants(ClassGroup(K));
  if not IsPrincipal(I^5) then
    print d,d_squarefree,a^2 eq d,K;
  end if;


end for;


C:=HyperellipticCurve(f);
//FF<X,Y> := FunctionField(C);
FF<x,y> := FunctionField(C);
p1,p2 := Explode(Zeros(x-1)); p1,p2;
IsPrincipal(5*(p1-p2));
t := -10;
K<a> := QuadraticField(Evaluate(f,t));
OK := MaximalOrder(K);

I1 := ideal< OK | [t-1,a-2]>;
I2 := ideal< OK | [t-1,a+2]>;

_<x>:=PolynomialRing(Integers());
f := x*(x^2+1)*(x^2+2*x-1);

for t in [1..100] do
  if not IsSquarefree(t) then continue; end if;
  K := QuadraticField(Evaluate(f,t));
  print t, NarrowClassNumber(K)/ClassNumber(K);
end for;

A3<x,y,z>:=AffineSpace(Rationals(),3);

C := Curve(A3,[y^2-(x*(x^2+1)*(x^2+2*x-1)), z^5*(x^5 - 5*x^4 + 10*x^3 - 10*x^2 + 5*x - 1)-((-6*x^2 - 4*x + 2)*y + x^5 + 13*x^4 - 2*x^3 + 10*x^2 - 7*x + 1)]);

C1 := ProjectiveClosure(C);
G := AutomorphismGroup(C1);
auts := Automorphisms(C1);
print auts[2];
C2 := CurveQuotient(AutomorphismGroup(C1, [G ! auts[2]]));

print C2;

_<x> :=PolynomialRing(Integers());
f := x^6 - x^5 + 5*x^3 - x + 1;
C2 := HyperellipticCurve(f);
K := QuadraticField(Discriminant(f));
//f1,f2 := Explode([pe[1] : pe in Factorization(ChangeRing(f,K))]);
//L := SplittingField(f);
//fs :=  &cat [[pe[1] : pe in Factorization(ChangeRing(f,M[1])) |  Degree(pe[1]) eq 4] : M in Subfields(L)];

L1 :=  NumberField(x^3 - x^2 + 2*x + 2);
fs := [pe[1] : pe in Factorization(ChangeRing(f,L1)) |  Degree(pe[1]) eq 4];

C := C2;
  time Hk , AtoHk, expvecs, factorbase := TwoCoverDescent(C: Raw:=true);
  
  // Get the map that allows us to construct covers explicitly as elements in an algebra A
  HktoA := createHktoA(AtoHk, expvecs, factorbase);

  Hk := Setseq(Hk);
  E := EllipticChabauty(C, fs[1], Hk[1], HktoA);




  time Hk , AtoHk, expvecs, factorbase := TwoCoverDescent(C: Raw:=true);
  
  // Get the map that allows us to construct covers explicitly as elements in an algebra A
  HktoA := createHktoA(AtoHk, expvecs, factorbase);

  Hk := Setseq(Hk);
  points := {@ @};
  hk_info := [* *];
  for i in [1..#Hk] do
    hk := Hk[i];
    print "doing hk", i, hk;
    assert AtoHk(HktoA(hk)) eq hk;
    success := false;
    fE_info := [* *];
    
    // we have 8 Galois conjugacy classes of genus 1 curves to work with
    for fE in fs do
      print "======== doing f =====", fE;
      success, gamma_g, new_points, message := EllipticChabauty(C, fE, hk, HktoA);
      print "chabauty result", success, gamma_g, new_points, message;
      Append(~fE_info, <gamma_g,message>);
      if success then
        points := points join new_points;
        break;
      end if;
    end for;
    //Append(~hk_info,fE_info);
    if not success then
      print "=============================== failed ======================";
      print hk, i, #Hk;
      print hk_info;
      //break;
    end if;
  end for;
  if #points gt 2 then
    print "====================== extra points!!! ====================";
    print points;
  end if;


  time Hk , AtoHk, expvecs, factorbase := TwoCoverDescent(C: Raw:=true);
  
  // Get the map that allows us to construct covers explicitly as elements in an algebra A
  HktoA := createHktoA(AtoHk, expvecs, factorbase);

  Hk := Setseq(Hk);
  points := {@ @};
  hk_info := [* *];
    i := 2;
    hk := Hk[i];
    print "doing hk", i, hk, HktoA(hk);
    assert AtoHk(HktoA(hk)) eq hk;
    success := false;
    fE_info := [* *];
    

    for fE in [* f1, f2 *] do
      print "======== doing f =====", fE;
      success, gamma_g, new_points, message := EllipticChabauty(C, fE, hk, HktoA);
      print "chabauty result", success, gamma_g, new_points, message;
      Append(~fE_info, <gamma_g,message>);
      if success then
        points := points join new_points;
        break;
      end if;
    end for;
    //Append(~hk_info,fE_info);
    if not success then
      print "=============================== failed ======================";
      print hk, i, #Hk;
      print hk_info;
    end if;
