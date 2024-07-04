// The functions here make it easy to combine two cover descent on hyper elliptic curves
// together with elliptic curve Chabauty.

QQ := Rationals();

// Innocent helper function
function Zip(a, b)
  //similar to pythons zip function
  assert #a eq #b;
  n := #a;
  return [<a[i],b[i]> : i in [1..n]];
end function;

// Magma's TwoCoverDescent function returns a function AtoHk that maps elements from
// some algebra A to the fake 2-Selmer set Hk. This function returns the inverse of this
// map. Elements in the algebra A can be used to construct the two-cover of the
// corresponding element in the fake 2-Selmer set.
function createHktoA(AtoHk, expvecs, factorbase)
  A := Domain(AtoHk);
  A2 := Parent(factorbase[1]);
  phi := hom< A2 -> A | A.1 >;
  Hk := Domain(expvecs);
  HktoA := map< Hk -> A | h :-> phi(&*[ef[2]^ef[1] : ef in Zip(Eltseq(expvecs(h)),factorbase)])>;
  return HktoA;
end function;

// This function carries out the elliptic curve chabauty algorithm explained in Section
// 4.2 of the paper. It is inspired by the example shown in:
//   https://magma.maths.usyd.edu.au/magma/handbook/text/1566#18074
// If C is a hyper elliptic curve given by y^2=f(x) then g should be a factor of degree
// 3 or 4 of f. The factor g is allowed to be defined over a field extension.
// The element hk should be an element of the fake 2 selmer group. And corresponds to
// a two cover D -> C . D will map to the genus 1 curve gamma*y^2=g(x) for some
// twisting factor gamma. This function applies elliptic curve chabauty to find all
// points on C(Q) that come from a point on D(Q). The return value is a 4 tuple:
//
//    (success, gamma, points, message)
//
//  - success: a boolean that when true guarantees that all points on C(Q) coming from
//    D(Q) are contained in the set of points output by this function
//  - gamma: is the twisting factor corresponding to g and hk ensuring that D has a map to
//    gamma*y^2=g(x)
//  - points: a set of points on C
//  - message: an message that indicates why we failed to do elliptic curve chabauty in
//    case of failure
function EllipticChabauty(C, g, hk, HktoA : Bound:= 100)
  // Basic setup is as follows
  assert IsMonic(g);
  assert BaseRing(C) eq QQ;
  assert Degree(g) in [3,4];
  success := true;
  points := {@ @};
  P1:=ProjectiveSpace(QQ,1);
  CtoP1 := map< C -> P1 | [C.1,C.3] >;
  K := BaseRing(g);
  A := Codomain(HktoA);
  f := ChangeRing(Modulus(A), QQ);
  L := quo< Parent(g) | g>;
  assert Evaluate(f, L.1) eq 0;
  AtoL := hom< A -> L | L.1>;

  // For some gamma, we have that a 2-cover hk covers the genus 1 curve E:y^2=gamma * g.
  // The following obtains this gamma.

  gamma := Norm(AtoL(HktoA(hk)));
  for p in PrimeDivisors(Numerator(Norm(gamma))) do
    while IsIntegral(gamma/p^2) do
      gamma := gamma/p^2;
    end while;
  end for;
  print "gamma", gamma;
  gamma_g := <gamma,g>;

  // the promised elliptic curve above
  E := HyperellipticCurve(gamma*g);

  // However, we need it as an EllipticCurve object, so we do the following
  EtoP1:=map<E->P1|[E.1,E.3]>;
  iselliptic, E1raw, EtoE1raw, E1rawtoE := IsEllipticCurve(E);

  // If magma didn't find an elliptic curve model, we manually search for a rational
  // point to function as the point at infinity
  if not iselliptic then
    print "finding points";
    time Epoints := Points(E : Bound:=Bound);

    // If we can't find a point, we report failure. Theoretically we could still try
    // to prove that E is a pointless genus 1 curve. However this will be hard since
    // the two descent we did earlier implies that E is coverd by an everywhere
    // locally solvable curve. And hence there are no local obstructions.
    if #Epoints eq 0 then
      success := false;
      return success, gamma_g, points, "point search failed";
    end if;

    iselliptic := true;
    P0 := Setseq(Epoints)[1];
    E1raw,EtoE1raw:=EllipticCurve(E,P0);
    E1rawtoE := Inverse(EtoE1raw);
  end if;

  // Just thing of the following E1 as E but with a nicer model
  E1 := MinimalModel(E1raw);

  print "Elliptic curve", E1;

  // We need to compute the MW group of the elliptic curve. Return if that fails.
  time MW1, MW1toE1set, flag1, flag2 := MordellWeilGroup(E1);


  if not (flag1 and flag2) then;
    r := Degree(K);
    return success, gamma, points, <"mw failed",Invariants(MW1), r>;
  end if;

  // Some functions to make certain objects more explicit
  MW1toE1 := map<MW1 -> E1 | x :-> MW1toE1set(x)>;
  E1toE1raw := Isomorphism(E1,E1raw);
  E1toE := E1toE1raw*E1rawtoE;
  E1toP1 := Expand(E1toE1raw*E1rawtoE*EtoP1);

  // The Chabauty rank condition for elliptic curve Chabauty to work
  if TorsionFreeRank(MW1) ge Degree(K) then
    success := false;
    return success, gamma_g, points, <"rank to large",TorsionFreeRank(MW1),Degree(K)>;
  end if;

  // If we got here, then we can do elliptic curve Chabauty
  print "starting chabauty";
  pointset, R := Chabauty(MW1toE1, E1toP1);
  print R, #pointset;


  // At this point the strategy is successful, so we can determine all points on the twist
  for t in pointset do
    P := E1toE(MW1toE1(t));
    xP := EtoP1(P);
    xP := P1(QQ) ! xP;
    new_points := RationalPoints( xP@@CtoP1 );
    points := points join new_points;
  end for;

  return success, gamma_g, points, "found all points";
end function;