Attach("X1_N_equations.m");
import "X1_N_equations.m": equations, gonality_upperbound;


function FunctionDegrees(divisor)
//{On input a divisor returns a list containing the degrees of all the non constant 
//functions in the RiemannRochspace of the divisor
//}
    constantField := BaseRing(Curve(divisor));
    space,map := RiemannRochSpace(divisor);
    return [Degree(map(i)) : i in space | map(i) notin constantField];
end function;

function SmallestCoprimePrime(N);
//{Returns the smallest prime that is coprime to the input integer.}
    for i in [1..N] do;
        p := NthPrime(i);
        if not N mod p eq 0 then;
            return p;
        end if;
    end for;
end function;


function DivisorSumsOfDegreeType(degree_type,divisors_by_degree)
    divisors := &cat divisors_by_degree;
    ZZdiv := FreeAbelianGroup(#divisors);
    divisors_of_degtype := {ZZdiv ! 0};
    for degree in degree_type do;
        divisors_of_degtype := {d+ZZdiv.Index(divisors,i) : 
                               d in divisors_of_degtype, i in divisors_by_degree[degree]};
    end for;
    return divisors_of_degtype; 
end function;


function X_1(N,base_ring)
    //Input: N - an integer
    //       base_ring - a ring
    //Output: C - a curve
    //Returns an algebraic model C of the modular curve X_1(N) as a curve over base_ring
    //The model is such that FunctionField(C).1 = x and FunctionField(C).2 = y
    //where x,y are as in section 2.1 of http://arxiv.org/pdf/1307.5719v1.pdf
    C:=Curve(Spec(ChangeRing(Parent(equations[N]),base_ring)),equations[N]);
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

function Cusp_GalQ_orbit_to_signature(cusp);
    //Input: cusp - a Place of X_1(N) that is a cusp
    //Output: [n1,n2] - integers n1 and n2 are the valuations of F2 and F3 resp
    //[n1,n2] is called the signature of cusp since it uniquely determines the galois orbit
    //over Q of (a lift of) this cusp. Note that this function still works if cusps is defined over a field of characteristic > 0. 
    C := Curve(cusp);
    x,y,r,s,b,c,F2,F3 := Functions_xyrsbcF2F3(C);
    return [Valuation(F2,Support(1*cusp)[1]),Valuation(F3,Support(1*cusp)[1])];
end function;

function Cusps_X1(curve)
    //Input: curve - the modular curve X_1(N) as returned by the function X_1
    //Output: places - a Set containing all places of X_1(N) that are cusps
    x,y,r,s,b,c,F2,F3 := Functions_xyrsbcF2F3(curve);
    return SequenceToSet(Support(Divisor(F2)) cat Support(Divisor(F3)) cat Support(Divisor(x)) cat Support(Divisor(x)));
end function;

function IsCusp(P);
//Given a point on X_1(N) returns wether this point is a cusp
    X1N := Curve(P);
    x,y,r,s,b,c,F2,F3:=Functions_xyrsbcF2F3(X1N);
    if Valuation(b,P) lt 0 or Valuation(c,P) lt 0 then;
      return true;
    end if;
    bP:=Evaluate(b,P);
    cP:=Evaluate(c,P);
    return bP^3 * (cP^4 - 8*bP*cP^2 - 3*cP^3 + 16*bP^2 - 20*bP*cP + 3*cP^2 + bP - cP) eq 0;
end function;

function Cusp_GalQ_orbits_dict(curve,cusps)
    //Input: curve - the modular curve X_1(N) as returned by the function X_1
    //       cusps - the set of cusps as returned by Cusps_X1
    //Output: cusp_dict - a dictionary/AssociativeArray that contains
    x,y,r,s,b,c,F2,F3 := Functions_xyrsbcF2F3(curve);
    cusp_dict := AssociativeArray(Universe([[1]]));
    for cusp in cusps do;
        signature := [Valuation(F2,cusp),Valuation(F3,cusp)];
        if IsDefined(cusp_dict,signature) then;
            Append(~cusp_dict[signature],cusp);
        else
            cusp_dict[signature] := [cusp];
        end if;
    end for;
    return cusp_dict;
end function;

function Cusp_GalQ_orbits(curve,cusps)
    cusp_dict := Cusp_GalQ_orbits_dict(curve,cusps);
    return [1*(&+cusp_dict[i]) : i in Keys(cusp_dict)];
end function;

function ClassSubgroup(given_divisors,Grp,m1,m2);    
    M := FreeAbelianGroup( #given_divisors );
    f := hom< M -> Grp | [m2(i) : i in given_divisors] >;
    N := Kernel(f);
    N := Lattice( Matrix( [ Eltseq(M ! i) : i in Generators(N)] ) );
    subgroup := Image(f);
    return N,subgroup,f;
end function;

function DivisorAsSumOfGivenDivisors(D,N,subgroup,f,m2);
    //Grp,m1,m2:=ClassGroup(Curve(D));
    if m2(D) notin subgroup then;
        return [0];
    end if;
    fInvD:=Vector(Eltseq((f^(-1))(m2(D))));
    v:=ClosestVectors(N,Vector(fInvD) : Max:=1)[1];
    fInvD:=Eltseq(fInvD-v); //now fInvD should be shorter, making compuations faster
    return fInvD;
end function;

function DivisorsOfDegreeType(degree_type,curve);
    divisors_old:={DivisorGroup(curve) ! 0};
    divisors_new:=divisors_old;
    for d in degree_type do;
        divisors_new:={D1+D2 : D1 in divisors_old, D2 in Places(curve,d)};
        divisors_old:=divisors_new;
    end for;
    return divisors_new;
end function;





function DegreeTypesNonCuspidal_of_Degree(degree,cusps,curve)
    occurring_degrees := {i : i in [1..degree] | 
                    HasPlace(curve,i) and not &and[i in cusps : i in Places(curve,degree)] };
    return RestrictedPartitions(degree,occurring_degrees);
end function;

function DegreeTypesNonCuspidal_up_to_Degree(degree,cusps,curve)
    return &cat[DegreeTypesNonCuspidal_of_Degree(d,cusps,curve) : d in [1..degree]];
end function;

function DegreeTypesCuspidal_of_Degree(degree,cusps)
    occurring_degrees := {Degree(c) : c in cusps};
    return RestrictedPartitions(degree,occurring_degrees);
end function;

function DegreeTypesCuspidal_up_to_Degree(degree,cusps)
    return &cat[DegreeTypesCuspidal_of_Degree(d,cusps) : d in [1..degree]];
end function;

function NonCuspidalPlaces(degree,curve)
    return [p : p in Places(curve,degree) | not IsCusp(p)];
end function;

function CuspidalPlaces(degree,cusps)
    return [p : p in cusps | Degree(p) eq degree];
end function;

function NonCuspidalDivisorsOfDegreeType(degree_type,cusps,curve)
    divisors_old:={DivisorGroup(curve) ! 0};
    divisors_new:=divisors_old;
    for d in degree_type do;
        divisors_new:={D1+D2 : D1 in divisors_old, D2 in NonCuspidalPlaces(d,curve)};
        divisors_old:=divisors_new;
    end for;
    return divisors_new;
end function;

function CuspidalDivisorsOfDegreeType(degree_type,cusps,curve)
    divisors_old:={DivisorGroup(curve) ! 0};
    divisors_new:=divisors_old;
    for d in degree_type do;
        divisors_new:={D1+D2 : D1 in divisors_old, D2 in CuspidalPlaces(d,cusps)};
        divisors_old:=divisors_new;
    end for;
    return divisors_new;
end function;


function NonCuspidalDivisorsOfDegree(degree,cusps,curve)
    return &join[NonCuspidalDivisorsOfDegreeType(degree_type,cusps,curve) : degree_type in DegreeTypesNonCuspidal_of_Degree(degree,cusps,curve)];
end function;

function CuspidalDivisorsOfDegree(degree,cusps,curve)
    return &join[CuspidalDivisorsOfDegreeType(degree_type,cusps,curve) : degree_type in DegreeTypesCuspidal_of_Degree(degree,cusps)];
end function;
    



function CuspSignatureMultiplicitiesToDivisor(cusp_signatures,multiplicities,orbit_dict)
    return &+[multiplicities[i]*(&+orbit_dict[cusp_signatures[i]]) : i in [1..#cusp_signatures]];
end function;

function FilterCandidates(curve,cusp_signatures,candidates)
    orbit_dict := Cusp_GalQ_orbits_dict(curve,Cusps_X1(curve));
    //print Keys(orbit_dict);
    //print cusp_signatures;
    survivors := [];
    for multiplicities in candidates do;
        D := CuspSignatureMultiplicitiesToDivisor(cusp_signatures,multiplicities,orbit_dict);
        if Dimension(D) ge 1 then;
            Append(~survivors,multiplicities);
        end if;
    end for;
    return survivors;
end function;

function Main2(N : d:=0, write_results_to_file := false);
    p:=2;
    for i in PrimesInInterval(2,20) do;
        if (N mod i) ne 0 then;
            p:=i;
            break;
        end if;
    end for;
    X1modp := X_1(N,GF(p));
    if d eq 0 then;
        gon := gonality_upperbound[N];
    else;
        gon := d;
    end if;
    Grp,m1,m2 := ClassGroup(X1modp);
    cusps := Cusps_X1(X1modp);
    cusp_orbits := Cusp_GalQ_orbits(X1modp,cusps);
    cusp_signatures := [Cusp_GalQ_orbit_to_signature(cusp) : cusp in cusp_orbits];
    if write_results_to_file then;
        PrintFile("data/cusp_signatures_" cat IntegerToString(N), cusp_signatures : Overwrite:=true);
    end if;
    non_cusps_by_degree := [NonCuspidalPlaces(d,X1modp) : d in [1..gon-1]];
    non_cusps := &cat non_cusps_by_degree;
    //print "non_cusps", #non_cusps;
    if #non_cusps eq 0 then;
        return cusp_signatures,[[] : degree in [1..gon-1]];
    end if;
    non_cusp_degrees := [Degree(P) : P in non_cusps];
    
    modular_units,CuspGrp,f := ClassSubgroup(cusp_orbits,Grp,m1,m2);
    //non_cusp_functions,NonCuspGrp,f2 := ClassSubgroup(non_cusps,Grp,m1,m2);
    M := FreeAbelianGroup( #non_cusps );
    f2 := hom< M -> Grp | [m2(i) : i in non_cusps] >;
    
    ZZnoncusps := Domain(f2);
    candidates := [];
    for degree in [1..gon-1] do;
        Append(~candidates,[]);
        degree_types := RestrictedPartitions(degree,SequenceToSet(non_cusp_degrees));
        for degree_type in degree_types do;
            for D in DivisorSumsOfDegreeType(degree_type,non_cusps_by_degree) do;
                if f2(D) in CuspGrp then;
                    fInvD:=Vector(Eltseq((f^(-1))(f2(D))));
                    v:=ClosestVectors(modular_units,Vector(fInvD) : Max:=1)[1];
                    fInvD:=Eltseq(fInvD-v); //now fInvD should be shorter, making compuations faster
                    Append(~candidates[degree],fInvD);
                end if;
            end for;
        end for;
	if write_results_to_file then;
            PrintFile("data/candidates_" cat IntegerToString(N) cat "_" cat IntegerToString(degree), candidates[degree] : Overwrite:=true);
        end if;
    end for;
    
    return cusp_signatures,candidates;
end function;
