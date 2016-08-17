function Counts(list);
//    {returns a dictionary containing the counts of all elements in the list
//    
//    Input: a list
//    
//    Output: a dictionay d such that d[i] is equal to the number of times that i occors in the list
//    }
    if #list eq 0 then;
        return AssociativeArray();
    end if;
    
    list := Sort(list);
    counts := AssociativeArray(); 
    old := list[1];
    counts[old]=1;
    list := Remove(list,1);
    for i in list do;
        if i eq old then;
            counts[i] +:= 1;
        else;
            counts[i] := 1;
            old := i;
        end if;
    end for;
    return counts;
end function;



function DegreeTypes_of_Degree(degree,curve)
//{A degree type of degree d is a list of pairs of integers [(n1,d1),...,(nk,dk)] such that
//The sum n1*d1+...+nk*dk is d. The values of di are restriced to the numbers for which the curve
//has a place of that degree. And the tuples are also sorted such that di >= d(i+1) and if di = d(i+1) then ni >= n(i+1).
//This function returns all degree types satisfying the above restrictions.
//}
    occurring_degrees := [i : i in [1..degree] | HasPlace(curve,i) ];
    degree_types_old:= [<[<0,0>],degree>];
    degree_types_new:= [<[<0,0>],degree>];
    degree_types_done := [];
    for d in Reverse(occurring_degrees) do
        for n in Reverse([1..Floor(degree/d)]) do
            for degree_type in degree_types_old do
                for i in [1..Floor(degree_type[2]/(n*d))] do
                    d_t := <degree_type[1] cat [<n,d> : j in [1..i]],degree_type[2]-n*d*i>;
                    degree_types_new := Append(degree_types_new, d_t);
                end for;
            end for;
            degree_types_done := degree_types_done cat [d_t[1][2..#d_t[1]] : d_t in degree_types_new | d_t[2] eq 0];
            degree_types_new := [d_t : d_t in degree_types_new | d_t[2] gt 0];
            degree_types_old := degree_types_new;
        end for;
    end for;
    return degree_types_done;
end function;



function Divisors_of_DegreeType(degree_type,curve)
    divisors_old:={<DivisorGroup(curve) ! 0,[]>};
    divisors_new:=divisors_old;
    for d in degree_type do;
        divisors_new:={<D1[1]+d[1]*D2,Append(D1[2],D2)> : D1 in divisors_old, D2 in Places(curve,d[2]) | D2 notin D1[2]};
        //divisors_new:={D1+D2 : D1 in divisors_old, D2 in Places(curve,d)};
        divisors_old:=divisors_new;
    end for;
    return divisors_new;
end function;

//function DominatingDegreeTypes_naive()
//{Returns 
//}
//end function;


function Gonality_lowerbound(curve,bound : verbose:=false)
//{Computes the gonality of a curve.
// Input: curve - a projective curve over a finite field
//        bound - an integer
//
// Output: True,bound - if the gonality of the curve is >= bound,
//         False, gon - where gon is the gonality of the curve otherwise
//
// Note this is horribly slow, so it only works in practice over very small finite fields and very small gonalities.
//}
    for degree in [1..bound-1] do;
        if verbose then;
            print "Checking divisors of degree:",degree;
        end if;
        for degree_type in DegreeTypes_of_Degree(degree,curve) do;
            for divisor in Divisors_of_DegreeType(degree_type,curve) do;
                if Dimension(divisor[1]) gt 1 then;
                    return false,degree;
                end if;
            end for;
        end for;
    end for;
    return true,bound;
end function;


function Gonality_naive(curve : verbose := false)
//{Computes the gonality of a curve.
// Input: a projective curve over a finite field
// Output: the gonality
//
// Note this is horribly slow, so it only works in practice over very small finite fields and very small gonalities.
//}
    dummy,gonality:=Gonality_lowerbound(curve,2*Genus(curve)+4 : verbose:=verbose);
    return gonality;
end function;


function Gonality(curve : search_bound := 129, gonality_bound := 0, verbose := false, fall_back_to_naive := true)
//{Computes the gonality of a curve.
// Input: a projective curve over a finite field
// Output: the gonality
//
// Note this is slow, so it only works in practice over very small finite fields and reasonably small gonalities.
//}
    Fp := BaseRing(curve);
    p := Characteristic(Fp);
    g := Genus(curve);
    plc1 := Places(curve,1);
    sum_plc1 := &+ plc1;
    n := Ceiling(#plc1/(p+1));
    if verbose then
        print "p,#places,#places/(p+1)",p,#plc1,n;
    end if;
    if n lt 4 and fall_back_to_naive then
       if verbose then
           print "falling back to naive algorithm";
       end if;
       return Gonality_naive(curve : verbose :=verbose );
    end if;
    for degree in [0..2*g+1] do;
        if degree+n eq gonality_bound then
            return degree+n;
        end if;
        if verbose then
            print "Checking if there are functions of degree",degree+n;
        end if;
        for degree_type in DegreeTypes_of_Degree(degree,curve) do;
            for divisor in Divisors_of_DegreeType(degree_type,curve) do;
                divisor2 := divisor[1] + sum_plc1;
                H,m:=RiemannRochSpace(divisor2);
                if p^Dimension(H) gt search_bound then
                    return "fail";
                end if;
                if Dimension(H) gt 1 then
                    d := Min(FunctionDegrees(divisor2));
                    if d eq degree+n then
                        return d;
                    end if;
                    assert d gt degree+n;
                end if;
            end for;
        end for;
    end for;
end function;
