function NextPrimeInArithmeticProgression(p, a, b);
// the smallest prime q > p such that q mod a = b mod a
  assert GCD(a,b) eq 1;
  b := b mod a;
  q := p;
  repeat q := NextPrime(q); until q mod a eq b;
  return q;
end function;

function UnitExponent(N)
// the exponent of the group (Z/NZ)^*
  return Exponent(UnitGroup(Integers(N)));
end function;

function SuitableRootOfUnity(N : lowerbound := 0)
//
  r := UnitExponent(N);
  if lowerbound eq 0 then
    lowerbound := r;
  end if;
  p := NextPrimeInArithmeticProgression(lowerbound, r, 1);
  zeta := RootOfUnity(r, GF(p));
  return zeta, r;
end function;


function NextSuitableRootOfUnity(zeta, r)
//
  p := Characteristic(Parent(zeta));
  p := NextPrimeInArithmeticProgression(p, r, 1);
  zeta := RootOfUnity(r, GF(p));
  return zeta, r;
end function;

function IsStrange(label : zeta := 0, r := 0, new := false)
    chi := DirichletCharacter(label: zeta:=zeta, r:=r);
    S := CuspidalSubspace(ModularSymbols(chi, 2, 1));
    if new then
      S := NewSubspace(S);
    end if;
    if Dimension(S) eq 0 then return false, 0, 0; end if;
    W := WindingSubmodule(S);
    dimS := Dimension(S);
    dimW := Dimension(W);
    return dimS ne dimW, dimW, dimS;
end function;



function StrangeInfo(N : lowerbound := 2^42, tries := 10, new := false, proof := true)
    zeta, r := SuitableRootOfUnity(N : lowerbound := lowerbound);
    potentially_strange_characters := [];
    for label in ConreyCharacterOrbitReps(N) do
      // skip the trivial character
      if Order(label) eq 1 then continue; end if;

      success := false;
      for i in [1..tries] do
        try
          is_strange, dimW, dimS := IsStrange(label : zeta:=zeta, r:=r, new := new);
          success := true;
          break;
        catch e
          zeta, r := NextSuitableRootOfUnity(zeta, r);
        end try;
      end for;
      error if not success, "computation did not succeed, increase the number of tries";


      if is_strange and proof then
        is_strange, dimW, dimS := IsStrange(label : zeta:=0, r:=0, new := new);
      end if;

      if is_strange then
        Append(~potentially_strange_characters, <CharacterOrbitLabel(label), Order(label), dimS - dimW, dimS>);
      end if;

    end for;
    return potentially_strange_characters;
end function;

function Quote(o)
  return Sprintf("\"%o\"",o);
end function;

function AsListOfStrings(l)
  return [Quote(x) : x in l];
end function;

function PrintStrangeInfoRange(a, b : lowerbound := 2^42, tries := 10, new := false, proof := true)
  print "{";
  separator := ",";
  for N in [a..b] do
    strange_info := StrangeInfo(N : lowerbound := lowerbound, tries := tries, new := new, proof := proof);
    if N eq b then separator := ""; end if;
    print Quote(N), ":", [AsListOfStrings(l) : l in strange_info], separator;
  end for;
  print "}";
  return "";
end function;

if assigned start and assigned stop then
  start := StringToInteger(start);
  stop := StringToInteger(stop);
  if not assigned new then new := false; end if;
  PrintStrangeInfoRange(start, stop : new := new);
  exit;
end if;

