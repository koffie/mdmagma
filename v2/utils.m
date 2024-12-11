intrinsic MDAttributeAsString(S, name::MonStgElt, level::MonStgElt) -> MonStgElt
{ Returns Sprint(S"name, level) if S"name is assigned and "<unassigned>" otherwise  }
  // needs new magma version
  //if assigned S"name then
  //   return Sprint(S"name, level);
  //end if;
  return "<unassigned>";
end intrinsic;


intrinsic MDValues(A::Assoc) -> SeqEnum
{ Return the values of the associative array }
    return [A[key] : key in Keys(A)];
end intrinsic;


intrinsic MDMagmaSourceDir() -> MonStgElt
{ Returns the source directory where mdmagma is installed }
  filenames := GetFilenames(MDAttributeAsString);
  assert #filenames eq 1;
  source_dir := "/" cat Join(s[1..(#s - 1)], "/") where s := Split(filenames[1,1],"/");
  return source_dir;
end intrinsic;

intrinsic MDMultiset(S::Setq) -> SetMulti
{The multiset derived from S}
    if IsNull(S) then
        return {* *};
    end if;
    MS := {* Universe(S) | *};
    // the reason for this is that magma sometimes screws up hashing
    // so we need to explicitly compare
    for x in S do
      added := false;
      for y in MS do
        if x eq y then;
          Include(~MS, y);
          added := true;
          break;
        end if;
      end for;
      if not added then;
         Include(~MS, x);
      end if;
    end for;
    return MS;
end intrinsic;
