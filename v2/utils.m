intrinsic MDAttributeAsString(S, name::MonStgElt, level::MonStgElt) -> MonStgElt
{ Returns Sprint(S"name, level) if S"name is assigned and "<unassigned>" otherwise  }
  // needs new magma version
  //if assigned S"name then
  //   return Sprint(S"name, level);
  //end if;
  return "<unassigned>";
end intrinsic;