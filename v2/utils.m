intrinsic MDAttributeAsString(S, name::MonStgElt, level::MonStgElt) -> MonStgElt
{ Returns Sprint(S"name, level) if S"name is assigned and "<unassigned>" otherwise  }
  // needs new magma version
  //if assigned S"name then
  //   return Sprint(S"name, level);
  //end if;
  return "<unassigned>";
end intrinsic;

intrinsic MDMagmaSourceDir() -> MonStgElt
{ Returns the source directory where mdmagma is installed }
  filenames := GetFilenames(MDAttributeAsString);
  assert #filenames eq 1;
  source_dir := "/" cat Join(s[1..(#s - 1)], "/") where s := Split(filenames[1,1],"/");
  return source_dir;
end intrinsic;