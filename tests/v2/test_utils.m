X := New(MDCrvMod);
X`level := 1;

/*
procedure test_AttributeAsString_assigned()
  attr := MDAttributeAsString(X, "base_ring", "Default");
  TSTAssertEQ(attr, "1");
end procedure;
*/

procedure test_AttributeAsString_unassigned()
  attr := MDAttributeAsString(X, "", "Default");
  TSTAssertEQ(attr, "<unassigned>");
end procedure;

procedure test_MDMagmaSourceDir()
  source_dir := MDMagmaSourceDir();
  TSTAssertEQ(s[#s] where s := Split(source_dir,"/"), "v2");
end procedure;

//test_AttributeAsString_assigned();
test_AttributeAsString_unassigned();
test_MDMagmaSourceDir();