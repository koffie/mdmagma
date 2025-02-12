X := MDX0(38,GF(3));

procedure test_print()
  TSTAssertEQ(
    Sprint(X),
    "The modular curve X_0(38) over Finite field of size 3"
  );
end procedure;

procedure test_Level()
  TSTAssertEQ(Level(X), 38);
end procedure;

procedure test_Genus()
  TSTAssertEQ(Genus(X), 5);
end procedure;

test_print();
test_Level();
//test_Genus();