X := MDX11(2, 20, GF(7));

procedure test_print()
  TSTAssertEQ(
    Sprint(X),
    "The modular curve X_1(2, 20) over Finite field of size 7"
  );
end procedure;

procedure test_level()
  TSTAssertEQ(
    Level(X),
    20
  );
end procedure;

test_print();
test_level();