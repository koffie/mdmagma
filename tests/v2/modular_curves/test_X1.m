procedure test_print()
  TSTAssertEQ(
    Sprint(MDX1(10,GF(3))),
    "The modular curve X_1(10) over Finite field of size 3"
  );
end procedure;

test_print();