procedure test_print()
  TSTAssertEQ(
    Sprint(MDX11(2, 20, GF(5))),
    "The modular curve X_1(2, 20) over Finite field of size 5"
  );
end procedure;

test_print();