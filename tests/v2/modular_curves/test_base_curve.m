procedure test_print()
  TSTAssertEQ(
    Sprint(New(MDCrvMod)),
    "The MDCrvMod class shouldn't be used directly"
  );
end procedure;

test_print();