procedure TSTAssertEQ(computed, expected)
  if computed eq expected then return; end if;
  assert false;
end procedure;

procedure TSTAssertNE(computed, expected)
  if computed ne expected then return; end if;
  assert false;
end procedure;