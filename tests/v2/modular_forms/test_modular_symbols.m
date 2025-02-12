procedure test_MDCongruenceSubgroup()
  for G in [Gamma0(112), Gamma1(17), CongruenceSubgroup(11)] do
    M := ModularSymbols(G,2);
    S := NewSubspace(CuspidalSubspace(M));
    G1 := MDCongruenceSubgroup(S);
    TSTAssertEQ(G, G1);
  end for;
end procedure;

expected_al := [
    [[ 1 ], [ -1 ], [ 1, -1 ]],
    [[ -1 ], [ 1 ], [ -1, -1 ]]
];

expected_q_basis :=
[
    [ 1, -1, 1, 1, 0, -1, -1, -1, -2, 0 ],
    [ 1, 1, -1, 1, -4, -1, 3, 1, -2, -4 ],
    [ 1, 2, -2, -2, 3, -4, -1, -4, 1, 6 ],
    [ 1, -2, -2, -2, 3, 4, -1, 4, 1, -6 ]
];

procedure test_MDDiagonalALBasis()
  G := Gamma0(38);
  M := ModularSymbols(G,2);
  S := CuspidalSubspace(M);
  dec := NewformDecomposition(S);
  q_basis, al := MDDiagonalALBasis(dec);
  TSTAssertEQ(al, expected_al);
  q_basis := [[Coefficient(f, j) : j in [1..10]] : f in &cat q_basis];
  TSTAssertEQ(q_basis, expected_q_basis);
end procedure;

test_MDCongruenceSubgroup();
test_MDDiagonalALBasis();