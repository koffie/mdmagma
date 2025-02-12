intrinsic MDRestriction(A :: AlgMatElt, V :: ModTupRng) -> AlgMatElt
{
  Return the matrix A where the domain and codomain of A are restricted to
  V. This assumes that V*A is a subspace of V.
}
  vA := [v*A : v in Basis(V)];
  return Matrix([Coordinates(V,v) : v in vA]);
end intrinsic;