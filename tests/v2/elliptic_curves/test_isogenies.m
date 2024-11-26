procedure test_MDIsogenies()

	p:=5;
	E:=EllipticCurve([GF(11)!1,3]);
	j:=jInvariant(E);
	isogenies:=MDIsogenies(E,p);
	TSTAssertEQ(#isogenies,p+1);
	for i in [1..#isogenies] do
		TSTAssertEQ(jInvariant(Domain(isogenies[i])), j);
	end for;
	jset:={jInvariant(Codomain(isogenies[i])) : i in [1..#isogenies]};
	TSTAssertEQ(#jset,p+1); //this might fail some times ir one changes E or p - this is expected, this means there is "some CM".
	
end procedure;

test_MDIsogenies();