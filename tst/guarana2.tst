gap> START_TEST("Test 2 of guarana package");
gap> ll := GUARANA.Tr_n_O1( 3 );
[ Pcp-group with orders [ 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], 
  [ [ 1 .. 3 ], [ 4 .. 6 ], [ 7 .. 12 ] ], 
  Pcp-group with orders [ 0, 0, 0, 0, 0, 0 ], 
  Pcp-group with orders [ 0, 0, 0, 0, 0, 0 ], 
  Pcp-group with orders [ 0, 0, 0 ], Pcp-group with orders [ 0, 0, 0 ] ]
gap> malCol := MalcevCollectorConstruction( ll );
<<Malcev collector>>
  F : [ 2, 2, 2 ]
  C : <<Malcev object of dimension 3>>
  N : <<Malcev object of dimension 6>>
gap> exps_g := [ 1, 1, 1, -3, -2, 1, -2, -1, 0, 3, -1,3 ];
[ 1, 1, 1, -3, -2, 1, -2, -1, 0, 3, -1, 3 ]
gap> exps_h := [ 1, 0, 1, -1, 0, 2, 0, 4, -1, 5, 9,-5 ];
[ 1, 0, 1, -1, 0, 2, 0, 4, -1, 5, 9, -5 ]
gap> g := MalcevGElementByExponents( malCol, exps_g );
[ 1, 1, 1, -3, -2, 1, -2, -1, 0, 3, -1, 3 ]
gap> h := MalcevGElementByExponents( malCol, exps_h );
[ 1, 0, 1, -1, 0, 2, 0, 4, -1, 5, 9, -5 ]
gap> k := g*h;
[ 0, 1, 0, -4, -2, 3, -7, 0, -37, -16, -352, -212 ]
gap> STOP_TEST( "guarana2.tst", 20000);

