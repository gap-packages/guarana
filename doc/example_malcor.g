n := 2;
F := FreeGroup( n );
c := 3;
N := NilpotentQuotient( F, c );

mo := MalcevObjectByTGroup( N );

dim := Dimension( mo );
UnderlyingGroup( mo );
UnderlyingLieAlgebra( mo );

g := MalcevGrpElementByExponents( mo, [1,1,0,2,-1/2] );    
x := MalcevLieElementByCoefficients( mo, [1/2, 2, -1, 3, 5 ] );

h := RandomGrpElm( mo );
y := RandomLieElm( mo );

z := Log( g );
Exp( z ) = g;

k := Exp( y );
Log( k ) = y;

g*h;
Comm(g,h);
Comm(x,y);

indets := List( List( [1..dim], i->Concatenation( "a_", String(i) ) ),
                  x->Indeterminate( Rationals, x : new ) );
g_sym := MalcevSymbolicGrpElementByExponents( mo, indets );
x_sym := Log( g_sym );
g_sym * g;


# not in manual 

SetLogAndExpMethod( mo, "pols" );
SetStarMethod( mo, "pols" );

# check whether groups which are not given by Mal'cev basis are rejected.
K := NilpotentEngelQuotient( F, 3 );
mo2 := MalcevObjectByTGroup( K );

