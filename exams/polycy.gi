#############################################################################
##
#W polycy.gi              GUARANA package                     Bjoern Assmann
##
## Examples of infinite polycyclic groups that can be used
## to test the Mal'cev collection.
##
#H  @(#)$Id$
##
#Y 2006
##
##

GUARANA.Tr_n_O := function( n, pol )
    local R, H, pcpH, noD, noU, ind_f, counter, ind_diag, ind_c, ind_n, 
          ind, new_pcs, G, indeces, N, NN, C, CC, i;

    R := GUARANA.Triang_PresentTriang( n, pol );
    H := R.Tr;
    pcpH := Pcp( H );
    # number of generators on diagoal
    noD := R.noD;
    # number of generators of unit group
    noU := noD/n;

    # get indeces of elms of finite order 
    ind_f := [];
    counter := 1;
    for i in [1..n] do 
	Add( ind_f, counter );
	counter := counter + noU;
    od;

    # get indeces of elms on diagonal of infinite order 
    ind_diag := [1..noD];
    ind_c := Difference( ind_diag, ind_f );

    # get indeces of elms of upper triangular part
    ind_n := [noD+1..Length(pcpH)];

    ind := Concatenation( ind_f, ind_c, ind_n );
    new_pcs := pcpH{ind};

    # change collector
    G := GUARANA.PcpGroupByPcs( H, new_pcs ); 

    indeces := [ [1..n], [n+1..noD], [noD+1..Length(pcpH)]];
    N := Subgroup( G, Pcp(G){indeces[3]} );
    NN := GUARANA.PcpGroupByPcs( G, Pcp(G){indeces[3]} );
    C := Subgroup( G, Pcp(G){indeces[2]} );
    CC := GUARANA.PcpGroupByPcs( G, Pcp(G){indeces[2]} );
    return [G,indeces,N,NN,C,CC];
end;

GUARANA.Tr_n_O1 := function( n )
    local x, pol;
    x := Indeterminate( Rationals );
    pol := x^2 - 3;
    return GUARANA.Tr_n_O( n, pol );
end;

GUARANA.MalcevColl_Tr_n_O1 := function( n )
    local l;
    l := GUARANA.Tr_n_O1( n );
    return MalcevCollectorConstruction( l );
end;

GUARANA.Tr_n_O2 := function( n )
    local x, pol;
    x := Indeterminate( Rationals );
    pol := x^3 - x^2 + 4;
    return GUARANA.Tr_n_O( n, pol );
end;

GUARANA.MalcevColl_Tr_n_O2 := function( n )
    local l;
    l := GUARANA.Tr_n_O2( n );
    return MalcevCollectorConstruction( l );
end;

GUARANA.Tr_n_O3 := function( n )
    local x, pol;
    x := Indeterminate( Rationals );
    pol := x^4+5*x^3-x^2+x+3;
    return GUARANA.Tr_n_O( n, pol );
end;

## functions for free-nilpotent-by automorphism examples.
##
GUARANA.F_nc_Aut1 := function( n, c )
    local F, F_nc, hl, auts, a, G, indeces, N, NN, C, CC;

    F := FreeGroup( n );
    F_nc := GUARANA.NilpotentQuotient( F, c );
    hl := HirschLength( F_nc );
    auts := GUARANA.GetSomeAutomorphsimOfF_nc( F_nc, n );
    a := auts[1]*auts[2]*auts[3]^3;
    G := GUARANA.PcpGroupTGroupByAbelianGroup( F_nc, [a] );
    indeces := [ [], [1], [2..hl+1] ];
    N := Subgroup( G, Pcp(G){indeces[3]} );
    NN := GUARANA.PcpGroupByPcs( G, Pcp(G){indeces[3]} );
    C := Subgroup( G, Pcp(G){indeces[2]} );
    CC := GUARANA.PcpGroupByPcs( G, Pcp(G){indeces[2]} );
    return [G,indeces,N,NN,C,CC];
end;

GUARANA.F_2c_Aut1 := function( c )
    return GUARANA.F_nc_Aut1( 2, c );
end;

GUARANA.MalcevColl_F_nc_Aut1 := function( n, c )
    local l;
    l := GUARANA.F_nc_Aut1( n, c );
    return MalcevCollectorConstruction( l );
end;

GUARANA.F_nc_Aut2 := function( n, c )
    local F, F_nc, hl, auts, a, G, indeces, N, NN, C, CC;

    F := FreeGroup( n );
    F_nc := GUARANA.NilpotentQuotient( F, c );
    hl := HirschLength( F_nc );
    auts := GUARANA.GetSomeAutomorphsimOfF_nc( F_nc, n );
    a := auts[1]*auts[3]^3;
    G := GUARANA.PcpGroupTGroupByAbelianGroup( F_nc, [a] );
    indeces := [ [], [1], [2..hl+1] ];
    N := Subgroup( G, Pcp(G){indeces[3]} );
    NN := GUARANA.PcpGroupByPcs( G, Pcp(G){indeces[3]} );
    C := Subgroup( G, Pcp(G){indeces[2]} );
    CC := GUARANA.PcpGroupByPcs( G, Pcp(G){indeces[2]} );
    return [G,indeces,N,NN,C,CC];
end;

GUARANA.F_3c_Aut2 := function( c )
    return GUARANA.F_nc_Aut2( 3, c );
end;

GUARANA.MalcevColl_F_nc_Aut2 := function( n, c )
    local l;
    l := GUARANA.F_nc_Aut2( n, c );
    return MalcevCollectorConstruction( l );
end;

## 
## OUT
## Examples that can be used to test Malcev collection. 
##
GUARANA.SomePolyMalcevExams := function( n )
    local T, H, ll, pcs, G, indeces, N, NN, C, CC, h1, h2, h3, h4, 
          f1, f2, f3, g1, g2, g3, g4, i, pcpH;

    # first 8 examples as in St Andrews paper 
    if n = 1 then 
	return GUARANA.Tr_n_O1( 3 );
    elif n = 2 then 
	return GUARANA.Tr_n_O1( 4 );
    elif n = 3 then 
	return GUARANA.Tr_n_O1( 5 );
    elif n = 4 then 
	return GUARANA.Tr_n_O2( 3 );
    elif n = 5 then 
	return GUARANA.Tr_n_O2( 4 );
    elif n = 6 then 
        return GUARANA.F_nc_Aut1( 2, 4 ); 
    elif n = 7 then 
        return GUARANA.F_nc_Aut1( 2, 5 );
    elif n = 8 then 
        return GUARANA.F_nc_Aut2( 3, 4 );
    fi;
    if n = 9 then 
	# in this example C and N have non trivial intersection
	T := GUARANA.NilpotentByFreeAbelianExams(2);
	H := T.G;
	pcpH := Pcp( H );
	h1 := pcpH[4]*pcpH[19];
	h2 := pcpH[6]*pcpH[19];
	h3 := pcpH[19];
	h4 := pcpH[20];
	f1 := pcpH[1]*pcpH[9];
	f2 := pcpH[3]*pcpH[12];
	f3 := pcpH[5]*pcpH[14];
        pcs := [f1,f2,f3,h1,h2];
	for i in [1..Length( Pcp(T.N)) ] do 
	    Add( pcs, Pcp(T.N)[i] );
	od;
	G := GUARANA.PcpGroupByPcs( H, pcs );
	indeces := [[1..3],[4,5],[6..17]];
	N := Subgroup( G, Pcp(G){indeces[3]} );
	NN := GUARANA.PcpGroupByPcs( G, Pcp(G){indeces[3]} );
	g1 := Pcp(G)[4];
	g2 := Pcp(G)[5];
	g3 := Pcp(G)[16];
	g4 := Pcp(G)[17];
	C := Subgroup( G, [g1,g2,g3,g4] );
	CC := GUARANA.PcpGroupByPcs( G, [g1,g2,g3,g4] );
	return [G,indeces,N,NN,C,CC];
    fi;
end;

ExamplesOfSomeMalcevCollectors := function( n )
    local l;
    l := GUARANA.SomePolyMalcevExams( n );
    return MalcevCollectorConstruction( l );
end;

#############################################################################
##
#E 
