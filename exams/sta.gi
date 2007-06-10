#############################################################################
##
#W sta.gi                 GUARANA package                     Bjoern Assmann
##
## Examples that were used in the Groups St Andrews 2005 paper 
## "Algorithmic use of the Mal'cev correspondence". 
##
#H  @(#)$Id$
##
#Y 2006
##
##

#############################################################################
##
#F GUARANA.GetSomeAutomorphsimOfF_nc( F_nc, n )
##
## IN
## F_nc ..... free nilpotent group on n gens and class c
## n ........ number of generators.
##
## OUT
## A list of  autmorphsims of F_nc.  
## Some of them are induced by automorphisms of F_n.
## 
GUARANA.GetSomeAutomorphsimOfF_nc := function( F_nc, n ) 
    local F,A,gensA,gensN,k,auts,ll,gensF,imgs,g,beta,i,aut;

    # get automorphsim induced from free group
    F := FreeGroup( n );
    gensF := GeneratorsOfGroup( F );
    A := AutomorphismGroup( F );
    gensA := GeneratorsOfGroup( A );
    gensN := Pcp( F_nc );
    k := Length( gensA );
    auts := [];
    for i in [1..k] do
        aut := gensA[i];
        ll := List( gensF, x->x^aut );
        imgs := List( ll, x -> MappedWord( x, gensF, gensN{[1..n]} ) );
        Add( auts, GroupHomomorphismByImagesNC(F_nc,F_nc,gensN{[1..n]},imgs ));
    od;

    # get  automorphism of Bryant's paper 
    # "Automorphism groups of free nilpotent groups" 
    # Arch. Math. Basel. 52, 1989.
    g := gensN[1]*Comm( Comm( gensN[1], gensN[2] ), gensN[1] );
    beta := GroupHomomorphismByImagesNC( F_nc,F_nc,gensN{[1,2]}, [g,gensN[2]]);
    Add( auts, beta );

    return auts;
end;


#############################################################################
##
#F GUARANA.PcpTGroupByAbelianGroup( N, auts )
##
## IN
## N ..... T-group given by a finite presentation
## auts... List Automorphisms of N, where <auts> is abelian
## 
## OUT
## polycyclic presentation of the semidirect product H \rhd N
## where H is a free abelian group with Length(auts) generators, which
## acts like auts on N
##
GUARANA.PcpGroupTGroupByAbelianGroup := function( N, auts )
    local n,m,l,coll,pcpN,exp1,i,j,con,exp,exp2,genList,g;
    
    #get collector
    n := HirschLength( N );
    m := Length( auts );
    l := m+n;
    coll := FromTheLeftCollector( l );
    
    # define conjugation relations within N,   g_i^g_j
    pcpN := Pcp( N );
    exp1 := List( [1..m], x-> 0 );
    for i in [2..n] do 
        for j in [1..(i-1)] do
            # compute g_i^g_j
            con := pcpN[i]^pcpN[j];
            exp2 := Exponents( con ); 
            exp := Concatenation( exp1, exp2 );
            genList:=GUARANA.Exp2GenList(exp);
            SetConjugate( coll, i+m, j+m, genList);

            # compute g_i^(g_j^-1)
            con := pcpN[i]^(pcpN[j]^-1);
            exp2 := Exponents( con ); 
            exp := Concatenation( exp1, exp2 );
            genList:=GUARANA.Exp2GenList(exp);
            SetConjugate( coll, i+m, -(j+m), genList);
        od;
    od;

    #define relations coming from the action of H on N
    for j in [1..m] do
        for i in [1..n] do
            #action of g_j on g_{i+m}
            g := pcpN[i]^(auts[j]);
            exp2 := Exponents( g );
            exp := Concatenation( exp1, exp2 );
            genList:=GUARANA.Exp2GenList(exp);
            SetConjugate( coll, i+m, j, genList);

            #action of g_j^-1 onf g_{i+m}
            g := pcpN[i]^(auts[j]^-1);
            exp2 := Exponents( g );
            exp := Concatenation( exp1, exp2 );
            genList:=GUARANA.Exp2GenList(exp);
            SetConjugate( coll, i+m, -j, genList);
        od;
    od;

    UpdatePolycyclicCollector( coll);
    return PcpGroupByCollector( coll );
end;


#############################################################################
##
#F GUARANA.NilpotentByFreeAbelianExams( n )
##
## IN
## n ........... natural number between 1 and 8 for the 
##               Sta paper examples.
##               Further choices [-11..-1]
##               [10..29]. 
##
## OUT 
## Pcp of a nilpotent by free abelian group.
##
GUARANA.NilpotentByFreeAbelianExams := function( n )
    local x,pol,R,G,m,ll,N,n1,n2,N2,F,F_nc,hl,auts,T,a;

    x := Indeterminate( Rationals );

    # first 8 groups as in Groups St Andrews 2005 paper
    if n=1 then 
        pol := x^2-3;
        R := GUARANA.Triang_PresentTriang( 3, pol );
        return GUARANA.Triang_UpperTriangAndUnitriang( R );
    fi;
    if n=2 then 
        pol := x^2-3;
        R := GUARANA.Triang_PresentTriang( 4, pol );
        return GUARANA.Triang_UpperTriangAndUnitriang( R );
    fi;
    if n=3 then
        pol := x^2-3;
        R := GUARANA.Triang_PresentTriang( 5, pol );
        return GUARANA.Triang_UpperTriangAndUnitriang( R );
    fi;
    if n=4 then
        pol :=  x^3 - x^2 + 4; 
        R := GUARANA.Triang_PresentTriang( 3, pol );
        return GUARANA.Triang_UpperTriangAndUnitriang( R );
    fi;
    if n=5 then
        pol :=  x^3 - x^2 + 4; 
        R := GUARANA.Triang_PresentTriang( 4, pol );
        return GUARANA.Triang_UpperTriangAndUnitriang( R );
    fi;
    if n=6 then
        m := 2;
        F := FreeGroup( m );
        F_nc := GUARANA.NilpotentQuotient( F, 4 );
        hl := HirschLength( F_nc );
        auts := GUARANA.GetSomeAutomorphsimOfF_nc( F_nc, m );
        a := auts[1]*auts[2]*auts[3]^3;
        G := GUARANA.PcpGroupTGroupByAbelianGroup( F_nc, [a] );
        N := Subgroup( G, Pcp(G){[2..hl+1]} );
        T := rec( G := G, N := N );
        return T;
    fi;
    if n=7 then
        m := 2;
        F := FreeGroup( m );
        F_nc := GUARANA.NilpotentQuotient( F, 5 );
        hl := HirschLength( F_nc );
        auts := GUARANA.GetSomeAutomorphsimOfF_nc( F_nc, m );
        a := auts[1]*auts[2]*auts[3]^3;
        G := GUARANA.PcpGroupTGroupByAbelianGroup( F_nc, [a] );
        N := Subgroup( G, Pcp(G){[2..hl+1]} );
        T := rec( G := G, N := N );
        return T;
    fi;
    if n=8 then
        m := 3;
        F := FreeGroup( m );
        F_nc := GUARANA.NilpotentQuotient( F, 4 );
        hl := HirschLength( F_nc );
        auts := GUARANA.GetSomeAutomorphsimOfF_nc( F_nc, m );
        a := auts[1]*auts[3]^3;
        G := GUARANA.PcpGroupTGroupByAbelianGroup( F_nc, [a] );
        N := Subgroup( G, Pcp(G){[2..hl+1]} );
        T := rec( G := G, N := N );
        return T;
    fi;
    # other examples 
    if n=11 then
        pol :=  x^3 - x^2 + 4; 
        R := GUARANA.Triang_PresentTriang( 4, pol );
        return GUARANA.Triang_UpperTriangAndUnitriang( R );
    fi;
    if n=12 then 
        pol :=  x^3 - x^2 + 4; 
        R := GUARANA.Triang_PresentTriang( 3, pol );
        G := R.Tr;
        m := Length( Pcp(G) );
        ll := [R.noD+1..m];
        N := Subgroup( G, Pcp(G){ll} );     
        n1 := Pcp(N)[2]*Pcp(N)[4];
        n2 := Pcp(N)[5]*Pcp(N)[1]^-1;
        N2 := NormalClosure( G, Subgroup( G, [n1,n2]) );
        return rec( G := G, N := N2 );
    fi;
    if n=13 then
        pol := x^7-x^6+x^5-x^4-3*x^3+3*x^2-2*x+1;
        R := GUARANA.Triang_PresentTriang( 3, pol );
        return GUARANA.Triang_UpperTriangAndUnitriang( R );
    fi;
    if n=14 then
        pol :=  x^3 - x^2 + 4;
        R := GUARANA.Triang_PresentTriang( 4, pol );
        return GUARANA.Triang_UpperTriangAndUnitriang( R );
    fi;
    if n=15 then 
        pol :=  x^3 - x^2 + 4; 
        R := GUARANA.Triang_PresentTriang( 4, pol );
        G := R.Tr;
        m := Length( Pcp(G) );
        ll := [R.noD+1..m];
        N := Subgroup( G, Pcp(G){ll} );     
        n1 := Pcp(N)[2]*Pcp(N)[5]*Pcp(N)[13]^-1;
        n2 := Pcp(N)[3]*Pcp(N)[7]*Pcp(N)[11]^-1;
        N2 := NormalClosure( G, Subgroup( G, [n1,n2]) );
        return rec( G := G, N := N2 );
    fi;
    if n = 11 then
        pol := x^2 -3;
        R := GUARANA.Triang_PresentTriang( 5, pol );
        return GUARANA.Triang_UpperTriangAndUnitriang( R );
    fi;
    if n = 21 then 
        m := 2;
        F := FreeGroup( m );
        F_nc := GUARANA.NilpotentQuotient( F, 5 );
        hl := HirschLength( F_nc );
        auts := GUARANA.GetSomeAutomorphsimOfF_nc( F_nc, 2 );
        a := auts[1]*auts[2]*auts[3]^3;
        G := GUARANA.PcpGroupTGroupByAbelianGroup( F_nc, [a] );
        N := Subgroup( G, Pcp(G){[2..hl+1]} );
        T := rec( G := G, N := N );
        return T;
    fi;
    if n = 29 then 
        # the calculation of the presentation needs 3000 sec.
        # Dimension = 34
        F := FreeGroup( 2 );
        return GUARANA.NilpotentQuotient( F, 6 );
    fi;
    if n in [-10..-1] then 
        G := ExamplesOfSomePcpGroups(-n);
        N := FittingSubgroup( G );
        return rec( G := G, N := N );
    fi;
    if n = -11 then
       G := ExamplesOfSomePcpGroups(13);
       N := FittingSubgroup( G );
       return rec( G := G, N := N );
    fi; 
end;

#############################################################################
##
#E 
