#############################################################################
##
#W tgrps.gi               GUARANA package                     Bjoern Assmann
##
## Examples of T-groups.
##
#H  @(#)$Id$
##
#Y 2006
##
##
##

#############################################################################
##
#F Guarana.Examples_FreeNilpotentGrp( n, c )
##
##
## IN
## n................... number of generators
## c .................  nilpotency class 
## 
## OUT
## Pcp of a free nilpotent group on n generators of class c
##
Guarana.Examples_FreeNilpotentGrp := function( n, c )
    local F,N;
    F := FreeGroup( n );
    LoadPackage( "nq" );
    N := NilpotentQuotient( F, c );
    return N;
end;



#############################################################################
##
#F Guarana.Examples_Unitriangular( dim , degree )
##
## IN
## dim .............. dim of an example number field K over Q.
##                    currently this can be only 2 or 3.
## degree ........... degree of matrix group 
##
## OUT 
## Pcp of a upper unitriangular matrix Tr_1(n,O) where O is 
## the maximal order of K.
Guarana.Examples_Unitriangular := function( dim, degree  )
    local x,pol,R;
    x := Indeterminate( Rationals );
    if degree = 2 then
        pol := x^2-3;
    elif degree = 3 then
        pol :=  x^3 - x^2 + 4;
    else
        Error( "Sorry no appropriate polynomial\n" );
    fi;
    R := GURANA.Triang_PresentTriang( dim, pol );
    return GURANA.Triang_UpperTriangAndUnitriang( R ).N;
end;

#####################

# Engel groups of Werners paper, NilpotentEngelQuotient and then factor
# torsion out.
Guarana.Examples_Engel := function( n, c )
    local G,T,H,N;
    G := NilpotentEngelQuotient( FreeGroup(n), c );
    T := TorsionSubgroup( G );
    H := G/T;
    N := PcpGroupBySeries( UpperCentralSeries(H), "snf" );
    return N;
end;


#############################################################################
# additional  ideas to produce T-groups
#
# - take subgroups
# - Nilpotent quotient of other finitely presented groups torsion-free groups
# - Nilpotent quotient of other finitely presented groups and then
#   factor torsion out. 

Guarana.Get_FNG_TGroupRecords := function( n, c )
    local i,ll,N,r;
    ll := [[n,c]];
    for i in [1..c] do
        Print( "Free nilpotent group ", n, " ", i, "\n" );
        N := Guarana.Examples_FreeNilpotentGrp( n, i );
        r := Guarana.TGroupRec( N );
        Add( ll, r );
    od;
    return ll;
end;

# dim ..... upper bound for dimension
Guarana.Get_Unitriangular_TGroupRecords := function( dim , degree )
    local i,ll,N,r;
    ll := [[dim,degree]];
    for i in [2..dim] do
        Print( "Untriangular group ", degree, " ", i, "\n" );
        N := Guarana.Examples_Unitriangular( i, degree );
        r := Guarana.TGroupRec( N );
        Add( ll, r );
    od;
    return ll;
end;

