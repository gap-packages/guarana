# Code for testing the BCH Algorithms                         22.10.05
#
#

# 1. Free nilpotent groups on n generators and class c
BCH_Examples_FreeNilpotentGrp := function( n, c )
    local F,N;
    F := FreeGroup( n );
    LoadPackage( "nq" );
    N := NilpotentQuotient( F, c );
    return N;
end;

# 2. Tr_1(n,O) where O is the maximal order of some number field.
BCH_Examples_Unitriangular := function( dim, degree  )
    local x,pol,R;
    x := Indeterminate( Rationals );
    if degree = 2 then
        pol := x^2-3;
    elif degree = 3 then
        pol :=  x^3 - x^2 + 4;
    else
        Error( "Sorry no appropriate polynomial\n" );
    fi;

    R := PresentTriang( dim, pol );
    return SC_Exams_Help1( R ).N;
end;

# Engel groups of Werners paper, NilpotentEngelQuotient and then factor
# torsion out.
BCH_Examples_Engel := function( n, c )
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

BCH_Get_FNG_TGroupRecords := function( n, c )
    local i,ll,N,r;
    ll := [[n,c]];
    for i in [1..c] do
        Print( "Free nilpotent group ", n, " ", i, "\n" );
        N := BCH_Examples_FreeNilpotentGrp( n, i );
        r := BCH_TGroupRec( N );
        Add( ll, r );
    od;
    return ll;
end;

# dim ..... upper bound for dimension
BCH_Get_Unitriangular_TGroupRecords := function( dim , degree )
    local i,ll,N,r;
    ll := [[dim,degree]];
    for i in [2..dim] do
        Print( "Untriangular group ", degree, " ", i, "\n" );
        N := BCH_Examples_Unitriangular( i, degree );
        r := BCH_TGroupRec( N );
        Add( ll, r );
    od;
    return ll;
end;

