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
GUARANA.NilpotentQuotient := function( ) return 0; end;
GUARANA.NilpotentEngelQuotient := function( ) return 0; end;

if IsPackageMarkedForLoading( "nq" , "2.0" ) <> true then 
    GUARANA.NqAvailable := false;
else
    GUARANA.NqAvailable := true;
    GUARANA.NilpotentQuotient := NilpotentQuotient;
    GUARANA.NilpotentEngelQuotient := NilpotentEngelQuotient;
fi;


#############################################################################
##
#F GUARANA.Examples_FreeNilpotentGrp( n, c )
##
##
## IN
## n................... number of generators
## c .................  nilpotency class 
## 
## OUT
## Pcp of a free nilpotent group on n generators of class c
##
GUARANA.Examples_FreeNilpotentGrp := function( n, c )
    local F,N;
    F := FreeGroup( n );
    LoadPackage( "nq" );
    N := GUARANA.NilpotentQuotient( F, c );
    return N;
end;



#############################################################################
##
#F GUARANA.Examples_Unitriangular( dim , degree )
##
## IN
## degree ........... degree of an example number field K over Q.
##                    currently this can be only 2 or 3.
## dim   ............ degree/dim  of matrix group 
##
## OUT 
## Pcp of a upper unitriangular matrix Tr_1(n,O) where O is 
## the maximal order of K.
##
GUARANA.Examples_Unitriangular := function( dim, degree  )
    local x,pol,R;
    x := Indeterminate( Rationals );
    if degree = 2 then
        pol := x^2-3;
    elif degree = 3 then
        pol :=  x^3 - x^2 + 4;
    else
        Error( "Sorry no appropriate polynomial\n" );
    fi;
    R := GUARANA.Triang_PresentTriang( dim, pol );
    return GUARANA.Triang_UpperTriangAndUnitriang( R ).N;
end;

#############################################################################
##
#F GUARANA.Examples_Engel( n, c )
##
## IN
## n................... number of generators
## c .................  nilpotency class 
## 
## OUT
## Pcp of  Engel groups of Werners paper. 
## Use first NilpotentEngelQuotient and then factor
## torsion out.
##
GUARANA.Examples_Engel := function( n, c )
    local G,T,H,N;
    G := GUARANA.NilpotentEngelQuotient( FreeGroup(n), c );
    T := TorsionSubgroup( G );
    H := G/T;
    N := PcpGroupBySeries( UpperCentralSeries(H), "snf" );
    return N;
end;

#############################################################################
##
## Examples coming from the package "polycyclic"
##
GUARANA.ExamplesOfSomeTGroups := function()
    local G, FF, F,n;
    n := 10;
    G := List( [1..n], x-> ExamplesOfSomePcpGroups(x) );
    FF := List( G, FittingSubgroup );
    F := List( [1..n], x-> PcpGroupByPcp( Pcp( FF[x] ) ) );
    List( F, IsNilpotent );
    return rec( G := G, FF := FF, F := F );

end;


#############################################################################
# additional  ideas to produce T-groups
#
# - take subgroups
# - Nilpotent quotient of other finitely presented groups torsion-free groups
# - Nilpotent quotient of other finitely presented groups and then
#   factor torsion out. 

#############################################################################
##
#E
