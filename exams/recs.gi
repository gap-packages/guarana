#############################################################################
##
#W recs.gi                GUARANA package                     Bjoern Assmann
##
## Some Examples of Mal'cev records.  
##
#H  @(#)$Id$
##
#Y 2006
##
##

#############################################################################
##
#F GUARANA.Get_FNG_TGroupRecords( n, c )
##
## IN
## n ................ number of generators of free group F that will be used
## c ................ upper bound for the nilpotency class of quotients
##                    F that will be used. 
## OUT
## A Liste containing all free nilpotent groups of class [1..c]
## 
GUARANA.Get_FNG_TGroupRecords := function( n, c )
    local i,ll,N,r;
    ll := [[n,c]];
    for i in [1..c] do
        Print( "Free nilpotent group ", n, " ", i, "\n" );
        N := GUARANA.Examples_FreeNilpotentGrp( n, i );
        r := GUARANA.TGroupRec( N );
        Add( ll, r );
    od;
    return ll;
end;

#############################################################################
##
#F GUARANA.Get_Unitriangular_TGroupRecords( dim , degree )
##
## IN
## dim ..... upper bound for dimension of matrix groups that are going to 
##           used.
## degree .. degreee of polynomial that is going to be used. 
## 
GUARANA.Get_Unitriangular_TGroupRecords := function( dim , degree )
    local i,ll,N,r;
    ll := [[dim,degree]];
    for i in [2..dim] do
        Print( "Untriangular group ", degree, " ", i, "\n" );
        N := GUARANA.Examples_Unitriangular( i, degree );
        r := GUARANA.TGroupRec( N );
        Add( ll, r );
    od;
    return ll;
end;

#############################################################################
##
#E 
