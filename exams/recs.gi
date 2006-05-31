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
    ll := [];
    for i in [1..c] do
        Print( "Free nilpotent group ", n, " ", i, "\n" );
        N := GUARANA.Examples_FreeNilpotentGrp( n, i );
        r := GUARANA.TGroupRec( [N] );
        Add( ll, r );
    od;
    return ll;
end;

GUARANA.Get_FNG_TGroupRecord := function( n, c )
    local N, r;
    N := GUARANA.Examples_FreeNilpotentGrp( n, c );
    r := GUARANA.TGroupRec( [N] );
    return r; 
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
    ll := [];
    for i in [2..dim] do
        Print( "Untriangular group ", degree, " ", i, "\n" );
        if degree = 2 then 
            N := GUARANA.Tr_n_O1( i )[4];
        elif degree =3 then 
            N := GUARANA.Tr_n_O2( i )[4];
        else 
            Error( "Wrong degree " );
        fi;
	    r:= GUARANA.TGroupRec( [N] );
        Add( ll, r );
    od;
    return ll;
end;

GUARANA.Get_Unitriangular_TGroupRecord := function( dim , degree )
    local N, r;
    if degree = 2 then 
        N := GUARANA.Tr_n_O1( dim )[4];
    elif degree =3 then 
        N := GUARANA.Tr_n_O2( dim )[4];
    else 
        Error( "Wrong degree " );
    fi;
    r:= GUARANA.TGroupRec( [N] );
    return r;
end;

#############################################################################
##
#F GUARANA.Get_FNG_LieAlgRecords( n, c )
##
## IN
## n ................ number of generators of free group F that will be used
## c ................ upper bound for the nilpotency class of quotients
##                    F that will be used. 
## 
GUARANA.Get_FNG_LieAlgRecords := function( n, c )
    local ll,ll2; 
    ll := GUARANA.Get_FNG_TGroupRecords( n, c );
    ll2 := List( ll, x-> GUARANA.LieAlgebraByTGroupRec( [x] ) );
    return ll2;
end;

GUARANA.Get_FNG_MalcevObjects := function( n, c )
    local ll,ll2; 
    ll := GUARANA.Get_FNG_TGroupRecords( n, c );
    ll2 := List( ll, x-> MalcevObjectConstruction( x ) ); 
    return ll2;
end;

GUARANA.Get_FNG_MalcevObject := function( n, c )
    local rT, mo;
    rT := GUARANA.Get_FNG_TGroupRecord( n, c );
    mo := MalcevObjectConstruction( rT); 
    return mo;
end;

#############################################################################
##
#F GUARANA.Get_Unitriangular_LieAlgRecords( dim , degree )
##
## IN
## dim ..... upper bound for dimension of matrix groups that are going to 
##           used.
## degree .. degreee of polynomial that is going to be used. 
## 
GUARANA.Get_Unitriangular_LieAlgRecords := function( dim , degree )
    local ll,ll2;
    ll := GUARANA.Get_Unitriangular_TGroupRecords( dim, degree );
    ll2 := List( ll, x-> GUARANA.LieAlgebraByTGroupRec( [x] ) );
    return ll2;
end;

GUARANA.Get_Unitriangular_MalcevObjects := function( dim , degree )
    local ll,ll2;
    ll := GUARANA.Get_Unitriangular_TGroupRecords( dim, degree );
    ll2 := List( ll, x-> MalcevObjectConstruction( x ) );
    return ll2;
end;

GUARANA.Get_Unitriangular_MalcevObject := function( dim , degree )
    local rT, mo;
    rT :=GUARANA.Get_Unitriangular_TGroupRecord( dim, degree );
    mo := MalcevObjectConstruction( rT );
    return mo;
end;

#############################################################################
##
#E 
