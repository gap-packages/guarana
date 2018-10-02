#############################################################################
#W help.gi                 GUARANA package                     Bjoern Assmann
##
## Some help functions 
##
#H  @(#)$Id$
##
#Y 2006
##

############################################################################
##
GUARANA.Exp2GenList := function( exp )
    local  n, genList, i;
    n := Length( exp );
    genList := [  ];
    for i  in [ 1 .. n ]  do
	if exp[i] <> 0  then
	    Append( genList, [ i, exp[i] ] );
	fi;
    od;
    return genList;
end;

#############################################################################
##
GUARANA.CopyVectorList := function( list )
    local  i, j, k, list2;
    list2 := [  ];
    for i  in [ 1 .. Length( list ) ]  do
	Add( list2, [  ] );
	for j  in [ 1 .. Length( list[i] ) ]  do
	    Add( list2[i], [  ] );
	    list2[i][j] := list[i][j];
	od;
    od;
    return list2;
end;
#############################################################################
##
GUARANA.Random_IntegralExpVector := function( recLieAlg, range )
    local n,ll,vec;
    n := recLieAlg.dim;
    ll := [ - range .. range ];
    vec := List( [ 1 .. n ], function ( x )
            return RandomList( ll );
        end );
    return vec;
end;

#############################################################################
##
#F GUARANA.GrpElmByExpsAndPcs( pcs, exp )
## 
## IN
## pcs ............................ polycyclic sequence
## exp ............................ exponent vector 
##
## OUT
## pcs^exp
##
## Note that exp is allowed to be shorter then pcs 
##
GUARANA.GrpElmByExpsAndPcs := function( pcs, exp )
    local elm, i;
    
    elm := pcs[1]^exp[1];
    for i in [2..Length(exp)] do
	elm := elm * pcs[i]^exp[i];
    od;
    return elm;
end;

GUARANA.GrpElmByExpsAndCollNC := function( coll, exp )
    local word;
     word := ObjByExponents( coll, exp );
     return PcpElementByGenExpListNC( coll, word );
end;


GUARANA.CompleteRuntime1:= function( func, input )
     local rec1,rec2, user_time, user_time_child, system_time,
     system_time_child, sum;
     rec1 := Runtimes();
     func( input );
     rec2 := Runtimes();

     user_time := rec2.user_time -rec1.user_time;
     user_time_child := rec2.user_time_children -rec1.user_time_children;
     system_time := rec2.system_time - rec1.system_time;
     system_time_child := rec2.system_time_children - rec1.system_time_children;

     sum := user_time + user_time_child + system_time + system_time_child;

     return sum;
end;

GUARANA.CompleteRuntime2:= function( func, input )
    local rec1,rec2, user_time, user_time_child, system_time,
    system_time_child, sum, result;
    rec1 := Runtimes();
    result := func( input );
    rec2 := Runtimes();

    user_time := rec2.user_time -rec1.user_time;
    user_time_child := rec2.user_time_children -rec1.user_time_children;
    system_time := rec2.system_time - rec1.system_time;
    system_time_child := rec2.system_time_children - rec1.system_time_children;

    sum := user_time + user_time_child + system_time + system_time_child;

    return rec( time := sum, result := result );
end;

#########################################################################
##
## This function is only temporarily necessary untill 
## I can use IsWeightedCollector from "polycyclic" again.
## In gap.dev USE_COMBINATORIAL_COLLECTOR is false and thus 
## IsWeigthedCollector does not work.
## 
GUARANA.IsWeightedCollector :=  function( coll )

    # check wehter weights are already computed and compute them 
    # if necessary
    if not IsBound( coll![PC_WEIGHTS] ) then
        if FromTheLeftCollector_SetWeights( coll ) = fail then 
            Error( "Computation of weights not possible\n" );
            return false;
        fi;
    fi;
    SetFilterObj( coll, IsWeightedCollector );
    return true;
end;

