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
