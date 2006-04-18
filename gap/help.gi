
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
