#############################################################################
##
#W logexp.gi              GUARANA package                     Bjoern Assmann
##
## functions for going back and forwards between groups and 
## lie algebras. 
##
#H  @(#)$Id$
##
#Y 2006
##
##

#############################################################################
##
#F GUARANA.AbstractLog( args )
##
## IN
## recLieAlg ..................... Lie Algebra record.
## g    .......................... Group element, either given as 
##                                 exponent vector or as an element of 
##                                 pcp group.
## rep_info ...................... indicates in which way the group element
##                                 is given and in which way the lie algbra 
##                                 element should be given. 
##                                 There are four possibilities:
##                                 "vecByVec"
##                                 "vecByElm"
##                                 "elmByVec"
##                                 "elmByElm"
## log_method .................... optional parameter. Decides with 
##                                 which method the logarithm 
##                                 should be computed. 
##                                 There are so far 2 possibilities
##                                 "simple"
##                                 "symbolic"
##              
## OUT
## Log( g ) either given as coefficient vector or as an elment of 
## a Lie algebra. 
##
GUARANA.AbstractLog := function( args )
    local recLieAlg, g, rep_info, log_method;

    recLieAlg := args[1];
    g := args[2];
    rep_info := args[3];
    if IsBound( args[4] ) then 
	log_method := args[4];
    else
	log_method := recLieAlg.log_method;
    fi;

    if log_method = "simple" then 
	if rep_info = "elmByElm" then 
	    return GUARANA.AbstractLog_Simple_ByElm( recLieAlg, g );
	elif rep_info = "elmByVec" then 
	    return GUARANA.AbstractLog_Simple_ByExponent( recLieAlg, g );
	elif rep_info = "vecByVec" then
	    return GUARANA.AbstractLogCoeff_Simple_ByExponent( recLieAlg, g);
	elif rep_info = "vecByElm" then
	    return GUARANA.AbstractLogCoeff_Simple_ByElm( recLieAlg, g );
	else 
	    Error( "wrong rep_info" );
	fi;
    elif log_method = "symbolic" then
	if rep_info = "vecByVec" then 
	    return GUARANA.Logarithm_Symbolic( [recLieAlg,g] );
	else 
	    Error( "Sorry method not implemented yet" );
	fi;
    else
	Error( "Wrong log_method" );
    fi;
end;
	
#############################################################################
##
#E
