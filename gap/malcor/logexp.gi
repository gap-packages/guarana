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
## EXAMPLE
## GUARANA.AbstractLog( [malcevRec.recL_NN, g, "vecByVec"] ); 
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

    # catch trivial case
    if rep_info = "vecByVec" then 
	if Length( g ) = 0 then 
	    return [];
	fi;
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
#F GUARANA.AbstractExp( args )
##
## IN
## recLieAlg ......................... Lie algebra record
## x ................................. Lie algebra elment etiher 
##                                     via a coefficient vector or as an
##                                     element of a Lie algebra (in the 
##                                     GAP sense).
## rep_info .......................... indicate in which way the Lie 
##                                     algebra elment is given and 
##                                     in which way the group element should
##                                     be given. There are 4 possibilities:
##                                     "vecByVec"
##                                     "vecByElm"
##                                     "elmByVec"
##                                     "elmByElm"
## exp_method ........................ optional parameter. Decides with 
##                                     which method the exponential should
##                                     be computed. So far there 2 possib.
##                                     "simple"
##                                     "symbolic"
## 
## OUT
## Exp( x ) either given as epxonent vector or as an element of a Pcp group
##
GUARANA.AbstractExp := function( args )
    local recLieAlg, x, rep_info, exp_method;

    recLieAlg := args[1];
    x := args[2];
    rep_info := args[3];
    if IsBound( args[4] ) then 
	exp_method := args[4];
    else
	exp_method := recLieAlg.exp_method;
    fi;

    if exp_method = "simple" then 
	if rep_info = "elmByElm" then 
	    Error( "Sorry method not implemented yet" );
	elif rep_info = "elmByVec" then 
	    Error( "Sorry method not implemented yet" );
	elif rep_info = "vecByVec" then
	    return GUARANA.Abstract_Exponential_ByVector( recLieAlg, x);
	elif rep_info = "vecByElm" then
	    return GUARANA.Abstract_Exponential_ByElm( recLieAlg, x );
	else 
	    Error( "wrong rep_info" );
	fi;
    elif exp_method = "symbolic" then
	if rep_info = "vecByVec" then 
	    return GUARANA.Exponential_Symbolic( [recLieAlg,x] );
	else 
	    Error( "Sorry method not implemented yet" );
	fi;
    else
	Error( "Wrong exp_method" );
    fi;
end;

## IN
## recLieAlg .............................. Lie algebra record of a 
##                                          Lie algebra L
## x,y .................................... elements of L 
## rep_info .......................... indicate in which way the Lie 
##                                     algebra elments are given.
##                                     There are two  possibilities:
##                                     "vec"
##                                     "elm"
## star_method ........................ optional parameter. Decides with 
##                                     which method star should
##                                     be computed. So far there 2 possib.
##                                     "simple"
##                                     "symbolic"
## 
## OUT
## x*y either given in the form x,y are given.
##
GUARANA.Star := function( args )
    local recLieAlg, x, y, rep_info, star_method, basis, xx, yy, 
          wx, wy, max_w, xx_star_yy, x_star_y;

    recLieAlg := args[1];
    x := args[2];
    y := args[3];
    rep_info := args[4];
    if IsBound( args[5] ) then 
	star_method := args[5];
    else
	star_method := recLieAlg.star_method;
    fi;

    if star_method = "simple" then 
	if rep_info = "vec" then 
	    basis := Basis( recLieAlg.L );
	    xx := LinearCombination( basis, x );
	    yy := LinearCombination( basis, y );
            wx := GUARANA.WeightOfLieAlgElm( recLieAlg, xx );
	    wy := GUARANA.WeightOfLieAlgElm( recLieAlg, yy );
	    max_w := recLieAlg.max_weight;
	    xx_star_yy := GUARANA.Star_Simple( xx, yy, wx, wy, max_w, 
	                                     "strucConst"   );
	    return Coefficients( basis, xx_star_yy );
	elif rep_info = "elm" then 
            wx := GUARANA.WeightOfLieAlgElm( recLieAlg, x );
	    wy := GUARANA.WeightOfLieAlgElm( recLieAlg, y );
	    max_w := recLieAlg.max_weight;
	    x_star_y := GUARANA.Star_Simple( x, y, wx, wy, max_w, 
	                                     "strucConst"   );
	    return x_star_y;
	else 
	    Error( "wrong rep_info" );
	fi;
    elif star_method = "symbolic" then
        Error( "Sorry, method not implemented yet" );
    else
	Error( "Wrong exp_method" );
    fi;
end;

#############################################################################
##
#E
