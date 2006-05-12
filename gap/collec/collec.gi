#############################################################################
##
#W collec.gi              GUARANA package                     Bjoern Assmann
##
## Methods for collection in CN.
## Methods for collection in G using the collection in CN
## as a black box algorithm.
##
#H  @(#)$Id$
##
#Y 2006
##
##

#############################################################################
##
#F GUARANA.CutExpVector( malcevRec, exp )
##
## IN
## malcevRec ............................ malcev record
## exp .................................. exponent vector of an elment
##                                        g in G.
##
## OUT
## [exp_f,exp_c,exp_n] where exp is the concantenation of these and
## exp_f specifies the exponents of the finite part G/CN
## exp_fa specifies the exponents of the free abelian part CN/N
## exp_n specifies the exponents of the nilpotent part  N.
##
GUARANA.CutExpVector := function( malcevRec, exp )
    local indeces, exp_f, exp_fa, exp_n;
    indeces := malcevRec.indeces;
    exp_f := exp{indeces[1]};
    exp_fa := exp{indeces[2]};
    exp_n := exp{indeces[3]};
    return [exp_f, exp_fa, exp_n];
end;

#############################################################################
##
#F GUARANA.RandomGrpElm( args )
##
## IN
## args[1]=malcevRec ............... malcev record
## args[2]=range ................... range of random element
## args[4]=form .................... string that is either 
##                                   "exp" or "elm" depending on 
##                                   how the random element should be
##                                   given. 
## args[3]=grp ..................... optional string specifying whether
##                                   and elment in G,CN or N should be given.
##                                   Default is "G".
##
## OUT
## Random group elment of N,CN or G given by a exponent vector or 
## as actual group element
##
GUARANA.RandomGrpElm := function( args )
    local malcevRec, range, form, grp, hl, domain, indeces, rels, exp, i;

    malcevRec := args[1];
    range := args[2];
    if IsBound( args[3]  ) then 
	grp := args[3];
    else
	grp := "G";
    fi;
    if IsBound( args[4] ) then 
        form := args[4];
    else
	form := "exp";
    fi;

    hl := HirschLength( malcevRec.G  );
    domain := [-range..range];
    indeces := malcevRec.indeces;
    rels := RelativeOrdersOfPcp( Pcp( malcevRec.G ));

    exp := List( [1..hl], x-> 0 );
    if grp = "N" then 
	for i in indeces[3] do
	    exp[i] := Random( domain );
	od;
    elif grp = "CN" then 
	for i in Concatenation( indeces[2], indeces[3] ) do
	    exp[i] := Random( domain );
	od;
    elif grp = "G" then 
	for i in Concatenation( indeces[2], indeces[3] ) do
	    exp[i] := Random( domain );
	od;
	for i in indeces[1] do
	    exp[i] := Random( domain ) mod rels[i];
	od;
    else
	Error( " " );
    fi;
	    
    if form = "exp" then 
	return exp;
    else 
	return GUARANA.GrpElmByExpsAndPcs( exp );
    fi;
end;

## Example
## ll := GUARANA.SomePolyMalcevExams( 3 );
## R := GUARANA.InitialSetupCollecRecord( ll );
## GUARANA.AddCompleteMalcevInfo( R );
## g := GUARANA.RandomGrpElm( [R,10, "CN"] );
## h := GUARANA.RandomGrpElm( [R,10, "CN"] );

#############################################################################
##
#F GUARANA.Collection_CN( malcevRec, exp_g, exp_h )
#F GUARANA.Collection_CN_Star( malcevRec, exp_g, exp_h )
#F GUARANA.Collection_CN_DT( malcevRec, exp_g, exp_h )
##
## IN
## malcevRec ........................ malcev record
## exp_g,exp_h  ..................... exponent vectors of group elements
##                                    g,h in CN 
## 
## OUT 
## The exponent vector of g*h
##
## COMMENT
## g h = c(g)n(g) c(h)n(h)
##     = c(g)c(h) n(g)^c(h) n(h)
##     = c(gh) t n(g)^c(h) n(h)
##
## In Collection_CN_Star the computation in N are done in L(N) 
## via the operation star. 
## 
GUARANA.CN_Collection_Star := function( malcevRec, exp_g, exp_h )
    local exp_g_cut, exp_h_cut, c_g, c_h, n_g, n_h, log, id_C, c_g_C, 
          c_h_C, c_gh, c_gh_C, log_c_g, log_c_h, log_c_g_c_h, log_c_gh, 
	  log_t_LC, log_t, log_n_h, y, z, n_gh, f_gh, exp_gh, i;

    # test whether g,h in CN
    exp_g_cut := GUARANA.CutExpVector( malcevRec, exp_g );
    exp_h_cut := GUARANA.CutExpVector( malcevRec, exp_h );
    if exp_g_cut[1] <> exp_g_cut[1]*0 then
	Error( "g is not in CN" );
    elif exp_h_cut[1] <> exp_h_cut[1]*0 then
	Error( "h is not in CN" );
    fi;

    c_g := exp_g_cut[2];
    c_h := exp_h_cut[2];
    n_g := exp_g_cut[3];
    n_h := exp_h_cut[3];

    # map n(g) to L(N) and apply the automorphism corresponding to c(h)
    log := GUARANA.AbstractLog( [malcevRec.recL_NN, n_g, "vecByVec"] );
    for i in malcevRec.indeces[2] do
        log := log*malcevRec.lieAuts[i]^exp_h[i];
    od;

    # get exp vector of c(g), c(h), c(gh) with repsect to the pcp of CC
    id_C := List( [1..malcevRec.recL_CC.dim, x-> 0 );
    c_g_C := id_C + c_g;
    c_h_C := id_C + c_h;
    c_gh := c_g + c_h;
    c_gh_C := id_C + c_gh_C;

    # compute log(t) with respect to the basis of L(C)
    log_c_g := GUARANA.AbstractLog( [malcevRec.recL_CC, c_g_C, "vecByVec"] );
    log_c_h := GUARANA.AbstractLog( [malcevRec.recL_CC, c_h_C, "vecByVec"] );
    log_c_g_c_h := GUARANA.Star( [malcevRec.recL_CC, log_c_g, log_c_h,
                                  "vec" ] );
    log_c_gh := GUARANA.AbstractLog( [malcevRec.recL_CC, c_gh_C, "vecByVec"] );
    log_t_LC := GUARANA.Star( [malcevRec.recL_CC, -log_c_gh, log_c_g_ch,
                               "vec" ] );

    # check wheter log(t) has the right form
    if log_t_LC{[1..malcevRec.lengths[2]]} <> 
       0*log_t_LC{[1..malcevRec.lengths[2]} then 
        Error( "log(t) does not have the correct form" );
    fi;
    
    # compute log(t) with respect to the basis of L(N)
    log_t := GUARANA.MapFromLCcapNtoLN( malcevRec, log_t_LC );

    # compute log( n(h) )
    log_n_h := GUARANA.AbstractLog( [malcevRec.recL_NN, n_h, "vecByVec"] );

    # compute z=log(t)*log(n(g)^c(h))*log(n(h))
    y := GUARANA.Star( [malcevRec.recL_NN, log_t, log, "vec" ] );
    z := GUARANA.Star( [malcevRec.recL_NN, y, log_n_h, "vec" ] );
    
    # compute n(gh)
    n_gh := GUARANA.AbstractExp( [malcevRec.recL_NN, z, "vecByVec"] );

    # get exp of finite part
    f_gh := List( [1..malcevRec.lengths[1]], x-> 0 );

    exp_gh := Concatenation( f_gh,c_gh, n_gh );
    return exp_gh;
end;

GUARANA.CN_Collection_DT := function( malcevRec, exp_g, exp_h )

end;

GUARANA.CN_Collection := function( malcevRec, exp_g, exp_h )
    local method;

    method := malcevRec.collCN_method;
    if method = "star" then 
	return GUARANA.Collection_CN_Star( malcevRec, exp_g, exp_h );
    elif method = "deepThought" then 
	Error( "Sorry not implemented yet" );
    else
	Error( " " );
    fi;
end;

#############################################################################
##
#E
