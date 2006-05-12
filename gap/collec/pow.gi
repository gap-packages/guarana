#############################################################################
##
#W pow.gi              GUARANA package                     Bjoern Assmann
##
## Several methods for computing with  powers of automorphisms of N
## and powers of group elements of CN.
##
#H  @(#)$Id$
##
#Y 2006
##
##

#############################################################################
##
#F GUARANA.ConPow_Compute_Pi_q( recLieAlg, x, Phi, q )
##
## IN
## recLieAlg ............................ lie algebra record of a Lie algebra
##                                        L
## x .................................... element of the Lie algebra
##                                        given as coefficient vector
## Phi .................................. autmorphisms of the lie algebra
## q .................................... natural number ge 0
##
## OUT
## The coefficients of Pi_q where
## Pi_q = x^(Phi^(q-1)) * x^(Phi^(q-2)) * ... * x
##
GUARANA.ConPow_Compute_Pi_q := function( recLieAlg, x, Phi, q )
    local dim, bin, Phi_p, Pi_p, k, b_i, left, right, i;

    dim := recLieAlg.dim;

    # catch trivial case
    if q = 0 then 
	return List( [1..dim], x-> 0 );
    fi;
    if q < 0 then 
	Error( "q has to be greater equal 0" );
    fi;

    # binary representation of q 
    bin := CoefficientsQadic( q, 2 );
    bin := Reversed( bin );

    # compute Pi_p by variation of repeated squaring.
    # We use repeated squaring for Pi_p and for Phi_p = Phi^p
    Phi_p := Phi;
    Pi_p := x;
    k := Length( bin );
    for i in [2..k] do 
	b_i := bin[i];
	# use Pi_2p =( Pi_p Phi^p) * Pi_p 
	left := Pi_p*Phi_p;
	right := Pi_p;
	Pi_p := GUARANA.Star( [recLieAlg, left, right, "vec" ] );
	Phi_p := Phi_p * Phi_p;
	
	if b_i = 1 then 
	    ## use Pi_2p+1 = (Pi_2p Phi)* x 
	    left := Pi_p*Phi;
	    right := x;
	    Pi_p := GUARANA.Star( [recLieAlg, left,right, "vec"] );
	fi;
    od;

    return Pi_p;
end;

#############################################################################
##
#F GUARANA.ConPow_Compute_pi_q( recLieAlg, g, Phi, q )
##
## recLieAlg ........................ record of a Lie algebra L
## g ................................ group element of Exp( L ) given 
##                                    by exp vector.
## Phi .............................. Lie algebra automorphism of L
## q ................................ natural number ge 0
## 
## OUT 
## The exp vector of pi_q where 
## pi_q = n^(phi^(q-1)) * .... * n^phi * n 
## where phi is the group aut corresponding to Phi.
##
GUARANA.ConPow_Compute_pi_q := function( recLieAlg, g, Phi, q )
    local x, Pi_q, pi_q;

    # comput lie elment 
    x := GUARANA.AbstractLog( [recLieAlg, g, "vecByVec"] );

    # deal with the consectuive powers there
    Pi_q := GUARANA.ConPow_Compute_Pi_q( recLieAlg, x, Phi, q );

    # compute corresponding group elm
    pi_q := GUARANA.AbstractExp( [recLieAlg, Pi_q, "vecByVec"] );

    return pi_q;
end;

## COMMENT
## g^-1 = ( c(g) n(g) )^-1
##      =  c(g)^-1  ( n(g)^-1 )^-1
##
GUARANA.CN_Inversion := function( malcevRec, exp_g )
    local exp_g_cut, c_g, n_g, id_C, c_g_C, log_c_g, log_inv_c_g, 
          log_n_g, log_inv_n_g, log, c_inv_g_C, log_c_inv_g, log_t_LC, 
	  log_t, log_n_inv_g, n_inv_g, c_inv_g, i;

    # test whether g in CN
    exp_g_cut := GUARANA.CutExpVector( malcevRec, exp_g );
    if exp_g_cut[1] <> exp_g_cut[1]*0 then
	Error( "g is not in CN" );
    fi;

    c_g := exp_g_cut[2];
    n_g := exp_g_cut[3];

    # get exp of c(g) with respect to the pcp of CC
    id_C := List( [1..malcevRec.recL_CC.dim], x-> 0 );
    c_g_C := id_C + c_g;

    # map c(g) to L(C) and invert it, i.e. compute log( c^-1)
    log_c_g := GUARANA.AbstractLog( [malcevRec.recL_CC, c_g_C, "vecByVec"] );
    log_inv_c_g := - log_c_g;

    # map n(g) to L(N) and invert it, i.e. compute log( n^-1 )
    log_n_g := GUARANA.AbstractLog( [malcevRec.recL_NN, n_g, "vecByVec"] );
    log_inv_n_g := - log_n_g;

    # apply c^-1 to log( n^-1 )
    log := log_inv_n_g;
    for i in malcevRec.indeces[2] do
	log := log*malcevRec.lieAuts[i]^(-exp_g[i]);
    od;

    # compute tail of log( c^-1) with respect to basis of L(C)
    c_inv_g_C := - c_g_C;  # c(g^-1) with respect to pcs of CC
    log_c_inv_g := GUARANA.AbstractLog( [malcevRec.recL_CC, c_inv_g_C,
                                           "vecByVec" ] );
    log_t_LC := GUARANA.Star( [malcevRec.recL_CC, log_c_inv_g, log_inv_c_g,
                               "vec" ] );
    # check wheter log(t) has the right form
    if log_t_LC{[1..malcevRec.lengths[2]]} <> 
       0*log_t_LC{[1..malcevRec.lengths[2]]} then 
        Error( "log(t) does not have the correct form" );
    fi;
 
    # compute log(t) with respect to the basis of L(N)
    log_t := GUARANA.MapFromLCcapNtoLN( malcevRec, log_t_LC );
     
    # compute log( n(g^-1) )
    log_n_inv_g := GUARANA.Star( [malcevRec.recL_NN, log_t, log, "vecByVec"] );
    
    # compute n(g^-1)
    n_inv_g := GUARANA.AbstractExp( [malcevRec.recL_NN, log_n_inv_g, 
                                     "vecByVec" ] );

    # get exp of c(g^-1)
    c_inv_g := -c_g;

    return Concatenation( c_inv_g, n_inv_g );
end;

## Test function for computing inverses

GUARANA.CN_Power := function( malcevRec, g, pow )
    
end;

#############################################################################
##
#E
