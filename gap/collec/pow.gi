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
#F GUARANA.ConPow_Compute_Pi_q( recLieAlg, x, myPhi, q )
##
## IN
## recLieAlg ............................ lie algebra record of a Lie algebra
##                                        L
## x .................................... element of the Lie algebra
##                                        given as coefficient vector
## myPhi .................................. autmorphisms of the lie algebra
## q .................................... natural number ge 0
##
## OUT
## The coefficients of Pi_q where
## Pi_q = x^(Phi^(q-1)) * x^(Phi^(q-2)) * ... * x
##
GUARANA.ConPow_Compute_Pi_q := function( recLieAlg, x, myPhi, q )
    local dim, bin, myPhi_p, Pi_p, k, b_i, left, right, i;

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
    myPhi_p := myPhi;
    Pi_p := x;
    k := Length( bin );
    for i in [2..k] do 
	b_i := bin[i];
	# use Pi_2p =( Pi_p Phi^p) * Pi_p 
	left := Pi_p*myPhi_p;
	right := Pi_p;
	Pi_p := GUARANA.Star( [recLieAlg, left, right, "vec" ] );
	# use Phi_2p = ( Phi_p )^2
	myPhi_p := myPhi_p * myPhi_p;
	
	if b_i = 1 then 
	    # use Pi_2p+1 = (Pi_2p Phi)* x 
	    left := Pi_p*myPhi;
	    right := x;
	    Pi_p := GUARANA.Star( [recLieAlg, left,right, "vec"] );
	    # use Phi_2p+1 =  Phi_2p  * Phi
	    myPhi_p := myPhi_p * myPhi;
	fi;
    od;

    return Pi_p;
end;

#############################################################################
##
#F GUARANA.ConPow_Compute_pi_q( recLieAlg, g, myPhi, q )
##
## recLieAlg ........................ record of a Lie algebra L
## g ................................ group element of Exp( L ) given 
##                                    by exp vector.
## myPhi .............................. Lie algebra automorphism of L
## q ................................ natural number ge 0
## 
## OUT 
## The exp vector of pi_q where 
## pi_q = n^(phi^(q-1)) * .... * n^phi * n 
## where phi is the group aut corresponding to Phi.
##
GUARANA.ConPow_Compute_pi_q := function( recLieAlg, g, myPhi, q )
    local x, Pi_q, pi_q;

    # comput lie elment 
    x := GUARANA.AbstractLog( [recLieAlg, g, "vecByVec"] );

    # deal with the consectuive powers there
    Pi_q := GUARANA.ConPow_Compute_Pi_q( recLieAlg, x, myPhi, q );

    # compute corresponding group elm
    pi_q := GUARANA.AbstractExp( [recLieAlg, Pi_q, "vecByVec"] );

    return pi_q;
end;

#############################################################################
##
#F GUARANA.CN_Inversion( malcevRec, exp_g )
##
## malcevRec ............................... Malcev record
## exp_g ................................... exponent vector of g in CN
##
## OUT
## Exponent vector of g^-1
##
## COMMENT
## g^-1 = ( c(g) n(g) )^-1
##      =  c(g)^-1  ( n(g)^-1 )^-1
##
GUARANA.CN_Inversion := function( malcevRec, exp_g )
    local exp_g_cut, c_g, n_g, id_C, c_g_C, log_c_g, log_inv_c_g, 
          log_n_g, log_inv_n_g, log, c_inv_g_C, log_c_inv_g, log_t_LC, 
	  log_t, log_n_inv_g, n_inv_g, c_inv_g, i, f_inv_g;

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
    log_t_LC := GUARANA.Star( [malcevRec.recL_CC, -log_c_inv_g, log_inv_c_g,
                               "vec" ] );

    # check wheter log(t) has the right form
    if log_t_LC{[1..malcevRec.lengths[2]]} <> 
       0*log_t_LC{[1..malcevRec.lengths[2]]} then 
        Error( "log(t) does not have the correct form" );
    fi;
 
    # compute log(t) with respect to the basis of L(N)
    log_t := GUARANA.MapFromLCcapNtoLN( malcevRec, log_t_LC );
     
    # compute log( n(g^-1) )
    log_n_inv_g := GUARANA.Star( [malcevRec.recL_NN, log_t, log, "vec"] );
    
    # compute n(g^-1)
    n_inv_g := GUARANA.AbstractExp( [malcevRec.recL_NN, log_n_inv_g, 
                                     "vecByVec" ] );

    # get exp of c(g^-1)
    c_inv_g := -c_g;

    # get exp of finite part
    f_inv_g := List( [1..malcevRec.lengths[1]], x-> 0 );

    return Concatenation( f_inv_g, c_inv_g, n_inv_g );
end;

#############################################################################
##
#F GUARANA.CN_PositivePower( malcevRec, exp_g, q )
##
## IN
## malcevRec ...................... Malcev record
## exp_g  ......................... exponent vector of g in CN
## q .............................. natural number ge 0
##
## OUT
## Exponent vector of g^q.
##
## COMMENT
## g^q = ( c(g)n(g) )^q = c(g)^q  pi_q
##
GUARANA.CN_PositivePower := function( malcevRec, exp_g, q )
    local exp_g_cut, c_g, n_g, id_C, c_g_C, log_c_g, log_c_g_q, divider, 
          log_divider, log_t_LC, log_t, myPhi, log_n, Pi_q, log_n_g_q, 
	  n_g_q, c_g_q, f_g_q, i;

    # test input and catch trivial case
    if q < 0 then 
	Error( "Power q has to positive" );
    elif q = 0 then 
	return 0*exp_g;
    elif q = 1 then 
	return exp_g;
    fi;

    # test whether g in CN
    exp_g_cut := GUARANA.CutExpVector( malcevRec, exp_g );
    if exp_g_cut[1] <> exp_g_cut[1]*0 then
	Error( "g is not in CN" );
    fi;

    c_g := exp_g_cut[2];
    n_g := exp_g_cut[3];

    # get exponents of c(g) with respect to the pcp of CC
    id_C := List( [1..malcevRec.recL_CC.dim], x-> 0 );
    c_g_C := id_C + c_g;

    # map c(g) to L(C) and power it, i.e. compute log( c^q )
    log_c_g := GUARANA.AbstractLog( [malcevRec.recL_CC, c_g_C, "vecByVec"] );
    log_c_g_q := q * log_c_g;

    # compute tail of log( c^q) with respect to basis of L(C)
    # i.e. we divide off in the lie algebra the part that belongs to CN/N
    divider := q* c_g_C;  # c(g^q) with respect to pcs of CC
    log_divider := GUARANA.AbstractLog( [malcevRec.recL_CC, divider,
                                           "vecByVec" ] );
    log_t_LC := GUARANA.Star( [malcevRec.recL_CC, -log_divider, log_c_g_q,
                               "vec" ] );

    # check wheter log(t) has the right form
    if log_t_LC{[1..malcevRec.lengths[2]]} <> 
       0*log_t_LC{[1..malcevRec.lengths[2]]} then 
        Error( "log(t) does not have the correct form" );
    fi;
 
    # compute log(t) with respect to the basis of L(N)
    log_t := GUARANA.MapFromLCcapNtoLN( malcevRec, log_t_LC );

    # compute Phi
    myPhi := IdentityMat( malcevRec.recL_NN.dim );
    for i in malcevRec.indeces[2] do
	myPhi := myPhi*malcevRec.lieAuts[i]^(exp_g[i]);
    od;

    # compute log( n(g) )
    log_n := GUARANA.AbstractLog( [malcevRec.recL_NN, n_g, "vecByVec" ] );   

    # compute Pi_q
    Pi_q := GUARANA.ConPow_Compute_Pi_q( malcevRec.recL_NN, log_n, myPhi, q );

    # compute n( g^q ) 
    log_n_g_q := GUARANA.Star( [malcevRec.recL_NN, log_t, Pi_q, "vec" ] );
    n_g_q := GUARANA.AbstractExp( [malcevRec.recL_NN, log_n_g_q, "vecByVec"]);

    # get c(g^q)
    c_g_q := q * c_g;

    # get exp of finite part
    f_g_q := List( [1..malcevRec.lengths[1]], x-> 0 );

    return Concatenation( f_g_q, c_g_q , n_g_q );
end;

#############################################################################
##
#F GUARANA.CN_Power( malcevRec, exp_g, q )
##
## IN
## malcevRec ...................... Malcev record
## exp_g  ......................... exponent vector of g in CN
## q .............................. integer 
##
## OUT
## Exponent vector of g^q.
##
GUARANA.CN_Power := function( malcevRec, exp_g, q )
    local exp_inv_g;

    if q >= 0 then 
	return GUARANA.CN_PositivePower( malcevRec, exp_g, q );
    else 
	exp_inv_g := GUARANA.CN_Inversion( malcevRec, exp_g );
	return GUARANA.CN_PositivePower( malcevRec, exp_inv_g, -q );
    fi;
end;

#############################################################################
##
## Test functions
##
if false then 
    ll := GUARANA.SomePolyMalcevExams( 3 );
    R := GUARANA.InitialSetupCollecRecord( ll );;
    GUARANA.AddCompleteMalcevInfo( R );
    g := GUARANA.RandomGrpElm( [R,10, "CN"] );
    h := GUARANA.RandomGrpElm( [R,10, "CN"] );
fi;
GUARANA.Test_CN_Inversion := function( malcevRec, range )
    local g, exp_g, exp_inv_g, inv_g;

    # get random element in CN
    g := GUARANA.RandomGrpElm( [malcevRec, range, "CN", "elm" ]);
    exp_g := Exponents( g );

    # invert it with Mal'cev
    exp_inv_g := GUARANA.CN_Inversion( malcevRec, exp_g );

    # invert it in the group 
    inv_g := g^-1;

    # compare
    return Exponents( inv_g ) = exp_inv_g;
end;

GUARANA.Test_CN_Inversion2 := function( malcevRec, range )
    local exp_g, exp_inv_g, exp_inv_inv_g;

    # get random element in CN
    exp_g := GUARANA.RandomGrpElm( [malcevRec, range, "CN", "exp" ]);

    # invert it with Mal'cev
    exp_inv_g := GUARANA.CN_Inversion( malcevRec, exp_g );

    # invert it again 
    exp_inv_inv_g := GUARANA.CN_Inversion( malcevRec, exp_inv_g );

    # compare
    return  exp_inv_inv_g =  exp_g;
end;

if false then 
   exp_g := [ 0, 0, 0, 0, 1, 0, 0, 3, -1, 0, 0, 0, 0, 0, 0, 0, 0 ];
   g := GUARANA.GrpElmByExpsAndPcs( Pcp( R.G ), exp_g );
   q := 6;
   malcevRec := R;
fi;

GUARANA.Test_CN_PositivePower := function( malcevRec, range, q )
    local g, exp_g, exp_g_q, g_q;

    # get random element in CN
    g := GUARANA.RandomGrpElm( [malcevRec, range, "CN", "elm" ]);
    exp_g := Exponents( g );

    # power it with Mal'cev
    exp_g_q := GUARANA.CN_PositivePower( malcevRec, exp_g, q );
    Print( exp_g_q, "\n"  );

    # power it in the group 
    g_q  := g^q;
    Print( Exponents( g_q ), "\n" );

    # compare
    return Exponents( g_q ) = exp_g_q;
end;

GUARANA.Test_CN_Power := function( malcevRec, range, q )
    local g, exp_g, exp_g_q, g_q;

    # get random element in CN
    g := GUARANA.RandomGrpElm( [malcevRec, range, "CN", "elm" ]);
    exp_g := Exponents( g );

    # power it with Mal'cev
    exp_g_q := GUARANA.CN_Power( malcevRec, exp_g, q );
    Print( exp_g_q, "\n"  );

    # power it in the group 
    g_q  := g^q;
    Print( Exponents( g_q ), "\n" );

    # compare
    return Exponents( g_q ) = exp_g_q;
end;

#############################################################################
##
#E
