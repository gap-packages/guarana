#############################################################################
##
#W pow.gi              GUARANA package                     Bjoern Assmann
##
## Several methods for computing with  powers of automorPhisms of N
##
#H  @(#)$Id$
##
#Y 2006
##
##


## IN
## recLieAlg ............................ lie algebra record of a Lie algebra
##                                        L
## x .................................... element of the Lie algebra
##                                        given as coefficient vector
## Phi .................................. autmorphisms of the lie algebra
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

GUARANA.ConPow_Compute_pi_q := function( recLieAlg, g, Phi, q )

end;

# TODO
# It might be useful to have a function that computes powers
# of g in CN but leaves the N-part as coefficient vector with respect to 
# L(N).
# Further it would be good to have a multiplication function in CN
# that can deal with this representation.

#############################################################################
##
#E
