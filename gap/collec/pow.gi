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

## IN
## x .................................... Malcev Gen element        
## myPhi .................................autmorphisms of the lie algebra
## q .................................... natural number ge 0
##
## OUT
## The coefficients of Pi_q where
## Pi_q = x^(Phi^(q-1)) * x^(Phi^(q-2)) * ... * x
##
GUARANA.Compute_Pi_q := function( x, myPhi, q )
    local bin, myPhi_p, Pi_p, k, b_i, left, right, i;

    # catch trivial case
    if q = 0 then 
        return 0*x;
    elif q < 0 then 
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
        Pi_p := BCHStar( left, right );
	    # use Phi_2p = ( Phi_p )^2
	    myPhi_p := myPhi_p * myPhi_p;
	
	    if b_i = 1 then 
	        # use Pi_2p+1 = (Pi_2p Phi)* x 
	        left := Pi_p*myPhi;
	        right := x;
            Pi_p := BCHStar( left, right );
	        # use Phi_2p+1 =  Phi_2p  * Phi
	        myPhi_p := myPhi_p * myPhi;
	    fi;
    od;

    return Pi_p;
end;

GUARANA.LieAutBy_C_Elm := function( malCol, c )
    local lieAuts, mat, exps_c, i;
    
    lieAuts := malCol!.C_lieAuts;
    if Length( lieAuts ) = 0 then 
        return EmptyMatrix( 0 );
    else 
        mat := lieAuts[1]^0; 
        exps_c := Exponents( c );
        for i in [1..Length(exps_c) ] do
            mat := mat * (lieAuts[i]^exps_c[i]);
        od;
        return mat;
    fi; 
end;

##
## (cn)^q = c^q *  Pi_q
GUARANA.CN_PosPower := function( g, q )
    local c, n, c_new, n_new, myPhi;

    # catch trivial cases
    if q < 0 then 
	    Error( "Power q has to positive" );
    elif q = 0 then 
	    return 0*g;
    elif q = 1 then 
	    return g;
    fi;

    c := g!.c;
    n := g!.n;
    if c = 0*c then 
        c_new := c;
        n_new := q*n;
    elif n = 0*n then
        c_new := q*c;
        n_new := n;
    else
        # get lie algebra automorphism
        myPhi := GUARANA.LieAutBy_C_Elm( g!.malCol, c );
        n_new := GUARANA.Compute_Pi_q( n, myPhi, q );
        c_new := c^q;
    fi;
    return MalcevCNElementBy2GenElements( g!.malCol, c_new, n_new );
end; 

InstallOtherMethod( \^, 
               "for Malcev CN elments (Guarana)",
	           true,
	        [IsMalcevCNElement, IsInt ],
		0, 
function( g, n )
    local g_inv;
    if n < 0 then 
        g_inv := GUARANA.CN_Inverse( g );
        return GUARANA.CN_PosPower( g_inv, -n );
    else
        return GUARANA.CN_PosPower( g, n );
    fi;
end);

#############################################################################
##
## Test code

GUARANA.Test_CNInversion := function( malCol, range )
    local a, a_i, a_ii;

    # get random element in CN
    a := Random( malCol, "CN2", range );

    # invert it with Mal'cev
    a_i := a^-1;

    # invert it again 
    a_ii := a_i^-1;

    # compare
    return  a = a_ii;
end;

#############################################################################
##
#E
