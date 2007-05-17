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


## IN
## malCol ..................... Malcev collector
##    n ....................... Malcev gen elment of N
##    c ....................... Malcev gen elment of C
## 
## n^c given as Malcev gen elment of N
##
GUARANA.N_ConjugationByC_Elm := function( malCol, n, c )
    local coeffs, exps_c, lieAuts, i;

    # catch trivial case 
    exps_c := Exponents( c );
    if exps_c = 0* exps_c then 
        return n;
    fi;
  
    coeffs := Coefficients( n );
    n := Length( exps_c );
    lieAuts := malCol!.C_lieAuts;
    for i in [1..n] do
        coeffs := coeffs*( lieAuts[i]^exps_c[i] );
    od;
    return MalcevGenElementByCoefficients( malCol!.mo_NN, coeffs );
end;

## IN
## malCol
## n ........................ Malcev gen elment of N
## exps_g .................... exponent vector of group elment g of G.
##                            the vector is allowed to be shorter then
##                            the lengths of the pcs of G.
##
## n^g given as Malcev gen elment 
##
GUARANA.N_ConjugationByExp := function( malCol, n, exps_g )
    local coeffs, lieAuts, i;
  
    coeffs := Coefficients( n );
    n := Length( exps_g );
    lieAuts := malCol!.G_lieAuts;
    for i in [1..n] do
        coeffs := coeffs*( lieAuts[i]^exps_g[i] );
    od;
    return MalcevGenElementByCoefficients( malCol!.mo_NN, coeffs );
end;

## IN 
## g,h ..................... Malcev CN elemnts
## 
## OUT
## g*h computed with the normal malcev collector
##
GUARANA.CN_Multiplication := function( g, h )
    local malCol, c_new, n_c, n_new;

    malCol := g!.malCol;
    c_new := g!.c * h!.c;
    n_c := GUARANA.N_ConjugationByC_Elm( malCol, g!.n, h!.c );
    n_new := n_c * h!.n; 
    return MalcevCNElementBy2GenElements( malCol, c_new, n_new ); 
end;

#############################################################################
##
#M g* h  ................................... Product of Malcev CN elements
##
## Comment: 
## g h = c(g)n(g) c(h)n(h)
##     = c(g)c(h0 n(g)^c(h) n(h)
##
## TODO
## Install a variation for this for symbolic CN elements
##
InstallOtherMethod( \*, 
               "for Malcev CN elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevCNElement, IsMalcevCNElement ],
		0, 
function( g, h )
    local malCol;

    malCol := g!.malCol;
    if malCol!.mult_method = "standard" then 
        return GUARANA.CN_Multiplication( g, h );
    elif malCol!.mult_method = "symbolic" then 
        return GUARANA.SC_CN_Multiplication( g, h );
    else
        Error( "Wrong mulitplication method is specified" );
    fi;
end);

## Note that if "symbol" should be used as a collection method 
## for the whole group, then we need symbolic inversion in CN as well.
##
## So far only symbolic mulitplication in CN is implemented.
##
InstallMethod( SetMultiplicationMethod,
               "for Malcev collectors and strings (Guarana)", 
               true, 
               [IsMalcevCollectorRep, IsString ], 
               0,
function( malCol, s)
    local possible_methods;

    possible_methods := [ "standard", 
                          "symbolic" ];
    if s in possible_methods then 
        if s = "symbolic" then 
            if not IsBound( malCol!.symCol ) then 
                #TODO: timing would be nice.
                Info( InfoGuarana, 1, "Computing symbolic collector ...\n" );
                AddSymbolicCollector( malCol );
            fi;
        fi;
	    malCol!.mult_method := s;
    else
	    Error( "Wrong Multiplication method specified\n" );
    fi;
end );

##
## (cn)^-1 = n^-1 c^-1 = c^-1 (n^-1)^(c^-1)
##
GUARANA.CN_Inverse := function( g )
    local c, n, c_inv, n_inv, n_inv_c_inv;

    c := g!.c;
    n := g!.n;
    c_inv := c^-1;
    n_inv := n^-1;

    n_inv_c_inv := GUARANA.N_ConjugationByC_Elm( g!.malCol, n_inv, c_inv );
    return MalcevCNElementBy2GenElements( g!.malCol, c_inv, n_inv_c_inv );
end;

#############################################################################
##
#F GUARANA.MO_CconjF_ElmConjByFiniteElm( malCol, exp_c, exp_f )
##
## IN
## malCol ........................ Malcev collector
## exp_c ............................ short exponent vector (with respect 
##                                    to pcs of CN/N) of an element 
##                                    in C.
##                                    c = c_1^x_1 ... c_k^x_k 
##                                    where k is the rank of CN/N.
## exp_f  ................................ short expvector (with respect 
##                                         to the pcs of G/CN ) of an elment
##                                         f in G, that is a product of 
##                                         generators of the finite part of the
##                                         pcs. 
## OUT
## Exponent vector of c^f
##
## COMMENT
## c^f = (c_1^x_1 ... c_k^x_k )^f
##     = (c_1^x_1)^f ... (c_k^x_k)^f
## 
GUARANA.MO_CconjF_ElmConjByFiniteElm := function( malCol, exp_c, exp_f )
    local res, x_i, c_i_f, i;

    res := GUARANA.CN_Identity( malCol );   
    for i in [1..Length( exp_c )] do
	    x_i := exp_c[i];
	    if x_i <> 0 then 
	        # look up c_i^f
            c_i_f := GUARANA.MO_CconjF_LookUp( malCol, i, exp_f );

            res := res*(c_i_f^x_i);
	    fi;
    od;
    return res;
end;

#############################################################################
##
#F GUARANA.MO_CN_ConjugationByFiniteElm( malCol, exp_g, exp_f )
##
## IN
## g      ................................ Malcev CN Element 
## exp_f  ................................ short expvector (with respect 
##                                         to the pcs of G/CN ) of an elment
##                                         f in G, that is a product of 
##                                         generators of the finite part of the
##                                         pcs. 
##
## OUT
## Malcev CN element  g^f
##
## COMMENT
## g^f = ( c(g) n(g) )^f
##     =  c(g)^f n(g)^f
##
GUARANA.MO_CN_ConjugationByFiniteElm := function( g, exp_f )
    local malCol, exps, hlCN, exps_c_g, c_g_f, n_g_f, c_new, n_new;

    # catch trivial case G/CN = 1
    if Length( exp_f ) = 0 then 
        return g;
    fi;

    # get part of the exponent vector of g that corresponds to CN/N
    malCol := g!.malCol;
    exps := Exponents( g );
    hlCN := malCol!.lengths[2];
    exps_c_g := exps{[1..hlCN]};

    # compute c(g)^f
    c_g_f := GUARANA.MO_CconjF_ElmConjByFiniteElm( malCol, exps_c_g, exp_f);

    # compute n(g)^f
    n_g_f := GUARANA.N_ConjugationByExp( malCol, g!.n, exp_f );

    c_new := c_g_f!.c;
    n_new := c_g_f!.n * n_g_f;
    return MalcevCNElementBy2GenElements( malCol, c_new, n_new );
end;

#############################################################################
##
#F GUARANA.MO_G_Collection( g, h )
##
## COMMENT
## g h = f(g)c(g)n(g) f(h)c(h)n(h) 
##     = f(g)f(h)        (c(g)n(g))^f(h)    c(h)n(h)
##     = f(gh) u         v                  w
## 
GUARANA.MO_G_Collection := function( g, h )
    local malCol, f_g, f_h, a, f_gh, u, b, v, w, uvw;

    # setup
    malCol := g!.malCol;

    # get f(gh) and u
    f_g := g!.exps_f;
    f_h := h!.exps_f;
    a := GUARANA.MO_G_CN_LookUpProduct( malCol, f_g, f_h );
    f_gh := a!.exps_f;
    u := a!.cn_elm;

    # compute v
    b := g!.cn_elm;
    v := GUARANA.MO_CN_ConjugationByFiniteElm( b, f_h );

    # get w 
    w := h!.cn_elm;

    # compute uvw
    uvw := u*v*w;

    return MalcevGElementByCNElmAndExps( malCol, f_gh, uvw );
end;

#############################################################################
##
#M g* h  ................................... Product of Malcev G elements
##
##
InstallOtherMethod( \*, 
               "for Malcev G elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevGElement, IsMalcevGElement ],
		0, 
function( g, h )
    return GUARANA.MO_G_Collection( g, h );
end);


## COMMENT
## g^-1 = (fcn)^-1 
##      = (cn)^-1 f^-1
##
## f^-1 can be precomputed.
##
GUARANA.G_Inversion := function( g )
    local malCol, cn, cn_inv, exps_f, f_inv, c_inv2, g_inv;

    malCol := g!.malCol;

    # compute (cn)^-1
    cn := g!.cn_elm;
    cn_inv := cn^-1; 

    # get f^-1
    exps_f := g!.exps_f;
    f_inv := GUARANA.MO_F_LookupInverse( malCol, exps_f ); 

    # compute (cn)^-1 f^-1 = 
    c_inv2 := MalcevGElementByCNElmAndExps( malCol, 0*exps_f, cn_inv );
    g_inv := c_inv2 * f_inv;

    return g_inv;
end;

InstallOtherMethod( Inverse, 
               "for Malcev G elments (Guarana)",
	       true,
	        [IsMalcevGElement ],
		0, 
function( g )
    return GUARANA.G_Inversion( g );
end);

InstallOtherMethod( \^, 
               "for Malcev G elments and integers (Guarana)",
	           true,
	        [IsMalcevGElement, IsInt ],
		0, 
function( g, d )
    local  res;

    # first catch the trivial cases
    if d = 0 then 
        return  GUARANA.G_Identity( g!.malCol );
    elif d = 1 then 
        return g;
    elif d = -1 then
        return Inverse(g);
    fi;

    # set up for computation
    if d < 0 then
        g := Inverse(g);
        d := -d;
    fi;

    # compute power
    res := g^0;
    while d > 0 do
        if d mod 2 = 1 then   res := res * g;   fi;

        d := QuoInt( d, 2 );
        if d <> 0 then   g := g * g;    fi;
    od;

    return res;
end);
#############################################################################
##
#E
