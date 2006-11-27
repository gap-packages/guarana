#############################################################################
##
#W collec.gi               GUARANA package                     Bjoern Assmann
##
## Methods for  symbolic collection ( in CN ).
##
#H  @(#)$Id$
##
#Y 2006
##
##

# TODO
# - Test the duSautoy function via the evaluation thing.
# - write a little function that gives you the functions 
#   and maybe some more information about it. 

GUARANA.SC_GetOmegeVars := function( symCol )
    local omega_var_recs, k, omega_vars, i;

    omega_var_recs := symCol.omega_var_recs;
    k := Length( omega_var_recs );
    omega_vars := [];
    for i in [1..k] do 
        Append( omega_vars, omega_var_recs[i].variables );
    od;
    return omega_vars;
end;

## OUT
## computes the values of the omega variables 
## for given (right) multiplying elment h.
##
GUARANA.SC_ValsOfOmegaVars := function( malCol, exp_h )
    local omega_var_recs, k, vals, omega_rec, eigenvalues, ll, i;

    omega_var_recs := malCol!.symCol.omega_var_recs;
    k := Length( omega_var_recs );
    vals := [];
    for i in [1..k] do 
        omega_rec := omega_var_recs[i];
        eigenvalues := omega_rec.eigenvalues;
        ll := List( eigenvalues, x-> x^exp_h[i] );
        Append( vals, ll );
    od;
    return vals;
end;

## IN:
## exp_g, exp_h .................... exponent vectors of grp elements
##                                   in CN.
## I think rational entries could be allowed.
## In the moment we have the problem of computing a^(1/2)
## where "a" is an elm of a number field.
## OUT:
## exponent vector of g*h, computed via the du Sautoy functions.
##
GUARANA.SC_Evaluate := function( malCol, exp_g, exp_h )
    local one, exp_g_emb, exp_h_emb, left_vars, right_vars, omega_vars, 
          indets, left_vals, right_vals, omega_vals, vals, n, 
          collFuncs, exp_gh, pol, r, exp_gh_rat, rep, l, i;

    # map the entries of exp_g, exp_h in the extension field
    one := One( malCol!.symCol.splitField );
    exp_g_emb := exp_g * one;
    exp_h_emb := exp_h * one;

    # get all indeterminants
    left_vars := malCol!.symCol.left_vars;
    right_vars := malCol!.symCol.right_vars;
    omega_vars := GUARANA.SC_GetOmegeVars( malCol!.symCol );
    indets := Concatenation( left_vars, right_vars, omega_vars );

    # get all values
    left_vals := exp_g_emb;
    right_vals := exp_h_emb;
    omega_vals := GUARANA.SC_ValsOfOmegaVars( malCol, exp_h );
    vals := Concatenation( left_vals, right_vals, omega_vals );
    

    # compute exp vector of g*h 
    n := Length( exp_g );
    collFuncs := malCol!.symCol.collFuncs;
    exp_gh := [];
    for i in [1..n] do 
        pol := collFuncs[i]; 
        r := Value( pol, indets, vals );
        Add( exp_gh, r ); 
    od;

    # transfrom the entries of exp_gh to GAP rationals
    exp_gh_rat := [];
    for i in [1..Length( exp_gh )] do 
        rep := ExtRepOfObj( exp_gh[i] );
        # test whether entry is purely rational
        l := Length( rep );
        if not rep{[2..l]} = 0* rep{[2..l]} then
            Error( "exponent vector is not purely rational" );
        fi;
        Add( exp_gh_rat, rep[1] );
    od;

    return exp_gh_rat;
end;

GUARANA.SC_CN_Multiplication := function( g, h )
    local malCol, exp_g, exp_h, exp_gh;

    malCol := g!.malCol;
    exp_g := Exponents( g );
    exp_h := Exponents( h );
    exp_gh := GUARANA.SC_Evaluate( malCol, exp_g, exp_h );
    return MalcevCNElementByExponents( malCol, exp_gh ); 
end;

GUARANA.SC_TestCollector := function( malCol, noTests, range )
    local g, h, gh, exp_gh_1, exp_gh_2, i;

    for i in [1..noTests] do 
        g := Random( malCol, "CN", range );
        h := Random( malCol, "CN", range );

        SetMultiplicationMethod( malCol, "symbolic" );
        gh := g*h;
        exp_gh_1 := Exponents( gh );

        SetMultiplicationMethod( malCol, "standard" );
        gh := g*h;
        exp_gh_2 := Exponents( gh );

        if exp_gh_1 <> exp_gh_2 then 
            Error( "exp vectors don't coincide" );
        fi;
    od;
end;


if false then 
    mc := ExamplesOfSomeMalcevCollectors(1 );   
    
    mc := ExamplesOfSomeMalcevCollectors(2 );   

    mc := ExamplesOfSomeMalcevCollectors(3 );   

    mc := GUARANA.MalcevColl_F_nc_Aut1( 2, 2 );

    mc := GUARANA.MalcevColl_F_nc_Aut1( 2, 3);

    mc := GUARANA.MalcevColl_F_nc_Aut1( 2, 5);


    mc := GUARANA.MalcevColl_F_nc_Aut1( 2, 2 );

    AddSymbolicCollector( mc );
    SetMultiplicationMethod( mc, "symbolic" );

    g := Random( mc, "CN" );
    h := Random( mc, "CN" );
    GUARANA.SC_CN_Multiplication( g, h );

    r := g*h;

    GUARANA.SC_TestCollector( mc, 20, 10 );
fi;


#############################################################################
##
#E
