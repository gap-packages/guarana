#############################################################################
##
#W setup.gi               GUARANA package                     Bjoern Assmann
##
## Methods for symbolic log and exp. 
##
#H  @(#)$Id$
##
#Y 2006
##
##

#############################################################################
##
##
InstallMethod( SetLogAndExpMethod, 
               "for Malcev objects and strings (Guarana)", 
               true, 
               [IsMalcevObjectRep, IsString ], 
               0,
function( malcevObject, s)
    local possible_methods;
    possible_methods := [ "pols", "simple" ];
    if s in possible_methods then 
	if s = "pols" then 
	    if not IsBound( malcevObject!.recLogPols) then 
        # we assume that if recLogPols is not set, then the exp pols are 
        # not known as well.
		Info( InfoGuarana, 1, "Computing Log and Exp Polynomials ...\n" );
		AddLogAndExpPolynomials( malcevObject );
	    fi;
	fi;
	malcevObject!.log_method := s;
    else
	Error( "Wrong Log method specified\n" );
    fi;
end );

InstallMethod( LogMethod, 
               "for Malcev objects (Guarana)", 
               true, 
               [IsMalcevObjectRep ], 
               0,
function( malcevObject)
    return malcevObject!.log_method;
end );

InstallMethod( ExpMethod, 
               "for Malcev objects (Guarana)", 
               true, 
               [IsMalcevObjectRep ], 
               0,
function( malcevObject)
    return malcevObject!.exp_method;
end );

InstallMethod( SetStarMethod,
               "for Malcev objects and strings (Guarana)", 
               true, 
               [IsMalcevObjectRep, IsString ], 
               0,
function( malcevObject, s)
    local possible_methods;
    possible_methods := [ "pols", "simple" ];
    if s in possible_methods then 
	    if s = "pols" then 
	        if not IsBound( malcevObject!.recStarPols) then 
		        Info( InfoGuarana, 1, "Computing Star Polynomials ...\n" );
		        AddStarPolynomials( malcevObject );
	        fi;
	    fi;
	    malcevObject!.star_method := s;
    else
	    Error( "Wrong Star method specified\n" );
    fi;
end );

InstallMethod( StarMethod, 
               "for Malcev objects (Guarana)", 
               true, 
               [IsMalcevObjectRep ], 
               0,
function( malcevObject)
    return malcevObject!.star_method;
end );

InstallMethod( SetMultiplicationMethod,
               "for Malcev objects and strings (Guarana)", 
               true, 
               [IsMalcevObjectRep, IsString ], 
               0,
function( malcevObject, s)
    local possible_methods;
    possible_methods := [ GUARANA.MultMethodIsStar, 
                          GUARANA.MultMethodIsCollection ];
    if s in possible_methods then 
	    malcevObject!.mult_method := s;
    else
	    Error( "Wrong Multplication method specified\n" );
    fi;
end );

InstallMethod( MultiplicationMethod, 
               "for Malcev objects (Guarana)", 
               true, 
               [IsMalcevObjectRep ], 
               0,
function( malcevObject)
    return malcevObject!.mult_method;
end );

#############################################################################
##
#F GUARANA.RationalVariableList( n, st )
##
## IN
## n ..... number of variables
## st ..... string for example  "x"
##
## OUT
## A list of rational n variables, called for example 
## x1,x2,...
##
GUARANA.RationalVariableList := function( n, st )
     return List(
                  List( [1..n], i->Concatenation( st, String(i) ) ),
                  x->Indeterminate( Rationals, x : new ) );
end;

#############################################################################
##
#F GUARANA.MO_AddStarPolynomialsSingleVersusGeneric( malcevObject )
##
## EFFECT
## recStarPolsSingVGen is added to the Malcev object
##
## For Log and Exp I just need 
## log x_i * sum_{j=i}^l beta_j log x_j
## So these are the polynomials I should save. It is NOT necessary to
## save x_i against a generic element of L.
## TODO
## compute polynomials which can used to speed up star operation
## How can I speed up this function ??
##
## Example
## 
##     malObjs := GUARANA.Get_FNG_MalcevObjects( 2, 9 );
##     malObj := malObjs[7];
##     GUARANA.MO_AddStarPolynomials( malObj );
##
GUARANA.MO_AddStarPolynomialsSingleVersusGeneric := function( malcevObject )
    local n, vars_x, vars_y, star_pols, x_i, elm_y, star, 
          recStarPolsSingVGen, i;

    # get variable for polynomials
    n := malcevObject!.dim;
    vars_x := GUARANA.RationalVariableList( n, "x" ); 
    vars_y := GUARANA.RationalVariableList( n, "y" );  

    # compute polynomials
    star_pols := [];
    for i in [1..n] do 
	x_i := MalcevSymbolicLieElementByWord( malcevObject, 
	                                       [ [i], vars_x{[i]} ] );
	elm_y := MalcevSymbolicLieElementByWord( malcevObject, 
	                                       [ [i..n], vars_y{[i..n]} ] );
        star := BCHStar( x_i, elm_y );
	star_pols[i] := Coefficients( star );
    od;

    recStarPolsSingVGen := rec( vars_x := vars_x,
                        vars_y := vars_y,
                        pols := star_pols );
    malcevObject!.recStarPolsSingVGen := recStarPolsSingVGen; 
    return 0;                  
end;

GUARANA.MO_AddStarPolynomials := function( malcevObject )
    local n, vars_x, vars_y, star_pols, x, y, star, recStarPols, i;

    # get variable for polynomials
    n := malcevObject!.dim;
    vars_x := GUARANA.RationalVariableList( n, "x" ); 
    vars_y := GUARANA.RationalVariableList( n, "y" );  

    # compute polynomials
    x := MalcevSymbolicLieElementByCoefficients( malcevObject, vars_x );
    y := MalcevSymbolicLieElementByCoefficients( malcevObject, vars_y );
    star := BCHStar( x, y );
    star_pols := Coefficients( star );

    recStarPols := rec( vars_x := vars_x,
                        vars_y := vars_y,
                        pols := star_pols );
    malcevObject!.recStarPols := recStarPols; 
    return 0;                  
end;

#############################################################################
##
#F GUARANA.MO_Star_Symbolic_SingleVersusGeneric(  x, y )
##
##
## IN
## x ............  Malcev lie element with word [ [i], [x_i] ]
## y ...........   Malcev lie elmeent with wore [ [i,...,n], [ y_i,...,y_n]] 
##
## OUT 
## Star  symbolically computed, for the input
## x*y = x_i * \sum_{j=i}^n \alpha_j y_y.
## For symobolic Log and Exp these are the kind of star operation which will 
## be needed. 
##
GUARANA.MO_Star_Symbolic_SingleVersusGeneric := function( x, y  )
    local malcevObject, word_x, word_y, n, index_x, vars_x, vars_y, 
          indets, vals, pols, coeff_res, r, i, coeffs_y;

    # some simple tests
    malcevObject := x!.malcevObject;
    if not IsBound( malcevObject!.recStarPolsSingVGen ) then 
	Error( "recStarPolsSingVGen has to be computed first" );
    fi;
    word_x := x!.word;
    word_y := y!.word;
    coeffs_y := y!.coefficients;
    if Length( word_x[1] ) <> 1 then
	Error( "x has not the correct form" );
    fi;
    if word_x[1][1] > word_y[1][1] then 
        Error( "wrong start index of y\n" ); 
    fi;

    # setup
    n := malcevObject!.dim;
    index_x := word_x[1][1];
    vars_x := malcevObject!.recStarPolsSingVGen.vars_x;
    vars_y := malcevObject!.recStarPolsSingVGen.vars_y;
    indets := Concatenation( vars_x{word_x[1]}, vars_y{[index_x..n]} );
    vals := Concatenation( word_x[2], coeffs_y{[index_x..n]} );

    # get polynomials which are going to be used
    pols := malcevObject!.recStarPolsSingVGen.pols[index_x];

    # compute result
    coeff_res := [ ];
    for i in [1..n] do 
	if i < index_x  then 
	    Add( coeff_res, 0 );
	else 
            r := Value( pols[i], indets, vals );
            Add( coeff_res, r );
	fi;
    od;
    return MalcevLieElementByCoefficients( malcevObject, coeff_res ); 
end;

GUARANA.MO_GetOneOfUnderlyingFieldOfPol := function( pol )
    local ext, coeff;

    ext := ExtRepPolynomialRatFun( pol );
    coeff := ext[2];
    return coeff^0;
end;

GUARANA.MO_Star_Symbolic := function( x, y  )
    local malcevObject, coeffs_x, coeffs_y, n, vars_x, vars_y, 
          indets, vals, pols, coeff_res, one, r, i;

    # some simple tests
    malcevObject := x!.malcevObject;
    if not IsBound( malcevObject!.recStarPols ) then 
	Error( "recStarPols has to be computed first" );
    fi;
    coeffs_x := x!.coefficients;
    coeffs_y := y!.coefficients;

    # setup
    n := malcevObject!.dim;
    vars_x := malcevObject!.recStarPols.vars_x;
    vars_y := malcevObject!.recStarPols.vars_y;
    indets := Concatenation( vars_x, vars_y);
    vals := Concatenation( coeffs_x, coeffs_y );

    # get polynomials which are going to be used
    pols := malcevObject!.recStarPols.pols;

    # compute result
    coeff_res := [ ];
    # check global flag to see whether the coefficients could polynomials 
    # over a number field(which is the case when we compute duSautoy functions)
    if GUARANA.COMP_OVER_EXT_FIELD then
        one := One( GUARANA.EXT_FIELD  );
        for i in [1..n] do
            r := Value( pols[i], indets, vals, 1, one );
            Add( coeff_res, r );
        od;
     else
        for i in [1..n] do
            r := Value( pols[i], indets, vals );
            Add( coeff_res, r );
        od;
    fi;

    if IsSymbolicElement( x ) or IsSymbolicElement( y ) then 
	    return MalcevSymbolicLieElementByCoefficients( malcevObject,coeff_res );
    else
        return MalcevLieElementByCoefficients( malcevObject, coeff_res );
    fi;
end;

#############################################################################
##
#F GUARANA.MO_AddLogPolynomialsByStarPols( malcevObject )
##
## IN
## malcevObject
##
## EFFECT 
## computes the polynomials describing log and stores them in the 
## malcevObject. 
## For this purpose the already computed star polynomials are used.    
## 
GUARANA.MO_AddLogPolynomialsByStarPols:=function( malcevObject )
    local n, vars_e, tail, x_i, i;

    # get variable for polynomials
    n := malcevObject!.dim;
    vars_e := GUARANA.RationalVariableList( n, "e" ); 

    # catch trivial case 
    if n = 0 then 
        malcevObject!.recLogPols := rec( pols := [],
                                         vars_e := [] );
        return 0;
    fi;

    # compute recursively polynomials as follows
    # log( g_1^e_1...g_n^e_n ) = log( g_1^e_1 )*(log( g_2^e_2...g_n^e_n ))
    #                          = e_1 log g_1 *  tail  
    tail := MalcevSymbolicLieElementByWord( malcevObject, 
	                                    [ [n], vars_e{[n]} ] );
    for i in Reversed( [1..(n-1)] ) do 
	x_i := MalcevSymbolicLieElementByWord( malcevObject, 
	                                       [ [i], vars_e{[i]} ] );
        tail := GUARANA.MO_Star_Symbolic_SingleVersusGeneric( x_i, tail );
    od;
    malcevObject!.recLogPols := rec( pols := Coefficients( tail ),
                                     vars_e := vars_e );
    #SetLogMethod( malcevObject, "pols" );
    return 0;
end;

GUARANA.LogByPols := function( g )
    local exp_g, malcevObject, n, indets, pols, coeffs, one, coeff, i;
    
    # setup    
    exp_g := Exponents( g );
    malcevObject := g!.malcevObject;
    n := malcevObject!.dim;
    indets := malcevObject!.recLogPols.vars_e;
    pols := malcevObject!.recLogPols.pols;
    
    # compute coeffs of result
    coeffs := []; 
    # check global flag to see whether the coefficients could polynomials 
    # over a number field(which is the case when we compute duSautoy functions)
    if GUARANA.COMP_OVER_EXT_FIELD then
        one := One( GUARANA.EXT_FIELD );
        for i in [1..n] do
            coeff := Value( pols[i], indets, exp_g, 1, one );
            Add( coeffs, coeff );
        od;
    else
        for i in [1..n] do
            coeff := Value( pols[i], indets, exp_g );
            Add( coeffs, coeff );
        od;
    fi;

    if IsSymbolicElement( g ) then 
	    return MalcevSymbolicLieElementByCoefficients( malcevObject,coeffs );
    else
        return MalcevLieElementByCoefficients( malcevObject, coeffs );
    fi;
end;

#############################################################################
##
#F GUARANA.MO_AddExpPolynomialsByStarPols( malcevObject )
##
## IN
## malcevObject
## 
## OUT
## The polynomials describing Exp are computed and stored in the 
## malcev object. Note that the star polynomials are used
## for this purpose.
## 
GUARANA.MO_AddExpPolynomialsByStarPols := function( malcevObject )
    local n, vars_a, exp_pols, tail, word_tail, a_bar, divider, i;

    # get variable for polynomials
    n := malcevObject!.dim;
    vars_a := GUARANA.RationalVariableList( n, "a" ); 

    # compute recursively polynomials as follows
    # exp( a_1 log g_1 + ... + a_n log g_n ) 
    # = exp( a_1 log g_1 ) * exp( -a_1 log_1 ) *
    #                               exp( a_1 log g_1 + ... + a_n log g_n )
    # = g_1^a_1 * exp( -a_1 log g_1 ) * exp( a_1 log g_1 + ... + a_n log g_n )
    # and then recurse
    exp_pols := [];
    tail := MalcevSymbolicLieElementByCoefficients( malcevObject, 
	                                            vars_a );
    for i in [1..n]  do 
        # get divider 
        word_tail := tail!.word;
        a_bar := word_tail[2][1];
	divider := MalcevSymbolicLieElementByWord( malcevObject, 
	                                    [ [i], [-1* a_bar ]] );
	tail := GUARANA.MO_Star_Symbolic_SingleVersusGeneric( divider, tail);
        Add( exp_pols, a_bar );
    od;
    malcevObject!.recExpPols := rec( pols := exp_pols , vars_a := vars_a );
    #SetExpMethod( malcevObject, "pols" );
    return 0;
end;

GUARANA.ExpByPols := function( x )
    local coeffs_x, malcevObject, indets, pols, coeffs, coeff, i,n;
    
    # setup    
    coeffs_x := Coefficients( x );
    malcevObject := x!.malcevObject;
    n := malcevObject!.dim;
    indets := malcevObject!.recExpPols.vars_a;
    pols := malcevObject!.recExpPols.pols;
    
    # compute coeffs of result
    coeffs := []; 
    for i in [1..n] do
        coeff := Value( pols[i], indets, coeffs_x );
        Add( coeffs, coeff );
    od;
    if IsSymbolicElement( x ) then 
	return MalcevSymbolicGrpElementByExponents( malcevObject,coeffs );
    else
        return MalcevGrpElementByExponents( malcevObject, coeffs );
    fi;
end;

GUARANA.MO_AddLogAndExpPols := function( malcevObject )
    GUARANA.MO_AddStarPolynomialsSingleVersusGeneric( malcevObject );
    GUARANA.MO_AddLogPolynomialsByStarPols( malcevObject );
    GUARANA.MO_AddExpPolynomialsByStarPols( malcevObject );
end;

InstallGlobalFunction( AddLogAndExpPolynomials,
function( malcevObject )
    GUARANA.MO_AddLogAndExpPols( malcevObject );
end);

InstallGlobalFunction( AddStarPolynomials,
function( malcevObject )
    GUARANA.MO_AddStarPolynomials( malcevObject );
end);

InstallGlobalFunction( AddDTPolynomials,
function( malcevObject )
    local T, coll;

    Info( InfoGuarana, 1, "Adding DT polynomials\n" );
    T := malcevObject!.recTGroup.T;
    coll := Collector( T );
    GUARANA.IsWeightedCollector( coll );
    AddHallPolynomials( coll );
end);


if false then 
    malObjs := GUARANA.Get_FNG_MalcevObjects( 2, 9 );
    GUARANA.MO_AddLogAndExpPols( malObjs[7] );

    malObj := malObjs[3];
    x := MalcevLieElementByCoefficients( malObj, [1,2,3,4,5] );
    y := MalcevLieElementByCoefficients( malObj, [1,2,3,4,5] );
    LieBracket( x, y );

    g := MalcevGrpElementByExponents( malObj, [1,2,3,4,5] );
    h := MalcevGrpElementByExponents( malObj, [1,2,3,4,5] );

    n := malObj!.dim;
    vars_a := GUARANA.RationalVariableList( n, "a" );
    vars_b := GUARANA.RationalVariableList( n, "b" );
    z := MalcevSymbolicLieElementByCoefficients( malObj, vars_a );
    o := MalcevSymbolicLieElementByCoefficients( malObj, vars_b );

    i_x := 2;
    x := GUARANA.MalcevBasisLieElement( malObj, i_x, vars_a[i_x] );
    y := MalcevLieElementByWord( malObj, [[i_x..n], vars_b{[i_x..n]}] );


    gg := MalcevSymbolicGrpElementByExponents( malObj, vars_x );

fi;

#############################################################################
##
#E
