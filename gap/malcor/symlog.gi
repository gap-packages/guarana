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
InstallMethod( SetLogMethod, 
               "for Malcev objects and strings (Guarana)", 
               true, 
               [IsMalcevObjectRep, IsString ], 
               0,
function( malcevObject, s)
    local possible_methods;
    possible_methods := [ "pols", "simple" ];
    if s in possible_methods then 
	malcevObject!.log_method := s;
    else
	Error( "Wrong Log method specified" );
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

InstallMethod( SetExpMethod, 
               "for Malcev objects and strings (Guarana)", 
               true, 
               [IsMalcevObjectRep, IsString ], 
               0,
function( malcevObject, s)
    local possible_methods;
    possible_methods := [ "pols", "simple" ];
    if s in possible_methods then 
	malcevObject!.exp_method := s;
    else
	Error( "Wrong Exp method specified" );
    fi;
end );

InstallMethod( ExpMethod, 
               "for Malcev objects (Guarana)", 
               true, 
               [IsMalcevObjectRep ], 
               0,
function( malcevObject)
    return malcevObject!.exp_method;
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
#F GUARANA.MO_AddStarPolynomials( malcevObject )
##
## EFFECT
## recStarPols is added to the Malcev object
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
GUARANA.MO_AddStarPolynomials := function( malcevObject )
    local n, vars_x, vars_y, star_pols, x_i, elm_y, star, recStarPols, i;

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

    recStarPols := rec( vars_x := vars_x,
                        vars_y := vars_y,
                        pols := star_pols );
    malcevObject!.recStarPols := recStarPols; 
    SetLogMethod( malcevObject, "pols" );
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
    if not IsBound( malcevObject!.recStarPols ) then 
	Error( "recStarPols has to be computed first" );
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
    vars_x := malcevObject!.recStarPols.vars_x;
    vars_y := malcevObject!.recStarPols.vars_y;
    indets := Concatenation( vars_x{word_x[1]}, vars_y{[index_x..n]} );
    vals := Concatenation( word_x[2], coeffs_y{[index_x..n]} );

    # get polynomials which are going to be used
    pols := malcevObject!.recStarPols.pols[index_x];

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
    return 0;
end;

GUARANA.LogByPols := function( g )
    local exp_g, malcevObject, indets, pols, coeffs, coeff, i,n;
    
    # setup    
    exp_g := Exponents( g );
    malcevObject := g!.malcevObject;
    n := malcevObject!.dim;
    indets := malcevObject!.recLogPols.vars_e;
    pols := malcevObject!.recLogPols.pols;
    
    # compute coeffs of result
    coeffs := []; 
    for i in [1..n] do
        coeff := Value( pols[i], indets, exp_g );
        Add( coeffs, coeff );
    od;
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
    SetExpMethod( malcevObject, "pols" );
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

if false then 
    malObjs := GUARANA.Get_FNG_MalcevObjects( 2, 4 );
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
## Old code
##


if false then
    malObjs := GUARANA.Get_FNG_MalcevObjects( 2, 4 );
    malObj := malObjs[3];
    dim := malObj!.dim;
    vars_x := GUARANA.RationalVariableList( dim, "x" ); 
    vars_y := GUARANA.RationalVariableList( dim, "y" );  
    x := MalcevLieElementByCoefficients( malObj, vars_x );
    y := MalcevLieElementByCoefficients( malObj, vars_y );
    
fi;
#############################################################################
##
#F GUARANA.GenericElement( n, vars )
##
## IN
## n .............. number of variables
## vars ........... variables that are used. 
##
## OUT
## A generic elment given as list l where 
## l[1] contains the indezes of the elment, i.e. [1..n] and
## l[2] contains the corresponding variables. 
## 
GUARANA.GenericElement := function( n, vars )
    return [[1..n], vars{[1..n]}];
end;

#############################################################################
##
#F GUARANA.GenericElementByDomain( vars, domain )
##
## IN 
## vars ................. list containg all used variables
## domain ..............  a list indezes 
##
## OUT 
## A symbolic element is always given by list containing the used
## indezes and a list containg the corresponding variables. 
##
GUARANA.GenericElementByDomain := function( vars, domain )
    return [ domain, vars{domain}];
end;

GUARANA.GenericBasisElement := function( i, vars )
    return [[i],[vars[i]]];    
end;

#############################################################################
##
#F GUARANA.Sum_Symbolic( list_x, list_y )
##
## IN
## list_x .... a list describing the lie algebra element x.
##              alpha_i log(n_i) is represented as [ [i],[alpha_i] ]
##              alpha_i log(n_i) + alpha_j log(n_j) is represented as
##              [ [i,j],[alpha_i,alpha_j] ] and so on 
## 
## OUT
## Sum of x,y given as symbolic element
##
GUARANA.Sum_Symbolic_new := function( list_x, list_y )
    local length_x,length_y,i_x,i_y,res;
    # catch trivial cases
    if Length( list_x[1] ) = 0 then
        return list_y;
    elif Length( list_y[1] ) = 0 then
        return list_x;
    fi;

    length_x := Length( list_x[1] );
    length_y := Length( list_y[1] );
    i_x := 1;
    i_y := 1;
    res := [[],[]];

    # go through list_x and list_y by using the counters i_x,i_y 
    # and add content to the result
    while i_x <= length_x +1 or i_y <= length_y +1 do
	if i_x = length_x +1 then 
	    # add remaining indezes of list_y
	    Append( res[1], list_y[1]{[i_y..length_y]} );
	    # add remaining content of list_y
	    Append( res[2], list_y[2]{[i_y..length_y]} );
            # return result
	    return res;
	elif i_y = length_y + 1 then
	    # add remaining indezes of list_x
	    Append( res[1], list_x[1]{[i_x..length_x]} );
	    # add remaining content of list_y
	    Append( res[2], list_x[2]{[i_x..length_x]} );
            # return result
	    return res;
	elif list_x[1][i_x] = list_y[1][i_y] then
	    # add index 
	    Add( res[1], list_x[1][i_x] );
	    # add content
	    Add( res[2], list_x[2][i_x] + list_y[2][i_y] );
	    # update counter
	    i_x := i_x+1;
	    i_y := i_y+1;
	elif list_x[1][i_x] < list_y[1][i_y] then
	    # add index
	    Add( res[1], list_x[1][i_x] );
	    # add content
	    Add( res[2], list_x[2][i_x] );
	    # update counter
	    i_x := i_x+1;
	else 
	    # add index
	    Add( res[1], list_y[1][i_y] );
	    # add content
	    Add( res[2], list_y[2][i_y] );
	    # update counter
	    i_y := i_y+1;
	fi;
    od;
    return res;
end;

GUARANA.Sum_Symbolic_old := function( list_x, list_y )
    local res,i,pos;
    # catch trivial cases
    if Length( list_x[1] ) = 0 then
       return list_y;
    elif Length( list_y[1] ) = 0 then
       return list_x;
    fi;

    res := [];
    res[1] := Union( list_x[1], list_y[1] );
    res[2] := List( [1..Length( res[1] )], x->0 );
    for i in [1..Length( list_x[1] )] do
	pos := Position( res[1], list_x[1][i] );
	res[2][pos] := list_x[2][i];
    od;
    for i in [1..Length( list_y[1] )] do
	pos := Position( res[1], list_y[1][i] );
	res[2][pos] := res[2][pos] + list_y[2][i];
    od;
    return res;
end;

GUARANA.Sum_Symbolic := function( list_x, list_y )
    local method;
    #method := "old";
    method := "new";
    if method = "old" then 
	return GUARANA.Sum_Symbolic_old( list_x, list_y );
    else 
	return GUARANA.Sum_Symbolic_new( list_x, list_y );
    fi;
end;

#############################################################################
##
#F GUARANA.EvaluateLieBracket_Symbolic( list_x, list_y, scTable )
## 
## IN
## list_x .... a list describing the lie algebra element x
##             alpha_i log(n_i) is represented as [ [i],[alpha_i] ]
##             alpha_i log(n_i) + alpha_j log(n_j) is represented as
##             [ [i,j],[alpha_i,alpha_j] ]
## scTable ....structure constant table of a Lie algebra.
## 
## OUT
## The lie bracket [x,y] evaluated symbolically
## 
GUARANA.EvaluateLieBracket_Symbolic_recursive := 
                                         function( list_x, list_y, scTable )
    local index_x,index_y,coeff_x,coeff_y,prod,list_x1,list_x2,
          list_y1,list_y2,sum1,sum2,l;
    
    # catch trivial case 
    if Length( list_x[1] )=0 then
        return [[],[]];
    fi;
    if Length( list_y[1] )=0 then
        return [[],[]];
    fi;

    # evaluate lie bracket recursively
    if Length( list_x[1]  ) = 1 and Length( list_y[1] )=1 then 
        index_x := list_x[1][1];
        index_y := list_y[1][1];
        coeff_x := list_x[2][1];
        coeff_y := list_y[2][1];
       
        prod := GUARANA.CopyVectorList( scTable[index_x][index_y] );
        prod[2] := coeff_x*coeff_y*prod[2];
        return prod;
    elif Length( list_x[1] ) > 1 then 
        l := Length( list_x[1] );
        # split x
        list_x1 := [ list_x[1]{[1]} , list_x[2]{[1]} ];
        list_x2 := [ list_x[1]{[2..l]} , list_x[2]{[2..l]} ];
        
        sum1 := GUARANA.EvaluateLieBracket_Symbolic_recursive( list_x1, list_y, scTable );
        sum2 := GUARANA.EvaluateLieBracket_Symbolic_recursive( list_x2, list_y, scTable );
        return GUARANA.Sum_Symbolic( sum1, sum2 );
    elif Length( list_y[1] ) > 1 then
        l := Length( list_y[1] );
        # split y
        list_y1 := [ list_y[1]{[1]} , list_y[2]{[1]} ];
        list_y2 := [ list_y[1]{[2..l]} , list_y[2]{[2..l]} ];
        
        sum1 := GUARANA.EvaluateLieBracket_Symbolic_recursive( list_x, list_y1, scTable );
        sum2 := GUARANA.EvaluateLieBracket_Symbolic_recursive( list_x, list_y2, scTable );
        return GUARANA.Sum_Symbolic( sum1, sum2 );
    fi; 
    Error( "wrong input" );
end;

GUARANA.EvaluateLieBracket_Symbolic_iter := function( list_x, list_y, scTable )
    local res,length_x,length_y,i_x,i_y,index_x,index_y,coeff_x,coeff_y,prod,
          dim,vec,i;
    # catch trivial case 
    if Length( list_x[1] )=0 then
        return [[],[]];
    fi;
    if Length( list_y[1] )=0 then
        return [[],[]];
    fi;

    # setup 
    length_x := Length( list_x[1] );
    length_y := Length( list_y[1] );
    # construct long vector that will be used for summing up
    dim := Length( scTable ) -2;
    vec := List( [1..dim], x-> 0 );

    for i_x in [1..length_x] do
	for i_y in [1..length_y] do
            index_x := list_x[1][i_x];
            index_y := list_y[1][i_y];
	    if index_x <> index_y then
                coeff_x := list_x[2][i_x];
	        coeff_y := list_y[2][i_y];
		prod := GUARANA.CopyVectorList( scTable[index_x][index_y] );
		prod[2] := coeff_x*coeff_y*prod[2];
		# sum up 
		for i in [1..Length(prod[1])] do
		    vec[prod[1][i]] := vec[prod[1][i]] + prod[2][i];
		od;
	     fi;
	od; 
    od;
    
    # transfrom vec to [indezes,pols] form
    res := [[],[]];
    for i in [1..dim] do 
        if vec[i] <> 0 then 
	    Add( res[1], i );
	    Add( res[2], vec[i] );
	fi; 
    od;
    return res;
end;

GUARANA.EvaluateLieBracket_Symbolic := function( list_x,list_y,scTable )
    local method;
    # choose method that should be used
    method := "iter";
    #method := "recursive";
    if method = "iter" then 
	return GUARANA.EvaluateLieBracket_Symbolic_iter( list_x,list_y,scTable);
    else
	return GUARANA.EvaluateLieBracket_Symbolic_recursive(list_x,
	          					     list_y,
							     scTable );
    fi;
end;
	    
#############################################################################
##
#F GUARANA.EvaluateLongLieBracket_Symbolic( list_x, list_y, com, scTable )
##
## IN 
## list_x,list_y ........... symbolic elements represented as lists
## com ..................... commutator, for example [1,2,1] 
## recLieAlg ............... lie algebra record
##
## OUT
## com( x,y) evaluated symbolically.
##
GUARANA.EvaluateLongLieBracket_Symbolic 
:= function( list_x, list_y, com, recLieAlg )
    local r,l,tmp,i;
##      Print( "com: ", com, "\n" );
##      Print( "x :", list_x, "\n" );
##      Print( "y :", list_y, "\n" );
    tmp := [list_x,list_y];
    r := tmp[com[1]];

    l := Length( com );
    for i in [2..l] do
            r := GUARANA.EvaluateLieBracket_Symbolic( r, tmp[com[i]], 
						      recLieAlg.scTable );
    od;
##      Print( "r :", r, "\n" );
    return r;
end;


## The following functions are probably not good enought to be used.
## 
## start

## IN
## recLieAlg ................ lie algebra record
##
## OUT 
## A record containing the information to compute [u,v] for 
## two generic elements u,v.
##
GUARANA.LieBracketPols := function( recLieAlg )
  local n, vars_u, vars_v, u, v, com;

    # get new variables
    n := recLieAlg.dim;
    vars_u := GUARANA.RationalVariableList( n, "u" ); 
    vars_v := GUARANA.RationalVariableList( n, "v" );  

    # get u and v
    u := GUARANA.GenericElement( n, vars_u );
    v := GUARANA.GenericElement( n, vars_v );

    # compute symbolic Lie bracket [u,v]
    com := GUARANA.EvaluateLieBracket_Symbolic( u,v, recLieAlg.scTable ); 

    return rec( vars_u := vars_u,
                vars_v := vars_v, 
		com := com );
end;

##
## EFFECT
## Add the record describing a Lie bracket between two generic symbolic 
## elements to the Lie algebra record.
##
GUARANA.AddLieBracketPols := function( recLieAlg )
  local recLieBracketPols;
    recLieBracketPols := GUARANA.LieBracketPols( recLieAlg );
    recLieAlg.recLieBracketPols := recLieBracketPols;
end;

## IN
##
## OUT
## [x,y] computed symbolically
##
## TODO
## This function does not work correctly if the input 
## has not the full range
GUARANA.LieBracket_Symbolic_GenericVersGeneric
:= function( x, y, recLieAlg )
  local n, vars_u, vars_v, indets, vals, pols, range_result, result, res, i;

    # setup
    n := recLieAlg.dim;
    vars_u := recLieAlg.recLieBracketPols.vars_u;
    vars_v := recLieAlg.recLieBracketPols.vars_v;
    indets := Concatenation( vars_u{x[1]}, vars_v{y[1]} );
    vals := Concatenation( x[2], y[2] );

    # get polynomials
    pols := recLieAlg.recLieBracketPols.com[2];
    range_result := StructuralCopy( recLieAlg.recLieBracketPols.com[1] );

    # compute result
    result := [];
    for i in [1..Length(range_result)] do
        res := Value( pols[i], indets, vals );
        Add( result, res );
    od;
    return [ range_result, result ];
end;

GUARANA.EvaluateLieBracket_Symbolic_UsingPols := function( x,y,recLieAlg )

    end;
## TODO
## Not correct yet if input is not generic.
##
## This method is a little bit quicker for long lie brackets
## and two generic elements in L(F_2,8)
GUARANA.EvaluateLongLieBracket_Symbolic2
:= function( x,y, com, recLieAlg )
  local tmp, r, l, i;

    tmp := [x,y];
    r := tmp[com[1]];

    l := Length( com );
    for i in [2..l] do
       r := GUARANA.LieBracket_Symbolic_GenericVersGeneric( r, tmp[com[i]], 
       						            recLieAlg );
    od;
##      Print( "r :", r, "\n" );
    return r;
end;
## end

## Creating symbolic objects


#############################################################################
##
#F GUARANA.ComputeStarPolys( list_x, list_y, w_x, w_y, max_weight, recLieAlg )
##
## IN 
## list_x .... a list describing the lie algebra element x
## list_y      alpha_i log(n_i) is represented as [ [i],[alpha_i] ]
##             alpha_i log(n_i) + alpha_j log(n_j) is represented as
##             [ [i,j],[alpha_i,alpha_j] ], etc...
## wx,wy ..... weight of x,y 
## max_weight ..... maximal weight of basis elms of the Lie algebra
## recLieAlg.. structure constant table
##
## OUT 
## x * y evaluated symbolically
## 
GUARANA.ComputeStarPolys 
:= function( list_x, list_y, wx, wy, max_weight, recLieAlg )
    local i,r,bchSers,com,a,term,max,min,bound;
    bchSers := GUARANA.recBCH.bchSers;
    
    # start with terms which are not given by Lie brackets
    r := GUARANA.Sum_Symbolic( list_x, list_y );

    # trivial check 
    if Length( list_x[1] ) = 0  or Length( list_y[1] ) =  0 then
        return r;
    fi;

    # compute upper bound for the Length of commutators, which 
    # can be involved
    max := Maximum( wx,wy );
    min := Minimum( wx,wy );
    # max + min* (bound-1 ) <= max_weight
    bound := Int( (max_weight-max)/min + 1 );

    # up to bound  compute the commutators and add them.
    # Note that the list contains commutators of length i at position i-1.
    for i in [1..bound-1] do
        for term in bchSers[i] do
            com := term[2];
            # check if weight of commutator is not to big
            if GUARANA.CheckWeightOfCommutator( com, wx, wy, max_weight ) then
                a := GUARANA.EvaluateLongLieBracket_Symbolic( 
		     list_x, list_y, com, recLieAlg );
                r := GUARANA.Sum_Symbolic( r,[a[1],term[1]*a[2]] );
                #r := r + term[1]*a; 
            fi;
         od;
    od;
    return r;
end;

#############################################################################
##
#F GUARANA.AddStarPolynomialsToRecLieAlg( recLieAlg )
##
## IN
## recLieAlg .............. lie algebra record
## 
## EFFECT
## recStarPols is added to the Lie algebra record
##
## For Log and Exp I just need 
## log x_i * sum_{j=i}^l beta_j log x_j
## So these are the polynomials I should save. It is NOT necessary to
## save x_i against a generic element of L.
## TODO
## compute polynomials which can used to speed up star operation
## How can I speed up this function ??
##
## Example:
## recLieAlgs_bch_F2c := GUARANA.Get_FNG_LieAlgRecords( 2, 9 );
## GUARANA.AddStarPolynomialsToRecLieAlg( recLieAlgs_bch_F2c[7] );
##
GUARANA.AddStarPolynomialsToRecLieAlg := function( recLieAlg )
    local i,n,vars_x,vars_y,x_i,elm_y,wx,wy,recStarPols,c,star_pols;

    # get variable for polynomials
    n := recLieAlg.dim;
    vars_x := GUARANA.RationalVariableList( n, "x" ); 
    vars_y := GUARANA.RationalVariableList( n, "y" );  

    # compute polynomials
    c := recLieAlg.max_weight;
    star_pols := [];
    for i in [1..n] do 
        x_i := GUARANA.GenericBasisElement( i, vars_x );
        elm_y := GUARANA.GenericElementByDomain( vars_y, [i..n] );
        wx := recLieAlg.weights[i];
        wy := recLieAlg.weights[i];
        star_pols[i] := GUARANA.ComputeStarPolys( 
                  x_i, elm_y,  wx, wy, c, recLieAlg );
    od;

    recStarPols := rec( vars_x := vars_x,
                        vars_y := vars_y,
                        pols := star_pols );
    recLieAlg.recStarPols := recStarPols; 
    return 0;                  
end;

# profiling the computation of star pols.
GUARANA.Profile := function( func, input,  subfuncs )
    local res;
    ProfileFunctions( subfuncs );
    ClearProfile();
    res := func( input );
    DisplayProfile();
    ClearProfile();
    return 0;
end;

if false then
    recLieAlgs_bch_F2c := GUARANA.Get_FNG_LieAlgRecords( 2, 9 );;
    input :=  recLieAlgs_bch_F2c[8];;
    func := GUARANA.AddStarPolynomialsToRecLieAlg;
    subfuncs := [ GUARANA.Sum_Symbolic, 
		  GUARANA.CopyVectorList,
		  GUARANA.EvaluateLieBracket_Symbolic ,
		  GUARANA.EvaluateLongLieBracket_Symbolic,
		  GUARANA.AddStarPolynomialsToRecLieAlg ];
    GUARANA.Profile( func, input, subfuncs );
    
fi;
#############################################################################
##
#F GUARANA.Star_Symbolic_SingleVersusGeneric( recLieAlg, x, y )
##
##
## IN
## recLieAlg ......Lie algebra record
## x ............  [ [i], [x_i] ] 
##                 where x_i in general will be a variable.
##                 It can be also a polynomial, or some other ring element
## y ...........   [ [i,...,n], [ y_i,...,y_n]] 
##
## OUT 
## Star  symbolically computed, for the input
## x*y = x_i * \sum_{j=i}^n \alpha_j y_y.
## For Log and Exp these are the kind of star operation which will 
## be needed. 
##
GUARANA.Star_Symbolic_SingleVersusGeneric := function( recLieAlg, x, y  )
    local n,vars_y,vars_x,indets,vals,pols,result,i,res,index_x,range_result;

    # simple test 
    if not x[1][1] = y[1][1] then Error( "wrong start index of y\n" ); fi;

    # setup
    n := recLieAlg.dim;
    index_x := x[1][1];
    vars_y := recLieAlg.recStarPols.vars_y;
    vars_x := recLieAlg.recStarPols.vars_x;
    indets := Concatenation( vars_x{x[1]}, vars_y{y[1]} );
    vals := Concatenation( x[2], y[2] );

    # get polynomials which are going to be used
    pols := recLieAlg.recStarPols.pols[index_x][2];
    range_result := StructuralCopy( recLieAlg.recStarPols.pols[index_x][1] );

    # compute result
    result := [];
    for i in [1..Length(range_result)] do
        res := Value( pols[i], indets, vals );
        Add( result, res );
    od;
    return [ range_result, result ];
end;

#############################################################################
##
#F GUARANA.ComputeSymbolicLogPolynomials( recLieAlg )
##
## IN
## recLieAlg ............... Lie algebra record
## 
## OUT
## For a given Lie algebra L(N) of dimension l, compute polynomials
## p_1,...,p_l such that 
## Log( g_1^e_1 ... g_l^e_l ) = p_1( e )Log g_1 + ... + p_l( e ) Log g_l 
## where (g_1,...,g_l) is a Malcev basis for N
## 
## Note: This can be speed up via using recStarPols
## 
GUARANA.ComputeSymbolicLogPolynomials := function( recLieAlg )
    local n, vars_e,c,log_pols,tail,x_i,w_x_i,w_tail,i;

    # get variables for polynomials
    n := recLieAlg.dim;
    vars_e := GUARANA.RationalVariableList( n, "e" ); 

    # compute recursively polynomials as follows
    # log( g_1^e_1...g_n^e_n ) = log( g_1^e_1 )*(log( g_2^e_2...g_n^e_n ))
    #                          = e_1 log g_1 *  tail  
    c := recLieAlg.max_weight;
    log_pols := [];
    tail := GUARANA.GenericBasisElement( n, vars_e );
    for i in Reversed( [1..(n-1)] ) do 
        x_i := GUARANA.GenericBasisElement( i, vars_e );
        w_x_i := recLieAlg.weights[i];
        w_tail := recLieAlg.weights[i+1];
        tail := GUARANA.ComputeStarPolys(  x_i, tail,  w_x_i, 
                                      w_tail, c, recLieAlg );
    od;
    return rec( pols := tail, vars_e := vars_e );
end;

# converts  [[i+1..n],[x_{i+1}..x_n]] to 
#           [[i..n],  [0,x_{i+1}..x_n]]
# 
GUARANA.Convert1 := function( x )
    local x_new,i,n;
    x_new := [];
    i := x[1][1] -1;
    n := x[1][Length(x[1])];
    x_new[1] := [i..n];
    x_new[2] := Concatenation( [ 0*x[2][1] ], x[2] );
    return x_new;
end;

#############################################################################
##
#F GUARANA.ComputeSymbolicLogPolynomialsByStarPols( recLieAlg )
##
## IN
## recLieAlg ................. Lie algebra record
##
## OUT 
## computes the polynomials describing log and stores them in the 
## lie algebra record. 
## For this purpose the already computed star polynomials are used.    
## 
## TODO
## Why is this so quick compared to the computation of the star pols ?
## 
GUARANA.ComputeSymbolicLogPolynomialsByStarPols := function( recLieAlg )
    local n,vars_e,c,log_pols,tail,x_i,i;

    # get variable for polynomials
    n := recLieAlg.dim;
    vars_e := GUARANA.RationalVariableList( n, "e" ); 

    # compute recursively polynomials as follows
    # log( g_1^e_1...g_n^e_n ) = log( g_1^e_1 )*(log( g_2^e_2...g_n^e_n ))
    #                          = e_1 log g_1 *  tail  
    c := recLieAlg.max_weight;
    tail := GUARANA.GenericBasisElement( n, vars_e );
    for i in Reversed( [1..(n-1)] ) do 
        x_i := GUARANA.GenericBasisElement( i, vars_e );
        tail := GUARANA.Convert1( tail );
        tail := GUARANA.Star_Symbolic_SingleVersusGeneric(recLieAlg,x_i,tail);
    od;
    return rec( pols := tail, vars_e := vars_e );
end;

GUARANA.AddLogPolynomialsToLieAlgRecord := function( recLieAlg )
    local pols, recLogPols;
    recLogPols := GUARANA.ComputeSymbolicLogPolynomialsByStarPols( 
                                                      recLieAlg );
    recLieAlg.recLogPols := recLogPols;
    recLieAlg.log_method := "symbolic";
    return 0;
end;

#############################################################################
##
#F GUARANA.Logarithm_Symbolic( args )
##
## IN
## args[1]=recLieAlg ...... Lie algebra record
## args[2[=exp_n .......... exponent vector of element n in \hat{N}
##    
## OUT: 
## coeffcients of Log n, computed with the polynomials describing Log
## 
##
##  Example:
##  recLieAlgs_bch_F2c := GUARANA.Get_FNG_LieAlgRecords( 2, 9 );
##  GUARANA.AddStarPolynomialsToRecLieAlg( recLieAlgs_bch_F2c[5] );
##  GUARANA.AddLogPolynomialsToLieAlgRecord( recLieAlgs_bch_F2c[5] );
##  exp_n := GUARANA.Random_IntegralExpVector( recLieAlgs_bch_F2c[5], 2^12 );  
##  GUARANA.Logarithm_Symbolic( [recLieAlgs_bch_F2c[5], exp_n] );
##
##  compare 
##  GUARANA.AbstractLog_Simple_ByExponent( recLieAlgs_bch_F2c[5], exp_n);
## 
GUARANA.Logarithm_Symbolic := function( args )
    local n,indets,pols,coeffs,i,coeff,recLieAlg,exp_n;
    # setup    
    recLieAlg := args[1];
    exp_n := args[2];
    n := recLieAlg.dim; 
    indets := recLieAlg.recLogPols.vars_e;
    pols := recLieAlg.recLogPols.pols[2];
    
    # compute coeffs of result
    coeffs := []; 
    for i in [1..n] do
        coeff := Value( pols[i], indets, exp_n );
        Add( coeffs, coeff );
    od;
    return coeffs;
end;

#############################################################################
##
#F GUARANA.ComputeSymbolicExpPolynomialsByStarPols( recLieAlg )
##
## IN
## recLieAlg ................. lie algebra record
## 
## OUT
## The polynomials describing Exp are computed and stored in the 
## Lie algebra record. Note that the star polynomials are used
## for this purpose.
## 
GUARANA.ComputeSymbolicExpPolynomialsByStarPols := function( recLieAlg )
    local i,n,vars_a,c,exp_pols,tail,a_bar,divider;

    # get variable for polynomials
    n := recLieAlg.dim;
    vars_a := GUARANA.RationalVariableList( n, "a" ); 

    # compute recursively polynomials as follows
    # exp( a_1 log g_1 + ... + a_n log g_n ) 
    # = exp( a_1 log g_1 ) * exp( -a_1 log_1 ) *
    #                               exp( a_1 log g_1 + ... + a_n log g_n )
    # = g_1^a_1 * exp( -a_1 log g_1 ) * exp( a_1 log g_1 + ... + a_n log g_n )
    # and then recurse
    c := recLieAlg.max_weight;
    exp_pols := [];
    tail := GUARANA.GenericElement( n , vars_a  );
    for i in [1..n]  do 
        # get divider 
        a_bar := tail[2][1];
        divider := [ [i], [ (-1)* a_bar ] ];
        tail:=GUARANA.Star_Symbolic_SingleVersusGeneric(recLieAlg,divider,tail);
        Remove( tail[1], 1 );
        Remove( tail[2], 1 );
        Add( exp_pols, a_bar );
    od;
    return rec( pols := [[1..n],exp_pols], vars_a := vars_a );
end;

GUARANA.AddExpPolynomialsToLieAlgRecord := function( recLieAlg )
    local pols, recExpPols;
    recExpPols := GUARANA.ComputeSymbolicExpPolynomialsByStarPols( recLieAlg );
    recLieAlg.recExpPols := recExpPols;

    recLieAlg.exp_method := "symbolic";
    return 0;
end;

#############################################################################
##
#F GUARANA.Exponential_Symbolic( args )
##
## IN
## args[1]=recLieAlg .............. lie algebra record
## args[2]=coeffs_x ............... coefficient vector of Lie algebra elm. x
##
## OUT 
## Exp(x) computed using polynomials.
## 
## Example:
## recLieAlgs_bch_F2c := GUARANA.Get_FNG_LieAlgRecords( 2, 9 );
## GUARANA.AddStarLogAndExpPols( recLieAlgs_bch_F2c[7] );
## coeffs_x :=GUARANA.Random_IntegralExpVector( recLieAlgs_bch_F2c[7], 2^12 );  
## GUARANA.Exponential_Symbolic( [recLieAlgs_bch_F2c[7], coeffs_x] );
##  
## compare 
## GUARANA.Abstract_Exponential_ByVector( recLieAlgs_bch_F2c[7], coeffs_x );
##
GUARANA.Exponential_Symbolic := function( args )
    local recLieAlg,n,indets,pols,exp,i,e,coeffs_x;
    
    # setup
    recLieAlg := args[1];
    coeffs_x := args[2];
    n := recLieAlg.dim; 
    indets := recLieAlg.recExpPols.vars_a;
    pols := recLieAlg.recExpPols.pols[2];
    
    # compute exponents of result
    exp := []; 
    for i in [1..n] do
        e := Value( pols[i], indets, coeffs_x );
        Add( exp , e );
    od;
    return exp;
end;

#############################################################################
##
#F GUARANA.AddStarLogAndExpPols( recLieAlg )
##
## IN
## recLieAlg ..................  Lie algebra record
## 
## EFFECT
## Star, Log and Exp polynomials are added to the Lie algebra record.
##
## Example usage:
## recLieAlgs_bch_F2c := GUARANA.Get_FNG_LieAlgRecords( 2, 9 );
## GUARANA.AddStarLogAndExpPols( recLieAlgs_bch_F2c[7] );
##
GUARANA.AddStarLogAndExpPols := function( recLieAlg )
    GUARANA.AddStarPolynomialsToRecLieAlg( recLieAlg ); 
    GUARANA.AddLogPolynomialsToLieAlgRecord( recLieAlg ); 
    GUARANA.AddExpPolynomialsToLieAlgRecord( recLieAlg ); 
    return 0;
end;

#############################################################################
##
#E
