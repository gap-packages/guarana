#############################################################################
##
#W malelm.gi              GUA_package                     Bjoern Assmann
##
#H  @(#)$Id$
##
#Y 2006
##
##

#############################################################################
##
## Methods for constructing Malcev elements
##
InstallGlobalFunction( MalcevGenElementByExponents, 
function( malcevObject, exps )
    local g, weight, x, elm, lie_elm;

   # get group element
   g := MalcevGrpElementByExponents( malcevObject, exps );

   # get weight 
   weight := Weight( g );

   # lie element unknown so far
   x := "unknown yet";

   elm := rec( malcevObject := malcevObject,
               grp_elm := Immutable( g ),
               lie_elm := x,
               weight := weight );
   return Objectify( malcevObject!.gen_elms_type , elm );
end);

InstallGlobalFunction( MalcevGenElementByCoefficients, 
function( malcevObject, coeffs )
    local g, weight, x, elm, lie_elm;

   # group element unknown so far
   g := "unknown yet";

   # get lie element 
   x := MalcevLieElementByCoefficients( malcevObject, coeffs ); 

   # get weight 
   weight := Weight( x );


   elm := rec( malcevObject := malcevObject,
               grp_elm := g,
               lie_elm := Immutable( x ),
               weight := weight );
   return Objectify( malcevObject!.gen_elms_type , elm );
end);


InstallGlobalFunction( MalcevGenElementByLieElement, 
function(  x )
    local malcevObject, g, weight, elm;

    malcevObject := x!.malcevObject;

    # get group element
    g := "unknown yet";

    # get weight 
    weight := Weight( x );


   elm := rec( malcevObject := malcevObject,
               grp_elm := g,
               lie_elm := Immutable( x ),
               weight := weight );
   return Objectify( malcevObject!.gen_elms_type , elm );
end);

InstallGlobalFunction( MalcevGenElementByGrpElement, 
function(  g )
    local malcevObject, x, weight, elm;

    malcevObject := g!.malcevObject;

    # lie element element
    x := "unknown yet";

    # get weight 
    weight := Weight( g );


   elm := rec( malcevObject := malcevObject,
               grp_elm := Immutable( g ),
               lie_elm := x,
               weight := weight );
   return Objectify( malcevObject!.gen_elms_type , elm );
end);

InstallGlobalFunction( MalcevSymbolicGenElementByExponents, 
function( malcevObject, exps )
    local g, weight, x, elm, lie_elm;

   # get group element
   g := MalcevSymbolicGrpElementByExponents( malcevObject, exps );

   # get weight 
   weight := Weight( g );

   # lie element unknown so far
   x := "unknown yet";

   elm := rec( malcevObject := malcevObject,
               grp_elm := Immutable( g ),
               lie_elm := x,
               weight := weight );
   elm := Objectify( malcevObject!.gen_elms_type , elm );
   Setter( IsSymbolicElement )( elm, true );
   return elm;
end);

InstallGlobalFunction( MalcevSymbolicGenElementByCoefficients, 
function( malcevObject, coeffs )
    local g, weight, x, elm, lie_elm;

   # group element unknown so far
   g := "unknown yet";

   # get lie element 
   x := MalcevSymbolicLieElementByCoefficients( malcevObject, coeffs ); 

   # get weight 
   weight := Weight( x );


   elm := rec( malcevObject := malcevObject,
               grp_elm := g,
               lie_elm := Immutable( x ),
               weight := weight );
   elm :=  Objectify( malcevObject!.gen_elms_type , elm );
   Setter( IsSymbolicElement )( elm, true );
   return elm;
end);

#############################################################################
##
#M Print Malcev gen elements
##
InstallMethod( PrintObj, 
               "for Malcev gen elements (Guarana)", 
               true, 
               [IsMalcevGenElement ], 
               0,
function( elm )
    Print( "Group element: ", elm!.grp_elm, "\n" );
    Print( "Lie element:   ", elm!.lie_elm );
end );

InstallMethod( LieElement,
               "for Malcev Gen element (Guarana)",
               [IsMalcevGenElement],
               0,
function( x )
    if IsString( x!.lie_elm ) then 
        x!.lie_elm := Immutable( Log( x!.grp_elm ) );
    fi;
    return x!.lie_elm;
end);
           
InstallMethod( GrpElement,
               "for Malcev Gen element (Guarana)",
               [IsMalcevGenElement],
               0,
function( x )
    if IsString( x!.grp_elm ) then 
        x!.grp_elm := Immutable( Exp( x!.lie_elm ) );
    fi;
    return x!.grp_elm;
end);

InstallOtherMethod( Exponents, 
               "for Malcev Gen element (Guarana)", 
	       [IsMalcevGenElement ],
	       0,
function( x )
    return Exponents( GrpElement( x ) );
end);

InstallOtherMethod( Coefficients, 
               "for Malcev Gen element (Guarana)", 
	       [IsMalcevGenElement ],
	       0,
function( x )
    return Coefficients( LieElement( x ) );
end);

#############################################################################
##
#M a * b .......................................Product of Malcev Gen elments
##
GUARANA.MultViaStar := function( a, b )
    local x_a, x_b, x_res;
    x_a := LieElement( a );
    x_b := LieElement( b );
    x_res := BCHStar( x_a, x_b );
    return MalcevGenElementByLieElement( x_res );
end;

GUARANA.MultViaCollection := function( a, b )
    local g_a, g_b, g_res;

    g_a := GrpElement( a );
    g_b := GrpElement( b );
    # note that in the group DT collection can be used.
    g_res := g_a * g_b; 
    return MalcevGenElementByGrpElement( g_res );
end;

InstallOtherMethod( \*, 
               "for Malcev Gen elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevGenElement, IsMalcevGenElement ],
		0, 
function( a, b  )
    local malcevObject;
    malcevObject := a!.malcevObject;
    if malcevObject!.mult_method = GUARANA.MultMethodIsStar then 
        return GUARANA.MultViaStar( a, b );
    elif malcevObject!.mult_method = GUARANA.MultMethodIsCollection then
        return GUARANA.MultViaCollection( a, b );
    else
        Error( " " );
    fi;
end);

#############################################################################
##
## Methods for constructing Malcev lie elements.
##
InstallGlobalFunction( MalcevLieElementConstruction, 
function( malcevObject, coefficients, word )
    local weight, i, elm, name;
    
    # determine weight
    if Length( word[1] ) = 0 then 
	weight := malcevObject!.max_weight + 1;
    else 
	i := word[1][1];
	weight := malcevObject!.weights[i];
    fi;

    elm := rec( malcevObject := malcevObject,
	        coefficients := Immutable( coefficients ),
		word := Immutable( word ),
		name := "x",
	        weight := weight );
    return Objectify( malcevObject!.lie_elms_type , elm );
end );

GUARANA.Coefficients2Word := function( malcevObject, coeffs  )
    local n, word, i;

    n := Length( coeffs );
    if n <> malcevObject!.dim then 
	Error( "Length of coefficient vector not correct" ); 
    fi;
    word := [[],[]];
    for i in [1..n] do 
	if coeffs[i] <> 0*coeffs[i] then 
	    Add( word[1], i );
	    Add( word[2], coeffs[i] );
	fi;
    od;
    return word;
end;

InstallGlobalFunction( MalcevLieElementByCoefficients,
function( malcevObject, coeffs )
    local word; 
    # check input
    if Length( coeffs ) <> Dimension( malcevObject ) then 
        Error( "Length of coefficient vector does not match dimension of Mal'cev object.");
    fi;
    word := GUARANA.Coefficients2Word( malcevObject, coeffs );
    return MalcevLieElementConstruction( malcevObject, coeffs, word );
end);

GUARANA.Word2Coefficients := function( malcevObject, word )
    local coeffs, i;
    coeffs := List( [1..malcevObject!.dim], x-> 0 );
    for i in [1..Length( word[1] )] do 
	coeffs[word[1][i]] := word[2][i];
    od;
    return coeffs;
end;

InstallGlobalFunction( MalcevLieElementByWord, 
function( malcevObject, word )
    local coeffs;
    coeffs := GUARANA.Word2Coefficients( malcevObject, word );
    return MalcevLieElementConstruction( malcevObject, coeffs, word );
end);

#############################################################################
##
## IN
## malcevObject
## i ................................. i in [1..dim] where dim
##                                     is the dimension of the Lie algebra.
## coeff ............................. coefficent
##
## OUT
## The Malcev lie elment with trivial coefficent vector except at the
## position i where we have the coefficient coeff.
##
GUARANA.MalcevBasisLieElement := function( malcevObject, i, coeff )
    local dim, coeffs;
    dim := malcevObject!.dim;
    coeffs := List( [1..dim], x-> 0 );
    coeffs[i] := coeff;
    return MalcevLieElementByCoefficients( malcevObject, coeffs );
end;

#############################################################################
##
## Methods for constructing symbolic Malcev Lie elements.
##
MalcevSymbolicLieElementByWord := function( malcevObject, word )
    local elm;
    elm := MalcevLieElementByWord( malcevObject, word );
    Setter( IsSymbolicElement )( elm, true );
    return elm;
end;

InstallGlobalFunction( MalcevSymbolicLieElementByCoefficients,
function( malcevObject, coeffs )
    local elm;
    elm := MalcevLieElementByCoefficients( malcevObject, coeffs );
    Setter( IsSymbolicElement )( elm, true );
    return elm;
end);

InstallMethod( IsSymbolicElement, 
               "for Malcev elements (Guarana)", 
	       true,
	       [IsMalcevElement],
	       0,
function( elm )
    return Tester( IsSymbolicElement )( elm );
end);

# TODO: What does SetFeatureObj( x, IsSymbolicElement, true ) do ? 
#

#############################################################################
##
#M Print Malcev Lie elements
##
InstallMethod( PrintObj, 
               "for Malcev elements (Guarana)", 
               true, 
               [IsMalcevLieElement ], 
               0,
function( elm )
    local coeffs;
    coeffs := elm!.coefficients;
    if Length( coeffs ) = 0 then 
        Print( "<zero of Lie algebra>" );
    else
        Print( coeffs );
    fi;
end );

InstallOtherMethod( Coefficients, 
               "for Malcev Lie element", 
	       [IsMalcevLieElement ],
	       0,
x -> x!.coefficients
);

InstallOtherMethod( Weight, 
               "for Malcev element (Guarana)", 
	       [IsMalcevElement ],
	       0,
x -> x!.weight
);

#############################################################################
##
## Methods for constructing Malcev group elements.
##
InstallGlobalFunction( MalcevGrpElementConstruction, 
function( malcevObject, exponents )
    local weight, i, elm, name;

    # check input
    if Length( exponents ) <> Dimension( malcevObject ) then 
        Error( "Length of exponent vector does not match dimension of Mal'cev object.");
    fi;
    
    # determine weight
    weight := malcevObject!.max_weight + 1;
    for i in [1..Length( exponents )] do 
        if exponents[i] <> 0 then 
            weight := malcevObject!.weights[i];
	    break;
	fi;
    od;

    elm := rec( malcevObject := malcevObject,
	        exponents := Immutable( exponents ),
		name := "g",
	        weight := weight );
    return Objectify( malcevObject!.grp_elms_type , elm );
end );

InstallGlobalFunction( MalcevGrpElementByExponents,
function( malcevObject, coeffs )
    return MalcevGrpElementConstruction( malcevObject, coeffs );
end);

#############################################################################
##
## Methods for constructing symbolic Malcev grp elements.
##
InstallGlobalFunction( MalcevSymbolicGrpElementByExponents,
function( malcevObject, exp )
    local elm;
    elm := MalcevGrpElementByExponents( malcevObject, exp );
    Setter( IsSymbolicElement )( elm, true );
    return elm;
end);

#############################################################################
##
## IN
## malcevObject
## i ................................. i in [1..dim] where dim
##                                     is the dimension of the Lie algebra.
## exp .............................   exponent
##
## OUT
## The Malcev grroup elment with trivial coefficent vector except at the
## position i where we have the exponent exp.
##
GUARANA.MalcevBasisGrpElement := function( malcevObject, i, exp )
    local dim, exps;
    dim := malcevObject!.dim;
    exps := List( [1..dim], x-> 0 );
    exps[i] := exp;
    return MalcevGrpElementByExponents( malcevObject, exps );
end;

#############################################################################
##
#M Print Malcev grp elements
##
InstallMethod( PrintObj, 
               "for Malcev grp elements (Guarana)", 
               true, 
               [IsMalcevGrpElement ], 
               0,
function( elm )
    local exps;
    exps := elm!.exponents;
    if Length( exps ) = 0 then 
        Print( "id" );
    else
        Print( exps );
    fi;
end );

InstallOtherMethod( Exponents, 
               "for Malcev Grp element", 
	       [IsMalcevGrpElement ],
	       0,
x -> x!.exponents
);

#############################################################################
##
#M x + y .......................................... Sum of Malcev Lie elments
##
InstallOtherMethod( \+, 
               "for Malcev Lie elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevLieElement, IsMalcevLieElement ],
		0, 
function( x, y )
    local malObj, coeff_x, coeff_y, coeff_res;

    malObj := x!.malcevObject;
    coeff_x := Coefficients( x );
    coeff_y := Coefficients( y );

    coeff_res := coeff_x + coeff_y;
    if IsSymbolicElement( x ) or IsSymbolicElement( y ) then
	return MalcevSymbolicLieElementByCoefficients( malObj, coeff_res );
    else 
	return MalcevLieElementByCoefficients( malObj, coeff_res );
    fi;
end);

#############################################################################
##
#M k * x ...................... Scalar mulitplication for Malcev Lie elements
##
InstallOtherMethod( \*, 
               "producto of scalar with  Malcev Lie elments (Guarana)",
	       true,
	        [IsScalar, IsMalcevLieElement ],
		0, 
function( k, x )
    local malObj, coeff_x, coeff_res;

    malObj := x!.malcevObject;
    coeff_x := Coefficients( x );

    coeff_res := k * coeff_x;
    if IsSymbolicElement( x ) or IsPolynomial( k ) then
	return MalcevSymbolicLieElementByCoefficients( malObj, coeff_res );
    else 
	return MalcevLieElementByCoefficients( malObj, coeff_res );
    fi;
end);

InstallOtherMethod( AdditiveInverseMutable, 
               "additive inverse for Malcev Lie elments (Guarana)",
	       true,
	        [IsMalcevLieElement ],
		0, 
function(  x )
    local malObj, coeff_x, coeff_res;

    malObj := x!.malcevObject;
    coeff_x := Coefficients( x );

    coeff_res := -coeff_x;
    if IsSymbolicElement( x ) then
	return MalcevSymbolicLieElementByCoefficients( malObj, coeff_res );
    else 
	return MalcevLieElementByCoefficients( malObj, coeff_res );
    fi;
end);

#############################################################################
##
#M x = y 
##
InstallOtherMethod( \=, 
               "for Malcev Lie elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevLieElement, IsMalcevLieElement ],
		0, 
function( x, y )
    return Coefficients( x ) = Coefficients( y );
end);


#############################################################################
##
#M LieBracket( x, y ) ............................... for Malcev lie elements
##
InstallOtherMethod( LieBracket, 
               "for Malcev Lie elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevLieElement, IsMalcevLieElement ],
		0, 
function( x, y )
    local malObj, basis, coeff_x, coeff_y, xx, yy, res, coeff_res;

    # get corresponding elms of lie algebra
    malObj := x!.malcevObject;
    basis := Basis( malObj!.L );
    coeff_x := Coefficients( x );
    coeff_y := Coefficients( y );
    xx := LinearCombination( basis, coeff_x );
    yy := LinearCombination( basis, coeff_y );

    # compute lie bracket 
    if malObj!.lieAlgebraType = "structureConstants" then 
        res := xx*yy;
        coeff_res := Coefficients( basis, res );
    elif malObj!.lieAlgebraType = "matrix" then 
	res := LieBracket( xx, yy );
        coeff_res := Coefficients( basis, res );
    fi;

    return MalcevLieElementByCoefficients( malObj, coeff_res );
end);

InstallOtherMethod( Comm,
               "for Malcev Lie elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevLieElement, IsMalcevLieElement ],
		0, 
function( x, y )
    return LieBracket( x,y );
end);

#############################################################################
##
#M LieBracket( x, y ) ...................... for symbolic Malcev lie elements
##
GUARANA.SymbolicLieBracket := function( x, y )
    local malObj, scTable, word_x, word_y, length_x, length_y, dim, 
          vec, index_x, index_y, coeff_x, coeff_y, prod, i_x, i_y, i;
 
    # catch trivial case 
    word_x := x!.word;
    word_y := y!.word;
    if Length( word_x[1] )=0 then
        return 0*x;
    fi;
    if Length( word_y[1] )=0 then
        return 0*x;
    fi;

    # set up
    malObj := x!.malcevObject;
    scTable := malObj!.scTable;
    length_x := Length( word_x[1] );
    length_y := Length( word_y[1] );

    # get coefficient vector that will be used for summing up
    dim := malObj!.dim;
    vec := List( [1..dim], x-> 0 );

    for i_x in [1..length_x] do
	for i_y in [1..length_y] do
            index_x := word_x[1][i_x];
            index_y := word_y[1][i_y];
	    if index_x <> index_y then
                coeff_x := word_x[2][i_x];
	        coeff_y := word_y[2][i_y];
		prod := GUARANA.CopyVectorList( scTable[index_x][index_y] );
		prod[2] := coeff_x*coeff_y*prod[2];
		# sum up 
		for i in [1..Length(prod[1])] do
		    vec[prod[1][i]] := vec[prod[1][i]] + prod[2][i];
		od;
	     fi;
	od; 
    od;
    return MalcevSymbolicLieElementByCoefficients( malObj, vec );
end;

InstallOtherMethod( LieBracket, 
               "for symbolic Malcev Lie elments (Guarana)",
	       IsIdenticalObj,
	       [IsMalcevLieElement and IsSymbolicElement, IsMalcevLieElement ],
		0, 
function( x, y )
    return GUARANA.SymbolicLieBracket( x,y );
end);

InstallOtherMethod( LieBracket, 
               "for symbolic Malcev Lie elments (Guarana)",
	       IsIdenticalObj,
	       [IsMalcevLieElement, IsMalcevLieElement and IsSymbolicElement],
		0, 
function( x, y )
    return GUARANA.SymbolicLieBracket( x, y );
end);

#############################################################################
##
#M g* h  ................................... Product of Malcev group elements
##
InstallOtherMethod( \*, 
               "for Malcev group elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevGrpElement, IsMalcevGrpElement ],
		0, 
function( g, h )
    local malObj, exp_g, exp_h, coll, gg, hh, exp_res;

    malObj := g!.malcevObject;
    exp_g := Exponents( g );
    exp_h  := Exponents( h );

    # get corresponing pcp elms
    coll := Collector( malObj!.recTGroup.T );
    gg := PcpElementByExponentsNC( coll, exp_g );
    hh := PcpElementByExponentsNC( coll, exp_h );

    exp_res := Exponents( gg*hh );
    if IsSymbolicElement( g ) or IsSymbolicElement( h ) then
	    return MalcevSymbolicGrpElementByExponents( malObj, exp_res );
    else 
	    return MalcevGrpElementByExponents( malObj, exp_res );
    fi;
end);

InstallOtherMethod( \*, 
               "for symbolic Malcev group elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevGrpElement and IsSymbolicElement, IsMalcevGrpElement ],
		0, 
function( g, h )
    local x, y, z;

     x := Log( g );
     y := Log( h );
     z := BCHStar( x, y );
     return Exp( z );
end);
     
InstallOtherMethod( \*, 
               "for symbolic Malcev group elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevGrpElement, IsMalcevGrpElement and IsSymbolicElement ],
		0, 
function( g, h )
    local x, y, z;

     x := Log( g );
     y := Log( h );
     z := BCHStar( x, y );
     return Exp( z );
end);

#############################################################################
##
#M g^n  ................................... power of Malcev group elements
##
InstallOtherMethod( \^, 
               "for Malcev group elments and an integer (Guarana)",
	        true,
	        [IsMalcevGrpElement, IsInt ],
		0, 
function( g, n )
    local malObj, exp_g, coll, gg, exp_res;

    malObj := g!.malcevObject;
    exp_g := Exponents( g );

    # get corresponing pcp elms
    coll := Collector( malObj!.recTGroup.T );
    gg := PcpElementByExponentsNC( coll, exp_g );

    exp_res := Exponents( gg^n );
	return MalcevGrpElementByExponents( malObj, exp_res );
end);

InstallOtherMethod( \^, 
               "for symbolic Malcev group elments and an integer (Guarana)",
	        true,
	        [IsMalcevGrpElement and IsSymbolicElement, IsInt ],
		0, 
function( g, n )
    local x;

    x := Log( g );
    return Exp( n*x );

end);

#############################################################################
##
#M Comm( g, h ) ........................ Commutator for Malcev group elements
##
InstallOtherMethod( COMM, 
               "for Malcev group elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevGrpElement, IsMalcevGrpElement ],
		0, 
function( g, h )
    return g^-1*h^-1*g*h;
end);

#############################################################################
##
#M g = h  .......................................... for Malcev grp elements
##
InstallOtherMethod( \=, 
               "for Malcev Grp elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevGrpElement, IsMalcevGrpElement ],
		0, 
function( x, y )
    return Exponents( x ) = Exponents( y );
end);

InstallOtherMethod( \*, 
               "for Malcev Gen elements and a matrix (Guarana)", 
               true, 
               [IsMalcevElement, IsMatrix ], 
               0,
function( elm, mat )
    local coeffs, coeffs_res;
    coeffs := Coefficients( elm );
    coeffs_res := coeffs * mat;
    return MalcevGenElementByCoefficients( elm!.malcevObject, coeffs_res );
end);

InstallOtherMethod( \^, 
               "for IsMalcevGenElement and a scalar (Guarana)",
	       true,
	        [IsMalcevGenElement, IsScalar ],
		0, 
function( elm, k )
    return k*elm;
end);

InstallOtherMethod( \*, 
               "product of scalar and Malcev Gen elments (Guarana)",
	       true,
	        [IsScalar, IsMalcevGenElement ],
		0, 
function( k, elm )
    local x, x_res;
    x := LieElement( elm );
    x_res := k*x;
    return MalcevGenElementByLieElement( x_res );
end);

#############################################################################
##
#M x = y 
##
InstallOtherMethod( \=, 
               "for Malcev Gen elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevGenElement, IsMalcevGenElement ],
		0, 
function( x, y )
    return Coefficients( x ) = Coefficients( y );
end);

if false then 
    malObjs := GUARANA.Get_FNG_MalcevObjects( 2, 4 );
    malObj := malObjs[3];
    x := MalcevLieElementByCoefficients( malObj, [1,2,3,4,5] );
    y := MalcevLieElementByCoefficients( malObj, [1,2,3,4,5] );
    LieBracket( x, y );

    n := 5;
    vars_x := GUARANA.RationalVariableList( n, "x" );
    vars_y := GUARANA.RationalVariableList( n, "y" );
    z := MalcevSymbolicLieElementByCoefficients( malObj, vars_x );
    o := MalcevSymbolicLieElementByCoefficients( malObj, vars_y );

    g := MalcevGrpElementByExponents( malObj, [1,2,3,4,5] );
    h := MalcevGrpElementByExponents( malObj, [1,0,0,0,0] );
    

fi;
#############################################################################
##
#E
