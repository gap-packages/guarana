#############################################################################
##
#W setup.gi               GUARANA package                     Bjoern Assmann
##
## Methods for the setup of the Malcev correspondence between the 
## radicable hull of a T-group and its corresponding Lie algebra. 
## For this setup we use the Bch formula (and not a unitriangular
## matrix representation). 
## We assume that certain informations about the T-group are given;
## in particular we assume that the pcp is given with respect to 
## a Mal'cev basis. These information about the T-group are stored
## in a record named recTGroup. 
##
#H  @(#)$Id$
##
#Y 2006
##
##

GUARANA.MultMethodIsStar := "Star";
GUARANA.MultMethodIsCollection := "Collec";

############################################################################
##
#F GUARANA.TGroupRec( args )
#F GUARANA.TGroupRec_GEN( N )
## 
## IN
## rgs[1]=N ................ a T-group given by a pcp
## args[2] .................  an optional string, that determines
##                            how the Mal'cev basis should be obtained.
##                            "gen" means that the Malcev basis is 
##                            assumed to be given. 
##                            "ucs" means that a Mal'cev basis going
##                            through the upper central series is computed
##
## OUT
## a record containing some information about N.
## For example weights of the elms of a Mal'cev basis of N.
## 
## Example of usage
## gap> ll:= GUARANA.SomePolyMalcevExams( 2 );
## gap> R := GUARANA.SetupCollecRecord( ll );
## gap> GUARANA.TGroupRec_GEN( R.T );
##
GUARANA.TGroupRec_GEN := function( N )
  local coll, weights, max_weight;

    # get weights 
    coll := Collector( N );
    if FromTheLeftCollector_SetWeights( coll ) <> fail then 
	weights := coll![PC_WEIGHTS];
    else
	Error( "Can not compute weights of Malcev basis of N" );
    fi;
    max_weight := Maximum( weights );

    return rec( T := N, weights := weights, 
                max_weight := max_weight,
                malcevBasisInfo := "gen" );
end;

GUARANA.TGroupRec := function( args )
    local N,malcevBasisInfo;

    N := args[1];
    
    # make choice how the record is set up.
    if IsBound( args[2] ) then 
	malcevBasisInfo := args[2];
    else
	#use default, i.e. the general method
	malcevBasisInfo := "gen";
    fi;

    if malcevBasisInfo = "gen" then
	return GUARANA.TGroupRec_GEN( N );
    elif malcevBasisInfo = "ucs" then 
        return GUARANA.TGroupRec_UCS( N );
    else
	Error( "wrong malcevBasisInfo specified" );
    fi;
end;

InstallGlobalFunction( MalcevObjectConstruction, 
function( recTGroup ) 
    local hl, T, L, malcevBasisInfo, lie_fam, lie_elms_type, grp_fam, 
          grp_elms_type, obj, gen_fam, gen_elms_type;
    
    # get dimension of algebra
    hl := HirschLength( recTGroup.T );

    # get prototype for structure constant table and Lie algebra 
    T:= EmptySCTable( hl, 0, "antisymmetric" );
    L:= LieAlgebraByStructureConstants( Rationals, T );

    # get information about Mal'cev basis.
    # This entry has influence how the Structure constants and Exp
    # are computed. 
    malcevBasisInfo := StructuralCopy( recTGroup.malcevBasisInfo );

    # create family and types for elements of Lie algebra and group
    gen_fam := NewFamily( "MalcevGenFamily",
                          IsMalcevGenElement,
                          IsMalcevGenElement );
    gen_elms_type := NewType( gen_fam, IsMalcevGenElementRep );
    lie_fam := NewFamily( "MalcevLieFamily", 
                           IsMalcevLieElement, 
			   IsMalcevLieElement );
    lie_elms_type := NewType( lie_fam, IsMalcevLieElementRep );
    grp_fam := NewFamily( "MalcevGrpFamily", 
                           IsMalcevGrpElement, 
			   IsMalcevGrpElement );
    grp_elms_type := NewType( grp_fam, IsMalcevGrpElementRep );

    obj := rec( L := L, 
                dim := hl,
                recTGroup := recTGroup, 
		scTable := T,
	        lieAlgebraType := "structureConstants", 
	        weights := recTGroup.weights,
		max_weight := recTGroup.max_weight,
	        malcevBasisInfo := malcevBasisInfo,
	        log_method := "simple",
		exp_method := "simple",
		star_method := "simple",
        mult_method := GUARANA.MultMethodIsStar,
        gen_fam := gen_fam,
        gen_elms_type := gen_elms_type,
		lie_fam := lie_fam,
		lie_elms_type := lie_elms_type,
		grp_fam := grp_fam,
		grp_elms_type := grp_elms_type );

    # compute structure constants
    obj :=  Objectify( NewType( MalcevObjectFamily, 
                               IsMalcevObjectRep and
			       IsMutable ),
		      obj );
    GUARANA.MO_ComputeStructureConstants(  obj );
    return obj;
end );

##
## some checks whether the input pcp group is defined with respect to 
## a Mal'cev basis
## 
GUARANA.IsGivenWithRespectToMalcevBasis := function( N )
    local C, rels, n, r, i, j, conj;

    if not IsPcpGroup( N ) then
        return fail;
    fi;

    C := Collector( N );

    # check wheter all relative orders are infinite
    rels := RelativeOrders( C );
    for r in rels do
        if r <> 0 then 
            return false;
        fi;
    od;

    # check whether nilpotent presentation
    n := NumberOfGenerators( C );
    for i in [1..n] do
        for j in [i+1..n] do
            conj := GetConjugate( C, j, i );
            if conj[1] <> j then
                return false;
            fi;
            if conj[2] <> 1 then 
                return false;
            fi;
        od;
    od;

    return true;
end;


## IN N ..................... T-group that is given by a pcp with respect 
##                            to Malcev basis.
##
InstallGlobalFunction( MalcevObjectByTGroup, 
function( N ) 
    local recT;
    if not GUARANA.IsGivenWithRespectToMalcevBasis( N ) then 
        return fail;
    fi;
    recT := GUARANA.TGroupRec( [N] );
    return MalcevObjectConstruction( recT );
end);

#############################################################################
##
##
InstallOtherMethod( UnderlyingLieAlgebra, "for Malcev objects", 
[ IsMalcevObjectRep ],
function( malObj )
    return malObj!.L;
end );

#############################################################################
##
##
InstallOtherMethod( UnderlyingGroup, "for Malcev objects", 
[ IsMalcevObjectRep ],
function( malObj )
    return malObj!.recTGroup.T;
end );

#############################################################################
##
##
InstallOtherMethod( Dimension, "for Malcev objects", 
[ IsMalcevObjectRep ],
function( malObj )
    return malObj!.dim;
end );

#############################################################################
##
## Funtions to view and print a Malcev object.
##
InstallMethod( ViewObj, "for Malcev object", [ IsMalcevObjectRep ],
function( malObj )
    Print( "<<Malcev object of dimension ",
           malObj!.dim,
	   ">>" );
end );

InstallMethod( PrintObj, "for Malcev object", [IsMalcevObjectRep ],
function( malObj )
    Print( "<<Malcev object of dimension ",
           malObj!.dim,
	   ">>" );
end );

#############################################################################
##
#F GUARANA.EvaluateGroupCommutator( g, h, com )
##
## IN
## g,h ................................... group elements
## com ..................................  list containing 1,2
##                                         g corresponds to 1
##                                         h corresponds to 2 
##
## OUT
## com(g,h)
##
GUARANA.EvaluateGroupCommutator := function( g, h, com )
    local tmp,r,l,i;
    tmp := [g,h];
    r := tmp[com[1]];

    l := Length( com );
    for i in [2..l] do
        r := Comm( r, tmp[com[i]] );
    od;
    return r;
end;

#############################################################################
##
#F GUARANA.WeightOfCommutator( com, wx, wy )
## 
## IN 
## com ................. commutator in x,y
## wx .................. weight of x
## wy .................  weight of y 
## 
## OUT
## Weight of com(x,y).
##
GUARANA.WeightOfCommutator := function( com, wx, wy )
    local weight,a;
    weight := 0;
    for a in com do
        if a = 1 then
            weight := weight + wx;
        elif a = 2 then
            weight := weight + wy;
        else
            Error( "Wrong input\n" );
        fi;
    od;
    return weight;
end;

#############################################################################
##
##
GUARANA.CheckWeightOfCommutator := function( com, wx, wy, max_weight )
    local w;
    w := GUARANA.WeightOfCommutator( com, wx, wy );
    if w > max_weight then
        return false;
    else
        return true;
    fi;
end;

GUARANA.CheckWeightOfMalcevCommutator := function( x, y, com )
    local wx, wy, max_weight, w;

    wx := Weight( x );
    wy := Weight( y );
    max_weight := x!.malcevObject!.max_weight;

    w := GUARANA.WeightOfCommutator( com, wx, wy );
    if w > max_weight then
        return false;
    else
        return true;
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

    n := 5;
    vars_x := GUARANA.RationalVariableList( n, "x" );
    vars_y := GUARANA.RationalVariableList( n, "y" );
    z := MalcevSymbolicLieElementByCoefficients( malObj, vars_x );
    o := MalcevSymbolicLieElementByCoefficients( malObj, vars_y );

    gg := MalcevSymbolicGrpElementByExponents( malObj, vars_x );

fi;

#############################################################################
##
#M Log  ............................................. for Malcev grp elments 
##
GUARANA.LogByStar := function( g )
    local exp, malcevObject, r, e, x, y, i;

    exp := Exponents( g );
    malcevObject := g!.malcevObject;

    # get zero element of lie algebra of correct type (symbolic/not symbolic)
    if IsSymbolicElement( g ) then 
	r := MalcevSymbolicLieElementByWord( malcevObject,  [[],[]] );
    else
	r := MalcevLieElementByWord( malcevObject, [[],[]] );
    fi;
	
    # go backwards through exponents
    for i in Reversed([1..Length( exp )] ) do
        e := exp[i];
        if e <> 0 then 
	    x := GUARANA.MalcevBasisLieElement( malcevObject, i, e );
            y := r;
	    r := BCHStar( x, y );
        fi;
    od;
    return r;
end;

InstallOtherMethod( Log, 
               "for Malcev group elments (Guarana)",
	       true,
	        [IsMalcevGrpElement ],
		0, 
function( g )
    if g!.malcevObject!.log_method = "pols" then 
	return GUARANA.LogByPols( g );
    else 
	return GUARANA.LogByStar( g );
    fi;
end);

#############################################################################
##
#F  GUARANA.MO_EvaluateLieBracketsInTermsOfLogarithmsSers 
##
## IN
## g,h .................................................. Malcev grp elments
##
## OUT 
## [Log g, Log h] 
## This is computed by using an  identity that expresses
## Lie brackets as a linear combination of logarithms of group
## commutators. 
##
GUARANA.MO_EvaluateLieBracketsInTermsOfLogarithmsSers 
                            := function(  g,h )
    local malcevObject, bchLBITOL, r, tree, max_weight, wg, wh, max,  
          min, bound, com, a, log_a, i, term;

    # catch trivial case
    malcevObject := g!.malcevObject;
    if Exponents( g ) = [] then 
	return MalcevLieElementByCoefficients( malcevObject, [] );
    fi;

    bchLBITOL := GUARANA.recBCH.bchLBITOL;
 
    # get zero elment of Lie algebra
    r := 0*GUARANA.MalcevBasisLieElement( malcevObject, 1, 0 );

    # set up tree 
    tree := [h];

    # compute upper bound for the Length of commutators, which 
    # can be involved
    max_weight := malcevObject!.max_weight;
    wg := Weight( g );
    wh := Weight( h );
    max := Maximum( wg,wh );
    min := Minimum( wg,wh );
    # max + min* (bound-1 ) <= max_weight
    bound := Int( (max_weight-max)/min + 1 );

    # up to bound  compute the commutators and add them.
    # Note that the list contains commutators of length i at position i-1.
    for i in [1..bound-1] do
        for term in bchLBITOL[i] do
            com := term[2];
            # check if weight of commutator is not to big
            if GUARANA.CheckWeightOfCommutator( com, wg, wh, max_weight ) then
                # evaluate commutator in the group
		a := Comm( g, h, com, tree );
                # map to the Lie algebra
		log_a := Log( a ); 
                r := r + term[1]*log_a; 
            fi;
        od;
    od;
    return r;
end;

#############################################################################
##
#F GUARANA.MO_LieAlgElm2CoeffGenList( x )
##
## IN
## x .............. malcev lie element
##
## OUT 
## a list, as required for SetEntrySCTable
##
GUARANA.MO_LieAlgElm2CoeffGenList := function( x )
    local coeffs, ll, i;
    coeffs := Coefficients(  x );
    ll := [];
    for i in [1..Length(coeffs)] do
        if coeffs[i] <> 0 then
            Append( ll, [coeffs[i],i] );
	fi;
    od;
    return ll;
end;

#############################################################################
##
#F GUARANA.MO_ComputeStructureConstants( malcevObject )
##
##
## This function computes the structure constants of the Lie algebra.
##
GUARANA.MO_ComputeStructureConstants := function( malcevObject )
    local dim, T, g, h, lie_elm, ll, index_y, index_x;

    # setup 
    dim := malcevObject!.dim;
    T := malcevObject!.scTable;

    # go through generators backwards and compute the 
    # structure constants.  
    for index_y in Reversed( [1..dim-1] ) do 
	for index_x in [1..index_y-1]  do
	    g := GUARANA.MalcevBasisGrpElement( malcevObject, index_x, 1 );
	    h := GUARANA.MalcevBasisGrpElement( malcevObject, index_y, 1 );

            # compute [log(g),log(h)]
            lie_elm := GUARANA.MO_EvaluateLieBracketsInTermsOfLogarithmsSers( 
                                                                       g,h );
            # if non trivial, then set entry in structure constant table 
	    if lie_elm <> 0*lie_elm then 
                ll := GUARANA.MO_LieAlgElm2CoeffGenList( lie_elm );
                SetEntrySCTable( T, index_x, index_y, ll );
            fi;
	od;
        # update Lie algebra with new structure constant table
        malcevObject!.L:= LieAlgebraByStructureConstants( Rationals, T );
    od;
    return 0;
end;

#############################################################################
##
#M Exp  ............................................. for Malcev lie elments 
##
GUARANA.ExpByStar := function( x )
    local malcevObject, coeffs, tail, largestAbelian, exp_x, divider, 
          l, exp_x_2ndPart, exp, j;

    # catch trivial case
    malcevObject := x!.malcevObject;
    coeffs := Coefficients( x );

    tail := x;
    # get smallest index i such that n_i,...,n_l generate an abelian group
    # TODO This can be done better. 
    largestAbelian := malcevObject!.dim -1 ; 

    exp_x := [];
    for j in [1..largestAbelian-1] do
        # get element to divide of
	divider := GUARANA.MalcevBasisLieElement( malcevObject, j, -coeffs[j]);

        # save exponent of divider
        Add( exp_x, coeffs[j] );

        # divide off
	tail := BCHStar( divider, tail );
        
        # set up coefficient vector
        coeffs := Coefficients( tail );
    od;

    # test intermediate result
    l := Length( exp_x );
    if not coeffs{[1..l]} = 0 *  coeffs{[1..l]} then
        Error( "Failure in the computation of Exp \n" );
    fi;

    # get the remaining coefficients 
    exp_x_2ndPart := coeffs{[l+1..Length(coeffs)]};
    
    exp := Concatenation( exp_x, exp_x_2ndPart );
    if IsSymbolicElement( x ) then 
	return MalcevSymbolicGrpElementByExponents( malcevObject, exp );
    else
        return MalcevGrpElementByExponents( malcevObject, exp );
    fi;
end;

InstallMethod( Exp,
               "for Malcev lie elments (Guarana)",
	       true,
	        [IsMalcevLieElement ],
		0, 
function( x )
    if x!.malcevObject!.exp_method = "pols" then 
	return GUARANA.ExpByPols( x );
    else
	return GUARANA.ExpByStar( x );
    fi;
end);


#############################################################################
##
#E 
