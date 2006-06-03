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
GUARANA.MultMethodIsDeepThought := "DeepThought";

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
## Old code 
## 

############################################################################
##
#F GUARANA.TGroupRec_UCS( N )
## 
## IN
## N ................ a T-group given by a pcp
##
## OUT
## a record containing some information about N. 
## It contains a Mal'cev basis of N  going through the 
## upper central series and a  
## pcp with respect to that basis. 
## 
## Example of usage
## gap> N := GUARANA.Examples_Unitriangular_UCS( 4, 2 );
## Pcp-group with orders [ 0, 0, 0, 0, 0, 0 ]
## gap> GUARANA.TGroupRec( N );
##
## Note:
## This function is not so nice it because it changes the
## underlying pcs of N.
## This should not happen at this stage. 
## Therefore this function we will be only used for test purposes. 
##
GUARANA.TGroupRec_UCS := function( N )
    local sers,i, NN,l,leFactors,s,indices,max_weight,
    wei, weights,a, largestAbelian;
    sers := UpperCentralSeries( N );
    NN := PcpGroupBySeries( sers ); 
    #basisNN := GeneratorsOfGroup( NN );
    #iso := NN!.bijection;
    #basisN := GeneratorsOfGroup( Image( iso ) );

    # get indices of generators of the factors of the upper central series
    l := Length( sers );
    leFactors := [];
    for i in [1..l-1] do
        Add( leFactors, HirschLength( sers[i]) - 
	HirschLength( sers[i+1]));
    od;
    s := 1;
    indices := [];
    for i in [1..l-1] do
	Add( indices, [s..s+leFactors[i]-1] );
	s := s + leFactors[i];
    od;

    # weights of generators
    wei := 1;
    weights := [];
    for a in leFactors do
	for i in [1..a] do
	    Add( weights, wei );
	od;
	wei := wei + 1;
    od;

    # max_weight of elms of the Malcev basis 
    max_weight := Length( leFactors );

    # get largest abelian group in the upper central series
    for i in [1..Length(sers)] do
	if IsAbelian( sers[i] ) then
	    largestAbelian := i;
	    break;
	fi;
    od;

    return rec( N := N, #basisN := basisN, 
                T := NN,
	        sers := sers, indices := indices,
	        max_weight := max_weight, weights := weights,
	        largestAbelian := largestAbelian,
		# information about how the Mal'cev basis was computed
	        malcevBasisInfo := "ucs" );
end;

#############################################################################
##
#F GUARANA.SetUpLieAlgebraRecordByMalcevbasis( recTGroup )
##
## IN
## recTGroup .............. record containing some information about
##                          a T-group. In particular a pcp with 
##                          respect to a Mal'cev basis. 
##
## OUT
## An empty data structure that contain be filled with information
## about the corresponding Lie algebra. 
##
GUARANA.SetUpLieAlgebraRecordByMalcevbasis := function( recTGroup )
    local hl, T, L, malcevBasisInfo, lie_fam, lie_elms_type;
    
    # get dimension of algebra
    hl := HirschLength( recTGroup.T );

    # get prototype for structure constant table and Lie algebra 
    T:= EmptySCTable( hl, 0, "antisymmetric" );
    L:= LieAlgebraByStructureConstants( Rationals, T );

    # get information about Mal'cev basis.
    # This entry has influence how the Structure constants and Exp
    # are computed. 
    malcevBasisInfo := StructuralCopy( recTGroup.malcevBasisInfo );

    # create family and types for elements of Lie algebra
    lie_fam := NewFamily( "MalcevLieFamily", 
                           IsMalcevLieElement, 
			   IsMalcevLieElement );
    lie_elms_type := NewType( lie_fam, IsMalcevLieElementRep );
    

    return rec( L := L, 
                dim := hl,
                recTGroup := recTGroup, 
		scTable := T,
	        weights := recTGroup.weights,
		max_weight := recTGroup.max_weight,
	        malcevBasisInfo := malcevBasisInfo,
	        log_method := "simple",
		exp_method := "simple",
		star_method := "simple",
		lie_fam := lie_fam,
		lie_elms_type := lie_elms_type );
end;

#############################################################################
##
#F GUARANA.EvaluateLieBracket( x,y, com, info )
##
## IN
## x,y .................. elements of Lie algebra.
## com .................  commutator given as list containing 1,2.
##                        x corresponds to 1
##                        y corresponds to 2 
##                        For example [1,2,1] corresponds to 
##                        the commutator [x,y,x]
## info ................  string that is either "strucConst" or 
##                        "matrix". In the first case the 
##                        x,y are elements of a Lie algebra given
##                        by strcuture constants. In the second 
##                        case x,y are matrices. 
## 
## OUT
## The result of Com( x,y ). 
##
GUARANA.EvaluateLieBracket := function( x, y, com, info )
    local r,l,tmp,i;
    tmp := [x,y];
    r := tmp[com[1]];

    l := Length( com );
    if info = "strucConst" then
        for i in [2..l] do
            r := r * tmp[com[i]];
        od; 
    elif info = "matrix" then 
        for i in [2..l] do
            r := LieBracket( r, tmp[com[i]] );
        od;
    else 
        Error( "wront input \n " );
    fi;
    return r;
end;

#############################################################################
##
#F GUARANA.Star_Simple( x, y, wx, wy, max_weight, info )
##
## IN
## x,y ................... elements of a Lie algebra.
##                         Brackets of Length max_weight + 1 are always 0.
## wx,wy ................. their weights. 
## max_weight .............max_weight of the basis elms of the Lie algebra. 
## info .................  strings which determines how lie brackets
##                         are evaluated
##
## OUT 
## x*y, where "*" is the BCH operation 
##
## This is a very simple implementation.
## TODO 
## Can this be done in a better way ? Maybe with a nice data 
## structure containing the information about the Bch-formula. 
## 
GUARANA.Star_Simple := function( x, y, wx, wy, max_weight, info  )
    local i,r,bchSers,com,a,term,max,min,bound;
    bchSers := GUARANA.recBCH.bchSers;
    
    # start with terms which are not given by Lie brackets
    r := x + y;

    # trivial check 
    if x = 0*x or y = 0*y then
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
                a := GUARANA.EvaluateLieBracket( x, y, com, info );
                r := r + term[1]*a; 
            fi;
        od;
    od;
    return r;
end;

#############################################################################
##
##
GUARANA.WeightOfLieAlgElm := function( recLieAlg, elm )
    local basisL,coeff,i,e;
    basisL := Basis( recLieAlg.L );
    coeff := Coefficients( basisL, elm );
    for i in [1..Length(coeff)] do
        e := coeff[i];
        if e <> 0 then
            return recLieAlg.weights[i];
        fi;
    od;
    return recLieAlg.max_weight;
end;
#############################################################################
##
#F GUARANA.AbstractLog_Simple_ByExponent( recLieAlg, exp )
## 
## IN
## recLieAlg ............... Lie algebra record
## exp ..................... exponent vector of a group element g. 
## 
## OUT 
## Log( g )
##
GUARANA.AbstractLog_Simple_ByExponent := function( recLieAlg, exp )
    local  gensL,r,i,e,x,y,wx,wy,max_weight;

    # setup
    gensL := GeneratorsOfAlgebra( recLieAlg.L );
    max_weight := recLieAlg.max_weight;

    r := Zero( recLieAlg.L );
    # go backwards through exponents. 
    # TODO  Which direction is more efficient ?
    for i in Reversed([1..Length( exp )] ) do
        e := exp[i];
        if e <> 0 then 
            x := e*gensL[i];
            y := r;
            # if y is the identity then x*y = x
            if y = Zero( recLieAlg.L ) then
                r := x;
            else
                # compute weights of x,y
                wx := recLieAlg.weights[i];
                wy := GUARANA.WeightOfLieAlgElm( recLieAlg, y );
                r := GUARANA.Star_Simple( x,y,wx,wy,max_weight, "strucConst" );
            fi;
        fi;
    od;
    return r;
end;

#############################################################################
##
#F GUARANA.AbstractLogCoeff_Simple_ByExponent( recLieAlg, exp )
## 
## IN
## recLieAlg ............... Lie algebra record
## exp ..................... exponent vector of a group element g. 
## 
## OUT 
## Log( g ) given as a coefficient vector. 
##
GUARANA.AbstractLogCoeff_Simple_ByExponent := function( recLieAlg,  exp )
    local r;
    r := GUARANA.AbstractLog_Simple_ByExponent( recLieAlg,   exp );
    return Coefficients( Basis( recLieAlg.L ), r);
end;

#############################################################################
##
#F GUARANA.AbstractLog_Simple_ByElm( recLieAlg, g )
## 
## IN
## recLieAlg ............... Lie algebra record
## g .....................   group element g. 
## 
## OUT 
## Log( g )
##
## TODO 
## adapt this function. Is it still necessary to extract the exponent 
## vector in this way ?
##
GUARANA.AbstractLog_Simple_ByElm := function( recLieAlg, g )
    local exp,hl,l,expNN;
    
    # get exponents with respect to gens of T group T 
    # (g might an elment of PcpGroup which contains T  as a normal sugroup,
    # then Exponents(g) would be longer hl(T ))
    exp := Exponents( g );
    #hl := HirschLength( recLieAlg.recTGroup.T );
    #l := Length( exp );
    #expNN := exp{[l-hl+1..l]};
    expNN := exp;
    return GUARANA.AbstractLog_Simple_ByExponent( recLieAlg, expNN );
end;

#############################################################################
##
#F GUARANA.AbstractLogCoeff_Simple_ByElm( recLieAlg, g )
## 
## IN
## recLieAlg ............... Lie algebra record
## g .....................   group element g. 
## 
## OUT 
## Log( g ) given by a coefficient vector.
##
GUARANA.AbstractLogCoeff_Simple_ByElm := function( recLieAlg, g )
    local r;
    r := GUARANA.AbstractLog_Simple_ByElm( recLieAlg, g );
    return Coefficients( Basis( recLieAlg.L ), r);
end;
#############################################################################
##
#F  GUARANA.EvaluateLieBracketsInTermsOfLogarithmsSers 
##
## IN
## recLieAlg ................... Lie algebra record
## g,h ......................... group elements 
## wg,wh ........................weights of g,h
##
## OUT 
## [Log g, Log h] 
## This is computed by using an  identity that expresses
## Lie brackets as a linear combination of logarithms of group
## commutators. 
##
GUARANA.EvaluateLieBracketsInTermsOfLogarithmsSers 
                            := function(  recLieAlg, g,h,wg,wh )
    local bchLBITOL,r,max,min,bound,max_weight,i,term,com,a,log_a;
    bchLBITOL := GUARANA.recBCH.bchLBITOL;
  
    r := Zero( recLieAlg.L );

    # compute upper bound for the Length of commutators, which 
    # can be involved
    max_weight := recLieAlg.max_weight;
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
                a := GUARANA.EvaluateGroupCommutator( g, h, com );
                # map to the Lie algebra
                log_a := GUARANA.AbstractLog_Simple_ByElm( recLieAlg,  a );
                r := r + term[1]*log_a; 
            fi;
        od;
    od;
    return r;
end;
#############################################################################
##
#F GUARANA.ComputeStructureConstants_UCS( recLieAlg )
##
## IN  
## recLieAlg ........... lie algebra record. 
##
## OUT
## 0
##
## This function computes the structure constants of the Lie algebra
## and stores them in the Lie algebra record.
## We assume here that the Mal'cev basis goes through the upper
## central series of N. This knowledge is used to reduce 
## the number of structure constants that has to be computed. 
##
GUARANA.ComputeStructureConstants_UCS := function( recLieAlg )
    local factors, index_x, index_y, indices,indicesNoCenter,l,g,h,wg,wh,
          lie_elm,ll,gens,T,out;

    # setup
    indices := recLieAlg.recTGroup.indices;
    l := Length( indices );
    indicesNoCenter := indices{[1..l-1]};
    gens := GeneratorsOfGroup( recLieAlg.recTGroup.T  );
    T := recLieAlg.scTable;

    # go through blocks of generators backwards, 
    # where a block corresponds to the generators of the factors
    # of the upper central series
    for factors in Reversed( indicesNoCenter ) do
        for index_y in Reversed( factors ) do
            for index_x in  [1..index_y-1]   do
                g := gens[index_x];
                h := gens[index_y];
                wg := recLieAlg.weights[index_x];
                wh := recLieAlg.weights[index_y];
                # compute [log(g),log(h)]
                lie_elm := GUARANA.EvaluateLieBracketsInTermsOfLogarithmsSers( 
                                             recLieAlg, g,h, wg, wh );
                # set entry in structure constant table 
                ll := GUARANA.LieAlgElm2CoeffGenList( recLieAlg.L, lie_elm );
                if Length( ll ) > 0 then
                    SetEntrySCTable( T, index_x, index_y, ll );
                fi;
            od;
        od;
        # update Lie algebra with new structure constant table
        recLieAlg.L:= LieAlgebraByStructureConstants( Rationals, T );
    od;
    return 0;
end;

#############################################################################
##
#F GUARANA.ComputeStructureConstants( recLieAlg )
##
## IN  
## recLieAlg ........... lie algebra record. 
##
## OUT
## 0
##
## This function computes the structure constants of the Lie algebra
## and stores them in the Lie algebra record.
##
##
GUARANA.ComputeStructureConstants := function( recLieAlg )
    # decide which method should be used. 
    if recLieAlg.malcevBasisInfo = "ucs" then
	return GUARANA.ComputeStructureConstants_UCS( recLieAlg );
    elif recLieAlg.malcevBasisInfo = "gen" then 
	return GUARANA.ComputeStructureConstants_GEN( recLieAlg );
    fi;   
    Error( "Unknown method" );
end;

## Example of usage
if false then 
    F := FreeGroup( 2 );
    N := NilpotentQuotient( F, 3 );
    recTGroup := GUARANA.TGroupRec( [N] );
    recLieAlg := GUARANA.SetUpLieAlgebraRecordByMalcevbasis( recTGroup );

    GUARANA.ComputeStructureConstants( recLieAlg );

    # try general method 
    recLieAlg.malcevBasisInfo := "gen";
    GUARANA.ComputeStructureConstants( recLieAlg );
    
    N := PcpGroupByMatGroup( PolExamples(4 ) );
    recTGroup := GUARANA.TGroupRec([ N] );  
    recLieAlg := GUARANA.SetUpLieAlgebraRecordByMalcevbasis( recTGroup );
fi;

#############################################################################
##
#F GUARANA.LieAlgElm2CoeffGenList( L, x )
##
## IN
## L .............. Lie algebra 
## x .............. elments of L
##
## OUT 
## a list, as required for SetEntrySCTable
##
GUARANA.LieAlgElm2CoeffGenList := function( L, x )
    local basis, coeff,i,ll;
    basis := Basis( L );
    coeff := Coefficients( basis, x );
    ll := [];
    for i in [1..Length(coeff)] do
        if coeff[i] <> 0 then
            Append( ll, [coeff[i],i] );
        fi;
    od;
    return ll;
end;
#############################################################################
##
#F GUARANA.ComputeStructureConstants_GEN( recLieAlg )
##
## IN  
## recLieAlg ........... lie algebra record. 
##
## OUT
## 0
##
## This function computes the structure constants of the Lie algebra
## and stores them in the Lie algebra record.
## It does not assume any any additonal knowledge about the 
## Mal'cev basis, apart from the fact that we need weights 
## of the elements of the Mal'cev basis. 
##
## TODO 
## check wether this function works fine.
##
GUARANA.ComputeStructureConstants_GEN := function( recLieAlg )
    local gens,n,T,index_y,index_x,g,h,wg,wh,lie_elm,ll; 

    # setup
    gens := GeneratorsOfGroup( recLieAlg.recTGroup.T );
    n := Length( gens );
    T := recLieAlg.scTable;

    # go through generators backwards and compute the 
    # structure constants.  
    for index_y in Reversed( [1..n-1] ) do 
	for index_x in [1..index_y-1]  do
            g := gens[index_x];
            h := gens[index_y];
            wg := recLieAlg.weights[index_x];
            wh := recLieAlg.weights[index_y];
            # compute [log(g),log(h)]
            lie_elm := GUARANA.EvaluateLieBracketsInTermsOfLogarithmsSers( 
                                             recLieAlg, g,h, wg, wh );
            # set entry in structure constant table 
            ll := GUARANA.LieAlgElm2CoeffGenList( recLieAlg.L, lie_elm );
            if Length( ll ) > 0 then
                SetEntrySCTable( T, index_x, index_y, ll );
            fi;
	od;
        # update Lie algebra with new structure constant table
        recLieAlg.L:= LieAlgebraByStructureConstants( Rationals, T );
    od;
    return 0;
end;
#############################################################################
##
#F GUARANA.Abstract_Exponential_ByElm_UCS( recLieAlg, x )
##
## IN 
## recLieAlg ................. Lie algebra record
## x  ........................ element of Lie algebra
##
## OUT
## Exp( x )
## Here we assume that the underlying Mal'cev basis goes through
## the upper central series.
##
GUARANA.Abstract_Exponential_ByElm_UCS := function(  recLieAlg, x )
    local indices,basis,max_weight,tail,coeffs,largestAbelian,exp_x,i,factor,
          divider,w_divider,w_tail,l,exp_x_2ndPart,j;

    indices := recLieAlg.recTGroup.indices;
    basis := Basis( recLieAlg.L );
    max_weight := recLieAlg.max_weight;

    # divide some exponents untill the remaining element lies in an
    # abelian group.
    tail := x;
    coeffs := Coefficients( basis, tail );
    largestAbelian := recLieAlg.recTGroup.largestAbelian;
    exp_x := [];
    for i in [1..largestAbelian-1] do
        factor := indices[i];
        for j in factor do 
            # get element to divide of
            divider := -coeffs[j]*basis[j];

            # save exponents of divider
            Add( exp_x, coeffs[j] );

            # divide off
            w_divider := GUARANA.WeightOfLieAlgElm ( recLieAlg, divider );
            w_tail := GUARANA.WeightOfLieAlgElm ( recLieAlg, tail );
            tail := GUARANA.Star_Simple(  divider, tail, w_divider, 
                                     w_tail, max_weight, "strucConst"  );
        
            # set up coefficient vector
            coeffs := Coefficients( basis, tail );
        od;
    od;

    # test intermediate result
    l := Length( exp_x );
    if not coeffs{[1..l]} = 0 *  coeffs{[1..l]} then
        Error( "Failure in Abstract_Exponential \n" );
    fi;

    # get the remaining coefficients 
    exp_x_2ndPart := coeffs{[l+1..Length(coeffs)]};
    
    return Concatenation( exp_x, exp_x_2ndPart );
    
end;
#############################################################################
##
#F GUARANA.Abstract_Exponential_ByElm_GEN( recLieAlg, x )
##
## IN 
## recLieAlg ................. Lie algebra record
## x  ........................ element of Lie algebra
##
## OUT
## Exp( x )
## GENeral method.
##
GUARANA.Abstract_Exponential_ByElm_GEN := function(  recLieAlg, x )
    local basis, max_weight, tail, coeffs, largestAbelian, exp_x, divider, 
          w_divider, w_tail, l, exp_x_2ndPart, j;

    basis := Basis( recLieAlg.L );
    max_weight := recLieAlg.max_weight;

    # divide some exponents untill the remaining element lies in an
    # abelian group.
    tail := x;
    coeffs := Coefficients( basis, tail );
    # smallest index i such that n_i,...,n_l generate an abelian group
    # TODO This can be done better. 
    largestAbelian := Length( basis ) -1 ; 
    exp_x := [];
    for j in [1..largestAbelian-1] do
        # get element to divide of
         divider := -coeffs[j]*basis[j];

        # save exponents of divider
        Add( exp_x, coeffs[j] );

        # divide off
        w_divider := GUARANA.WeightOfLieAlgElm ( recLieAlg, divider );
        w_tail := GUARANA.WeightOfLieAlgElm ( recLieAlg, tail );
        tail := GUARANA.Star_Simple(  divider, tail, w_divider, 
                                     w_tail, max_weight, "strucConst"  );
        
        # set up coefficient vector
        coeffs := Coefficients( basis, tail );
    od;

    # test intermediate result
    l := Length( exp_x );
    if not coeffs{[1..l]} = 0 *  coeffs{[1..l]} then
        Error( "Failure in Abstract_Exponential \n" );
    fi;

    # get the remaining coefficients 
    exp_x_2ndPart := coeffs{[l+1..Length(coeffs)]};
    
    return Concatenation( exp_x, exp_x_2ndPart );
end;

#############################################################################
##
#F GUARANA.Abstract_Exponential_ByElm( recLieAlg, x )
##
## IN 
## recLieAlg ................. Lie algebra record
## x  ........................ element of Lie algebra
##
## OUT
## Exp( x )
##
GUARANA.Abstract_Exponential_ByElm := function( recLieAlg, x )
    if recLieAlg.malcevBasisInfo = "ucs" then 
	return GUARANA.Abstract_Exponential_ByElm_UCS( recLieAlg, x );
    elif recLieAlg.malcevBasisInfo = "gen" then 
	return GUARANA.Abstract_Exponential_ByElm_GEN( recLieAlg, x );
    else
	Error( "Method not known" );
    fi;
end;

#############################################################################
##
#F GUARANA.Abstract_Exponential_ByVector( recLieAlg, x )
##
## IN 
## recLieAlg ................. Lie algebra record
## vec........................ coeff vector of an element of Lie algebra
##
## OUT
## Exp( x )
##
GUARANA.Abstract_Exponential_ByVector := function(  recLieAlg, vec )
    local basis,x;
    basis := Basis( recLieAlg.L );
    x := LinearCombination( basis, vec );
    return GUARANA.Abstract_Exponential_ByElm(  recLieAlg, x ); 
end;
#############################################################################
##
#F GUARANA.LieAlgebraByTGroupRec( args )
##
## IN
## args[1]=recTGroup .................... T-group record
## args[2]=malcevBasisInfo .............. String containing information
##                                        about the Malcev basis. 
##                                        This determines the way how
##                                        the structure constants and Exp
##                                        are computed. Currently 
##                                        there are 2 choices 
##                                        "gen" general
##                                        "ucs" upper central series. 
##
## OUT 
## Lie algebra record, that contains already the structure 
## constants of the Lie algebra. 
## 
GUARANA.LieAlgebraByTGroupRec := function( args )
    local recLieAlg,recTGroup;
    recTGroup := args[1];
    recLieAlg := GUARANA.SetUpLieAlgebraRecordByMalcevbasis( recTGroup );
    if IsBound( args[2] ) then
	recLieAlg.malcevBasisInfo := args[2];
    fi;
    GUARANA.ComputeStructureConstants(  recLieAlg );

    return recLieAlg;
end;
#############################################################################
##
#F GUARANA.AbstractLog( args )
##
## IN
## recLieAlg ..................... Lie Algebra record.
## g    .......................... Group element, either given as 
##                                 exponent vector or as an element of 
##                                 pcp group.
## rep_info ...................... indicates in which way the group element
##                                 is given and in which way the lie algbra 
##                                 element should be given. 
##                                 There are four possibilities:
##                                 "vecByVec"
##                                 "vecByElm"
##                                 "elmByVec"
##                                 "elmByElm"
## log_method .................... optional parameter. Decides with 
##                                 which method the logarithm 
##                                 should be computed. 
##                                 There are so far 2 possibilities
##                                 "simple"
##                                 "symbolic"
##              
## OUT
## Log( g ) either given as coefficient vector or as an elment of 
## a Lie algebra. 
##
## EXAMPLE
## GUARANA.AbstractLog( [malcevRec.recL_NN, g, "vecByVec"] ); 
##
GUARANA.AbstractLog := function( args )
    local recLieAlg, g, rep_info, log_method;

    recLieAlg := args[1];
    g := args[2];
    rep_info := args[3];
    if IsBound( args[4] ) then 
	log_method := args[4];
    else
	log_method := recLieAlg.log_method;
    fi;

    # catch trivial case
    if rep_info = "vecByVec" then 
	if Length( g ) = 0 then 
	    return [];
	fi;
    fi;

    if log_method = "simple" then 
	if rep_info = "elmByElm" then 
	    return GUARANA.AbstractLog_Simple_ByElm( recLieAlg, g );
	elif rep_info = "elmByVec" then 
	    return GUARANA.AbstractLog_Simple_ByExponent( recLieAlg, g );
	elif rep_info = "vecByVec" then
	    return GUARANA.AbstractLogCoeff_Simple_ByExponent( recLieAlg, g);
	elif rep_info = "vecByElm" then
	    return GUARANA.AbstractLogCoeff_Simple_ByElm( recLieAlg, g );
	else 
	    Error( "wrong rep_info" );
	fi;
    elif log_method = "symbolic" then
	if rep_info = "vecByVec" then 
	    return GUARANA.Logarithm_Symbolic( [recLieAlg,g] );
	else 
	    Error( "Sorry method not implemented yet" );
	fi;
    else
	Error( "Wrong log_method" );
    fi;
end;

#############################################################################
##
#F GUARANA.AbstractExp( args )
##
## IN
## recLieAlg ......................... Lie algebra record
## x ................................. Lie algebra elment etiher 
##                                     via a coefficient vector or as an
##                                     element of a Lie algebra (in the 
##                                     GAP sense).
## rep_info .......................... indicate in which way the Lie 
##                                     algebra elment is given and 
##                                     in which way the group element should
##                                     be given. There are 4 possibilities:
##                                     "vecByVec"
##                                     "vecByElm"
##                                     "elmByVec"
##                                     "elmByElm"
## exp_method ........................ optional parameter. Decides with 
##                                     which method the exponential should
##                                     be computed. So far there 2 possib.
##                                     "simple"
##                                     "symbolic"
## 
## OUT
## Exp( x ) either given as epxonent vector or as an element of a Pcp group
##
GUARANA.AbstractExp := function( args )
    local recLieAlg, x, rep_info, exp_method;

    recLieAlg := args[1];
    x := args[2];
    rep_info := args[3];
    if IsBound( args[4] ) then 
	exp_method := args[4];
    else
	exp_method := recLieAlg.exp_method;
    fi;

    if exp_method = "simple" then 
	if rep_info = "elmByElm" then 
	    Error( "Sorry method not implemented yet" );
	elif rep_info = "elmByVec" then 
	    Error( "Sorry method not implemented yet" );
	elif rep_info = "vecByVec" then
	    return GUARANA.Abstract_Exponential_ByVector( recLieAlg, x);
	elif rep_info = "vecByElm" then
	    return GUARANA.Abstract_Exponential_ByElm( recLieAlg, x );
	else 
	    Error( "wrong rep_info" );
	fi;
    elif exp_method = "symbolic" then
	if rep_info = "vecByVec" then 
	    return GUARANA.Exponential_Symbolic( [recLieAlg,x] );
	else 
	    Error( "Sorry method not implemented yet" );
	fi;
    else
	Error( "Wrong exp_method" );
    fi;
end;

## IN
## recLieAlg .............................. Lie algebra record of a 
##                                          Lie algebra L
## x,y .................................... elements of L 
## rep_info .......................... indicate in which way the Lie 
##                                     algebra elments are given.
##                                     There are two  possibilities:
##                                     "vec"
##                                     "elm"
## star_method ........................ optional parameter. Decides with 
##                                     which method star should
##                                     be computed. So far there 2 possib.
##                                     "simple"
##                                     "symbolic"
## 
## OUT
## x*y either given in the form x,y are given.
##
GUARANA.Star := function( args )
    local recLieAlg, x, y, rep_info, star_method, basis, xx, yy, 
          wx, wy, max_w, xx_star_yy, x_star_y;

    recLieAlg := args[1];
    x := args[2];
    y := args[3];
    rep_info := args[4];
    if IsBound( args[5] ) then 
	star_method := args[5];
    else
	star_method := recLieAlg.star_method;
    fi;

    if star_method = "simple" then 
	if rep_info = "vec" then 
	    basis := Basis( recLieAlg.L );
	    xx := LinearCombination( basis, x );
	    yy := LinearCombination( basis, y );
            wx := GUARANA.WeightOfLieAlgElm( recLieAlg, xx );
	    wy := GUARANA.WeightOfLieAlgElm( recLieAlg, yy );
	    max_w := recLieAlg.max_weight;
	    xx_star_yy := GUARANA.Star_Simple( xx, yy, wx, wy, max_w, 
	                                     "strucConst"   );
	    return Coefficients( basis, xx_star_yy );
	elif rep_info = "elm" then 
            wx := GUARANA.WeightOfLieAlgElm( recLieAlg, x );
	    wy := GUARANA.WeightOfLieAlgElm( recLieAlg, y );
	    max_w := recLieAlg.max_weight;
	    x_star_y := GUARANA.Star_Simple( x, y, wx, wy, max_w, 
	                                     "strucConst"   );
	    return x_star_y;
	else 
	    Error( "wrong rep_info" );
	fi;
    elif star_method = "symbolic" then
        Error( "Sorry, method not implemented yet" );
    else
	Error( "Wrong exp_method" );
    fi;
end;

#############################################################################
##
#E 
