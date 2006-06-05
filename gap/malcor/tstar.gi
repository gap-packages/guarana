#############################################################################
##
#W tstar.gi               GUARANA.package                     Bjoern Assmann
##
## Methods for computing the BCH star operation using a tree structure.
##
#H  @(#)$Id$
##
#Y 2006
##

#############################################################################
##
## IN
## com  ..................... list containing 1 and 2's. For example
##                            [2,1,2,1]
##                            com has to star with 2.
## OUT
## Number such that com is a "binary" representation of it. 
##
GUARANA.TS_Com2Num := function( com )
    local n, num, i;

    n := Length( com );
    num := 0;
    for i in com do 
	num := num*2;
	if i = 2 then 
	    num := num + 1;
	fi;
    od;
    return num;
end;

#############################################################################
##
#M Comm( x, y, com, tree ) ...............long lie or group commutator 
#M Comm( x, y, com )                      for Malcev elements
#M LieBracket( x,y, com , tree )
#M LieBracket( x,y, com )
##
## COMMENT
## com  ..................... list containing 1 and 2's. For example
##                            [2,1,2,1]
##                            com has to start with 2.
##
## com(x,y) = [ com_prev(x,y), x ] or [ com_prev(x,y,), y ]
## 
InstallOtherMethod( Comm, 
               "for Malcev elments (Guarana)",
	       true,
	        [IsMalcevElement, IsMalcevElement, IsList, IsList ],
		0, 
function( x, y, com, tree )
    local num_com, n, com_prev, com_prev_xy, tmp, com_xy;

    # check if com alreay known
    num_com := GUARANA.TS_Com2Num( com );
    if IsBound( tree[num_com] ) then 
	return tree[num_com];
    fi;

    # check if weight is not to high
    if not GUARANA.CheckWeightOfMalcevCommutator( x, y, com ) then
	return 0*x;
    fi;

    # get/compute com_prev(x,y)
    n := Length( com );
    com_prev := com{[1..n-1]};
    com_prev_xy := Comm(x, y, com_prev, tree );

    # compute com(x,y)
    tmp := [x,y];
    com_xy := Comm( com_prev_xy, tmp[com[n]] );

    # update tree 
    tree[num_com] := com_xy;

    return com_xy;
end);

InstallOtherMethod( LieBracket, 
               "for Malcev lie elments (Guarana)",
	       true,
	        [IsMalcevLieElement, IsMalcevLieElement, IsList, IsList ],
		0, 
function( x, y, com, tree )
    return Comm( x,y, com, tree );
end);

InstallOtherMethod( Comm, 
               "for Malcev elments (Guarana)",
	       true,
	        [IsMalcevElement, IsMalcevElement, IsList ],
		0, 
function( x, y, com )
    local tree, com_new, i;

    if com[1]=2 then 
	tree := [y];
	return Comm( x,y, com, tree );
    elif com[1] = 1 then 
	# switch the roles of x,y
	tree := [x];
	com_new := [];
	for i in com do 
	    if i = 1 then 
		Add( com_new, 2 );
	    else 
		Add( com_new, 1 );
	    fi;
	od;
	return Comm( y,x, com_new, tree );
    fi;
end);

InstallOtherMethod( LieBracket, 
               "for Malcev elments (Guarana)",
	       true,
	        [IsMalcevLieElement, IsMalcevLieElement, IsList ],
		0, 
function( x, y, com )
    return Comm( x,y,com );
end);

#############################################################################
##
#M BCHStar .......................................... for Malcev lie elements
##
GUARANA.MO_BCHStar_Simple := function( x, y )
    local wx, wy, malObj, max_weight, bchSers, r, tree, max, 
          min, bound, com, a, i, term;

    # setup
    wx := Weight( x );
    wy := Weight( y );
    malObj := x!.malcevObject;
    max_weight := malObj!.max_weight;
    bchSers := GUARANA.recBCH.bchSers;
    
    # start with terms which are not given by Lie brackets
    r := x + y;

    # trivial check 
    if x = 0*x or y = 0*y then
        return r;
    fi;

    # set up tree
    tree := [y];

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
            a := LieBracket( x, y, com, tree );
            r := r + term[1]*a; 
        od;
    od;
    return r;
end;

InstallMethod( BCHStar, 
               "for Malcev Lie elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevLieElement, IsMalcevLieElement ],
		0, 
function( x, y )
    local malcevObject;
    malcevObject := x!.malcevObject;
    if StarMethod( malcevObject ) = "pols" then 
        return GUARANA.MO_Star_Symbolic( x,y );
    else
        return GUARANA.MO_BCHStar_Simple( x, y );
    fi;
end);

InstallMethod( BCHStar, 
               "for Malcev Gen elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevGenElement, IsMalcevGenElement ],
		0, 
function( g, h )
    local l_g, l_h, l_res;
    l_g := LieElement( g );
    l_h := LieElement( h );
    l_res := BCHStar( l_g, l_h );
    return MalcevGenElementByLieElement( l_res );
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

fi;

#############################################################################
## 
## Old code
##
GUARANA.TS_EvaluateShortLieBracket := function( x,y, info )
    local res;

    if info = "strucConst" then
	res := x*y;
    elif info = "matrix" then 
        res := LieBracket( x, y );
    else 
        Error( "wront input \n " );
    fi;
    return res;
end;

## COMMENT
## com(x,y) = [ com_prev(x,y), x ] or [ com_prev(x,y,), y ]
##
GUARANA.TS_EvaluateLieBracket := function( x, y, com, tree, info )
    local num_com, n, com_prev, com_prev_xy, tmp, com_xy;

    # check if com alreay known
    num_com := GUARANA.TS_Com2Num( com );
    if IsBound( tree[num_com] ) then 
	return tree[num_com];
    fi;

    # get/compute com_prev(x,y)
    n := Length( com );
    com_prev := com{[1..n-1]};
    com_prev_xy := GUARANA.TS_EvaluateLieBracket(x,y,com_prev, tree, info );

    # compute com(x,y)
    tmp := [x,y];
    com_xy := GUARANA.TS_EvaluateShortLieBracket( com_prev_xy,tmp[com[n]],
                                                  info );

    # update tree 
    tree[num_com] := com_xy;

    return com_xy;
end;

## IN
## x,y ................... elements of a Lie algebra.
##                         Brackets of Length max_weight + 1 are always 0.
## wx,wy ................. their weights. 
## max_weight .............max_weight of the basis elms of the Lie algebra. 
## info .................  strings which determines how lie brackets
##                         are evaluated.
##                         Possibilities
##                         "strucConst"  ( Lie algebra given by 
##                                         structure constant )
##                         "matrix"      ( Lie algebra given by 
##                                         matrices )
##
## OUT 
## x*y, where "*" is the BCH operation 
##
GUARANA.TS_Star := function( args  )
    local x, y, wx, wy, max_weight, info, bchSers, r, tree, max, 
          min, bound, com, a, i, term;

    x := args[1];
    y := args[2];
    wx := args[3];
    wy := args[4];
    max_weight := args[5];
    info := args[6];

    bchSers := GUARANA.recBCH.bchSers;
    
    # start with terms which are not given by Lie brackets
    r := x + y;

    # trivial check 
    if x = 0*x or y = 0*y then
        return r;
    fi;

    # set up tree
    tree := [y];

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
                a := GUARANA.TS_EvaluateLieBracket( x, y, com, tree, info );
                r := r + term[1]*a; 
            fi;
        od;
    od;
    return r;
end;

if false then 
lieRecs := GUARANA.Get_FNG_LieAlgRecords( 2, 9  );;
lieRec := lieRecs[9];;
x := Random( lieRec.L );;
y := Random( lieRec.L );;
x_s_y := GUARANA.Star_Simple( x,y, 1,1, lieRec.max_weight, "strucConst" );
x_s_y := GUARANA.TS_Star( [x, y, 1, 1, lieRec.max_weight, "strucConst"]  );
fi;
   
# profiling

if false then
lieRecs := GUARANA.Get_FNG_LieAlgRecords( 2, 9  );;
lieRec := lieRecs[9];;
x := Random( lieRec.L );;
y := Random( lieRec.L );;
    input :=  [x, y, 1, 1, lieRec.max_weight, "strucConst"];; 
    func := GUARANA.TS_Star;
    subfuncs := [ GUARANA.TS_EvaluateShortLieBracket,
		  GUARANA.TS_EvaluateLieBracket, 
		  GUARANA.TS_Star];
    GUARANA.Profile( func, input, subfuncs );
fi;

##  gap>     GUARANA.Profile( func, input, subfuncs );
##    count  self/ms  chld/ms  function
##    217     1331        0  TS_EvaluateShortLieBracket
##    395        3     1331  GUARANA.TS_EvaluateLieBracket
##      1       36     1334  GUARANA.TS_Star
##           1370           TOTAL
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
GUARANA.TS_ComputeStarPolys 
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
end;
#############################################################################
##
#E
