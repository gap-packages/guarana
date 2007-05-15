#############################################################################
##
#W test.gi               GUARANA package                     Bjoern Assmann
##
## Methods for testing the setup of the Malcev correspondence.
##
##
#H  @(#)$Id$
##
#Y 2006
##

#############################################################################
#
# Test Log, Exp and the computation of structur constants
#
InstallOtherMethod( Random, 
               "for Malcev objects (Guarana)", 
               true, 
               [IsMalcevObjectRep], 
               0,
function( malcevObject)
    local n, range, exps, a, x;
    n := malcevObject!.dim;
    range := 10;
    exps := List( [1..n], x-> Random( [-range..range] ) );
    a := MalcevGenElementByExponents( malcevObject, exps );
    # compute also the lie element
    x := LieElement( a );
    return a;
end);

InstallOtherMethod( Random, 
               "for Malcev objects (Guarana)", 
               true, 
               [IsMalcevObjectRep, IsInt], 
               0,
function( malcevObject, range)
    local n, exps, a, x;
    n := malcevObject!.dim;
    exps := List( [1..n], x-> Random( [-range..range] ) );
    a := MalcevGenElementByExponents( malcevObject, exps );
    # compute also the lie element
    x := LieElement( a );
    return a;
end);

InstallMethod( RandomGrpElm, 
               "for Malcev objects (Guarana)", 
               true, 
               [IsMalcevObjectRep], 
               0,
function( malcevObject)
    local n, range, exps;
    n := malcevObject!.dim;
    range := 10;
    exps := List( [1..n], x-> Random( [-10..10] ) );
    return MalcevGrpElementByExponents( malcevObject, exps );
end);

InstallOtherMethod( RandomGrpElm, 
               "for Malcev objects and integers (Guarana)", 
               true, 
               [IsMalcevObjectRep, IsInt], 
               0,
function( malcevObject, range)
    local n, exps;
    n := malcevObject!.dim;
    exps := List( [1..n], x-> Random( [-range..range] ) );
    return MalcevGrpElementByExponents( malcevObject, exps );
end);

InstallMethod( RandomLieElm, 
               "for Malcev objects (Guarana)", 
               true, 
               [IsMalcevObjectRep], 
               0,
function( malcevObject)
    local n, range, exps;
    n := malcevObject!.dim;
    range := 10;
    exps := List( [1..n], x-> Random( [-10..10] ) );
    return MalcevLieElementByCoefficients( malcevObject, exps );
end);

InstallOtherMethod( RandomLieElm, 
               "for Malcev objects and integers (Guarana)", 
               true, 
               [IsMalcevObjectRep, IsInt], 
               0,
function( malcevObject, range)
    local n, exps;
    n := malcevObject!.dim;
    exps := List( [1..n], x-> Random( [-range..range] ) );
    return MalcevLieElementByCoefficients( malcevObject, exps );
end);

GUARANA.MO_Test_ExpOfLog := function(  malcevObject, noTests,range, method )
    local g, x, gg, i;
    for i in [1..noTests] do
	g := RandomGrpElm( malcevObject, range );
	SetLogAndExpMethod( malcevObject, method );
	#SetExpMethod( malcevObject, method );
	x := Log( g );
	gg := Exp( x );
	if not gg = g then 
	    Error( "Mist \n" );
	fi;
    od;
    return 0;
end;

if false then 
    noTests := 100;
    range := 2^4;
    malObjs := GUARANA.Get_FNG_MalcevObjects( 2, 9 );
    for malObj in malObjs do
	Print( "Testing ", malObj ,"\n" );
	GUARANA.MO_Test_ExpOfLog( malObj, noTests, range, "simple" );
    od;

    # big example
    GUARANA.MO_Test_ExpOfLog( malObj, 100, 2^10, "pols" );
fi;

#############################################################################
##
## old functions

GUARANA.Test_LogOfExp := function(  recLieAlg, noTests )
    local x,exp,x2,i;
    for i in [1..noTests] do
        x := Random( recLieAlg.L );
        exp := GUARANA.Abstract_Exponential_ByElm(  recLieAlg, x );
        x2 := GUARANA.AbstractLog_Simple_ByExponent( recLieAlg,  exp );
        if not x = x2 then
            Error( "Mist\n" );
        fi;
    od;
    return 0;
end;
# Testing some examples 
if false then
    recLieAlgs := GUARANA.Get_FNG_LieAlgRecords( 2, 5 );

    recLieAlg := recLieAlgs[5];
    GUARANA.Test_LogOfExp( recLieAlg, 10 );
    recLieAlg.malcevBasisInfo := "gen";
    GUARANA.Test_ExpOfLog( recLieAlg, 10 );

    recLieAlgs := GUARANA.Get_Unitriangular_LieAlgRecords( 6, 2 );
    for i in [1..Length( recLieAlgs )] do
	recLieAlg := recLieAlgs[i];
	GUARANA.Test_LogOfExp( recLieAlg, 10 );
	recLieAlg.malcevBasisInfo := "gen";
	GUARANA.Test_LogOfExp( recLieAlg, 10 );
    od;

fi;

GUARANA.Test_ExpOfLog := function(  recLieAlg, noTests,range )
    local i,dim,domain,exp,x,exp2;
    dim := recLieAlg.dim;
    domain := [-range..range];
    for i in [1..noTests] do
        exp := List( [1..dim], x -> Random( domain ) );
        x := GUARANA.AbstractLog_Simple_ByExponent( recLieAlg,  exp );
        exp2 := GUARANA.Abstract_Exponential_ByElm(  recLieAlg, x );
        if not exp = exp2 then
            Error( "Mist\n" );
        fi;
    od;
    return 0;
end;


# Testing some examples 
if false then
    recLieAlgs := GUARANA.Get_FNG_LieAlgRecords( 2, 5 );

    recLieAlg := recLieAlgs[5];
    GUARANA.Test_ExpOfLog( recLieAlg, 10, 100 );
    recLieAlg.malcevBasisInfo := "gen";
    GUARANA.Test_ExpOfLog( recLieAlg, 10, 100 );

    recLieAlgs := GUARANA.Get_Unitriangular_LieAlgRecords( 6, 2 );
    for i in [1..Length( recLieAlgs )] do
	recLieAlg := recLieAlgs[i];
	GUARANA.Test_ExpOfLog( recLieAlg, 10, 100 );
	recLieAlg.malcevBasisInfo := "gen";
	GUARANA.Test_ExpOfLog( recLieAlg, 10, 100 );
    od;

fi;

# produce random element of Tr_0(dim,Q)
GUARANA.RandomNilpotentMat := function( dim )
    local range,ll,kk,g,j,k;
    range := 7;
    ll := [ - range .. range ];
    kk := [ - range .. range ];
    ll := Filtered( ll, function ( x )
            return x <> 0;
        end );
    kk := Filtered( kk, function ( x )
            return x <> 0;
        end );

    g := NullMat( dim, dim, Rationals );
    for j  in [ 1 .. dim ]  do
        for k  in [ j .. dim ]  do
            if j < k  then
                g[j][k] := RandomList( ll ) / RandomList( kk );
            else
                g[j][k] := 0;
            fi;
        od;
    od;
    return g;
end;

# n ....... dim of matrices that are used.
GUARANA.TestBchSeriesOrdered := function(  n )
    local no_tests,i,x,y,wx,wy,x_star_y,exp_x_star_y,exp_x,exp_y;

    no_tests := 100;
    for i in [1..no_tests] do
        # produce two random matrices x,y in Tr_0(n,Q)
        x := GUARANA.RandomNilpotentMat( n );
        y := GUARANA.RandomNilpotentMat( n );
        wx := 1;
        wy := 1;

        # compute  exp(x*y) with BCH, 
        # note that we need terms of length at most n-1
        x_star_y := GUARANA.Star_Simple (  x, y, wx, wy, n-1, "matrix"  );
        exp_x_star_y := GUARANA.Exponential( Rationals, x_star_y );   

        # compute exp(x)exp(y) and compare
        exp_x := GUARANA.Exponential( Rationals, x );
        exp_y := GUARANA.Exponential( Rationals, y );

        if not exp_x_star_y = exp_x*exp_y then
            Error( "Mist \n " );
        fi;

    od;
    return 0;
end;

#############################################################################
##
#F GUARANA.Test_ComputeCommutatorSeries( recComSers, n, kappa )
##
## IN 
## recComSers ...............as computed by GUARANA.ComputeCommutatorSeries
## n ....................... dimension of matrices that are going to be used
## kappa ................... commutator given as list. for example
##                           [1,2,1,1]
##
## Tests if the computation of kappa(x,y), where x,y are two random
## matrices in Tr_0(n,Q), using a BCH like formula is correct.
##
## Example
## compute all terms of the commutator bch series up to weights 6 (wSers)
## of all commutators up to weight 3 (wCom)
## Note that n should be bigger wSers.
## recComSers := GUARANA.ComputeCommutatorSeries( 6, 3 );
## GUARANA.Test_ComputeCommutatorSeries( recComSers, 4, [1,2,1] );
##
GUARANA.Test_ComputeCommutatorSeries := function( recComSers, n, kappa )
    local  no_tests,i,x,y,wx,wy,exp_x,exp_y,exp_z,r,max_weight,
          sers,com,a,term,exp_z_bch,weight, pos;

    no_tests := 10;
    for i in [1..no_tests] do
        # produce two random matrices x,y in Tr_0(n,Q)
        x := GUARANA.RandomNilpotentMat( n );
        y := GUARANA.RandomNilpotentMat( n );
        wx := 1;
        wy := 1;
   
        # compute kappa( exp(x),exp(y) ) normally
        exp_x := GUARANA.Exponential( Rationals, x );
        exp_y := GUARANA.Exponential( Rationals, y );
        exp_z := GUARANA.EvaluateGroupCommutator( exp_x, exp_y, kappa );

        # compute z where exp(z)= kappa( exp(x),exp(y) ) with extended BCH
        r := NullMat( n,n );
        max_weight := n-1;
        weight := Length( kappa );
        pos := Position( recComSers.coms[weight], kappa );
        sers := recComSers.lie[weight][pos];
        for term in sers do
            com := term[2];
            #Print( "term ", term, "\n" );
            # check if weight of commutator is not to big
            if GUARANA.CheckWeightOfCommutator( com, wx, wy, max_weight ) then
                # evaluate commutator in the lie algebra
                a := GUARANA.EvaluateLieBracket( x, y, com, "matrix" );
                r := r + term[1]*a;
                #        Print( "log_a ", log_a, "\n" );
                #        Print( "r ", r, "\n\n" );
             fi;
        od; 
        exp_z_bch := GUARANA.Exponential( Rationals, r );

        # compare
        if not exp_z_bch = exp_z then 
            Error( "Mist\n" );
        fi; 
    od;
    return 0;
end;

GUARANA.TestSeriesLieBracketInTermsOfLogs := function(  n )
    local no_tests,i,x,y,x_bracket_y,g,h,wg,wh,r,max_weight,max,min,bound,j,
          term,a,log_a,x_bracket_y_2,bchLBITOL,com;

    no_tests := 10;
    for j in [1..no_tests] do
        # produce two random matrices x,y in Tr_0(n,Q)
        x := GUARANA.RandomNilpotentMat( n );
        y := GUARANA.RandomNilpotentMat( n );
              
        # compute [x,y] in normal way
        x_bracket_y := LieBracket( x, y );

        # compute [x,y] with series "liebrackets in terms of logs"
            g := GUARANA.Exponential( Rationals, x );
            h := GUARANA.Exponential( Rationals, y );
            wg := 1;
            wh := 1;

            bchLBITOL := GUARANA.recBCH.bchLBITOL;
            r := NullMat( n,n );
            # compute upper bound for the Length of commutators, which 
            # can be involved
            max_weight := n-1;
            max := Maximum( wg,wh );
            min := Minimum( wg,wh );
            # max + min* (bound-1 ) <= max_weight
            bound := Int( (max_weight-max)/min + 1 );

            # up to bound  compute the commutators and add them.
            # Note that the list contains comms of length i at position i-1.
            for i in [1..bound-1] do
                for term in bchLBITOL[i] do
                    com := term[2];
                    #Print( "term ", term, "\n" );
                    # check if weight of commutator is not to big
                    if GUARANA.CheckWeightOfCommutator( com, wg, wh,        
			max_weight ) then
                        # evaluate commutator in the group
                        a := GUARANA.EvaluateGroupCommutator( g, h, com );
                        # map to the Lie algebra
                        log_a := GUARANA.Logarithm( a );
                        r := r + term[1]*log_a;
                        #Print( "log_a ", log_a, "\n" );
                        #Print( "r ", r, "\n\n" );
                    fi;
                od;
            od;
            x_bracket_y_2 := r;
    
        # compare
        if not x_bracket_y = x_bracket_y_2 then
            Error( "Mist\n" );
        fi;
    od;
    return 0;
end;

#############################################################################
##
#E
