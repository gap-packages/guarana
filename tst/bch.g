
BCH_Abstract_Exponential_ByElm := function( recBCH, recLieAlg, x )
    local indices,basis,class,tail,coeffs,largestAbelian,exp_x,i,factor,
          divider,w_divider,w_tail,l,exp_x_2ndPart,j;

    indices := recLieAlg.recTGroup.indices;
    basis := Basis( recLieAlg.L );
    class := recLieAlg.recTGroup.class;

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
            w_divider := BCH_WeightOfLieAlgElm ( recLieAlg, divider );
            w_tail := BCH_WeightOfLieAlgElm ( recLieAlg, tail );
            tail := BCH_Star_Simple( recBCH, divider, tail, w_divider, 
                                     w_tail, class, "lieAlgebra"  );
        
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

BCH_Abstract_Exponential_ByVector := function( recBCH, recLieAlg, vec )
    local basis,x;
    basis := Basis( recLieAlg.L );
    x := LinearCombination( basis, vec );
    return BCH_Abstract_Exponential_ByElm( recBCH, recLieAlg, x ); 
end;

BCH_LieAlgebraByTGroupRec := function( recBCH, recTGroup )
    local recLieAlg;

    recLieAlg := BCH_SetUpLieAlgebraRecordByMalcevbasis( recTGroup );
    BCH_ComputeStructureConstants( [recBCH, recLieAlg] );
   
    return recLieAlg;
end;

#############################################################################
#
# TEST Functions
#
#############################################################################


BCH_Test_LogOfExp := function( recBCH, recLieAlg, noTests )
    local x,exp,x2,i;

    for i in [1..noTests] do
        x := Random( recLieAlg.L );
        exp := BCH_Abstract_Exponential_ByElm( recBCH, recLieAlg, x );
        x2 := BCH_AbstractLog_Simple_ByExponent( recLieAlg, recBCH, exp );
        if not x = x2 then
            Error( "Mist\n" );
        fi;
    od;
    return 0;
end;

BCH_Test_ExpOfLog := function( recBCH, recLieAlg, noTests,range )
    local i,hl,domain,exp,x,exp2;
    hl := HirschLength( recLieAlg.recTGroup.NN );
    domain := [-range..range];
    for i in [1..noTests] do
        exp := List( [1..hl], x -> Random( domain ) );
        x := BCH_AbstractLog_Simple_ByExponent( recLieAlg, recBCH, exp );
        exp2 := BCH_Abstract_Exponential_ByElm( recBCH, recLieAlg, x );
        if not exp = exp2 then
            Error( "Mist\n" );
        fi;
    od;
    return 0;
end;


BCH_Test := function( n )
    local P,F,G,elm,string,A,ll,lie,kk,kk2;        
 
    P := PolynomialRing( Rationals, n );
    F := BCH_Compute_F( n );;
    G := BCH_Compute_G( n, P );;
    elm := BCH_Logarithm( F*G )[1][n+1];
    string := String( elm );
    A := FreeAlgebraWithOne( Rationals, 2 );
    ll := BCH_SumOfSigmaExpresions2List( n, string );
    lie := BCH_MonomialList2LieBracketList( ll );
    Length( ll );
    Length( lie );
    kk := List( lie, x->x[2] );
    kk2 := AsSet( kk );
    return [ Length( lie ), Length( AsSet( kk ) )];
end;


BCH_List2Number := function( base, ll )
    local n,l,numb,i;
    n := Length( ll );
    #throw away the first two entries 1 and 0
    l := ll{[3..n]};

    numb := 0;
    for i in Reversed( [3..n]) do
        if ll[i]<>0 then
            numb := numb + base^(n-i);
        fi;
    od;
    return numb;
end;

# produce random element of Tr_0(dim,Q)
BCH_RandomNilpotentMat := function( dim )
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

BCH_TestBchSeriesOrdered := function( recBCH, n )
    local no_tests,i,x,y,wx,wy,x_star_y,exp_x_star_y,exp_x,exp_y;

    no_tests := 100;
    for i in [1..no_tests] do
        # produce two random matrices x,y in Tr_0(n,Q)
        x := BCH_RandomNilpotentMat( n );
        y := BCH_RandomNilpotentMat( n );
        wx := 1;
        wy := 1;

        # compute  exp(x*y) with BCH, 
        # note that we need terms of length at most n-1
        x_star_y := BCH_Star_Simple ( recBCH, x, y, wx, wy, n-1, "matrix"  );
        exp_x_star_y := BCH_Exponential( Rationals, x_star_y );   

        # compute exp(x)exp(y) and compare
        exp_x := BCH_Exponential( Rationals, x );
        exp_y := BCH_Exponential( Rationals, y );

        if not exp_x_star_y = exp_x*exp_y then
            Error( "Mist \n " );
        fi;

    od;
    return 0;
end;

# recComSers as computed by BCH_ComputeCommutatorSeries
# recComSers := BCH_ComputeCommutatorSeries( 6, 3 );
BCH_Test_ComputeCommutatorSeries := function( recComSers, n, kappa )
    local  no_tests,i,x,y,wx,wy,exp_x,exp_y,exp_z,r,class,
          sers,com,a,term,exp_z_bch,weight, pos;

    no_tests := 10;
    for i in [1..no_tests] do
        # produce two random matrices x,y in Tr_0(n,Q)
        x := BCH_RandomNilpotentMat( n );
        y := BCH_RandomNilpotentMat( n );
        wx := 1;
        wy := 1;
   
        # compute kappa( exp(x),exp(y) ) normally
        exp_x := BCH_Exponential( Rationals, x );
        exp_y := BCH_Exponential( Rationals, y );
        exp_z := BCH_EvaluateGroupCommutator( exp_x, exp_y, kappa );

        # compute z where exp(z)= kappa( exp(x),exp(y) ) with extended BCH
        r := NullMat( n,n );
        class := n-1;
        weight := Length( kappa );
        pos := Position( recComSers.coms[weight], kappa );
        sers := recComSers.lie[weight][pos];
        for term in sers do
            com := term[2];
            #Print( "term ", term, "\n" );
            # check if weight of commutator is not to big
            if BCH_CheckWeightOfCommutator( com, wx, wy, class ) then
                # evaluate commutator in the lie algebra
                a := BCH_EvaluateLieBracket( x, y, com, "matrix" );
                r := r + term[1]*a;
                #        Print( "log_a ", log_a, "\n" );
                #        Print( "r ", r, "\n\n" );
            fi;
        od;
        exp_z_bch := BCH_Exponential( Rationals, r );

        # compare
        if not exp_z_bch = exp_z then 
            Error( "Mist\n" );
        fi;
    od;
    return 0;
end;

BCH_TestSeriesLieBracketInTermsOfLogs := function( recBCH, n )
    local no_tests,i,x,y,x_bracket_y,g,h,wg,wh,r,class,max,min,bound,j,
          term,a,log_a,x_bracket_y_2,bchLBITOL,com;

    no_tests := 10;
    for j in [1..no_tests] do
        # produce two random matrices x,y in Tr_0(n,Q)
        x := BCH_RandomNilpotentMat( n );
        y := BCH_RandomNilpotentMat( n );
              
        # compute [x,y] in normal way
        x_bracket_y := LieBracket( x, y );

        # compute [x,y] with series "liebrackets in terms of logs"
            g := BCH_Exponential( Rationals, x );
            h := BCH_Exponential( Rationals, y );
            wg := 1;
            wh := 1;

            bchLBITOL := recBCH.bchLBITOL;
            r := NullMat( n,n );
            # compute upper bound for the Length of commutators, which 
            # can be involved
            class := n-1;
            max := Maximum( wg,wh );
            min := Minimum( wg,wh );
            # max + min* (bound-1 ) <= class
            bound := Int( (class-max)/min + 1 );

            # up to bound  compute the commutators and add them.
            # Note that the list contains comms of length i at position i-1.
            for i in [1..bound-1] do
                for term in bchLBITOL[i] do
                    com := term[2];
                    #Print( "term ", term, "\n" );
                    # check if weight of commutator is not to big
                    if BCH_CheckWeightOfCommutator( com, wg, wh, class ) then
                        # evaluate commutator in the group
                        a := BCH_EvaluateGroupCommutator( g, h, com );
                        # map to the Lie algebra
                        log_a := BCH_Logarithm( a );
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
