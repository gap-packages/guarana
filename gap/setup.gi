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


BCH_TGroupRec := function( N )
    local sers,i, NN,l,leFactors,s,indices,class,
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
    
    # nilpotency class of group
    class := Length( leFactors );

    # weights of generators
    wei := 1;
    weights := [];
    for a in leFactors do
        for i in [1..a] do
            Add( weights, wei );
        od;
        wei := wei + 1;
    od;

    # get larges abelian group in the upper central series
    for i in [1..Length(sers)] do
        if IsAbelian( sers[i] ) then
            largestAbelian := i;
            break;
        fi;
    od;

    return rec( N := N, #basisN := basisN, 
                NN := NN,
                sers := sers, indices := indices,
                class := class, weights := weights,
                largestAbelian := largestAbelian );
end;

BCH_SetUpLieAlgebraRecordByMalcevbasis := function( recTGroup )
    local NN,hl,T,L;
    
    # get dimension of algebra
    NN := recTGroup.NN;
    hl := HirschLength( NN );

    # get prototype for structure constant table and Lie algebra 
    T:= EmptySCTable( hl, 0, "antisymmetric" );
    L:= LieAlgebraByStructureConstants( Rationals, T );

    

    # set the structure constants we already know,
    # i.e. lie brackets involving log g_i where g_i in Z(NN)
        # not necessary, because empty entries are interpreted 
        # as commuting elements

    return rec( L := L, recTGroup := recTGroup, scTable := T  );
end;


# x corresponds to 1
# y corresponds to 2 
BCH_EvaluateLieBracket := function( x, y, com, info )
    local r,l,tmp,i;
    tmp := [x,y];
    r := tmp[com[1]];

    l := Length( com );
    if info = "lieAlgebra" then
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

# g corresponds to 1
# h corresponds to 2 
BCH_EvaluateGroupCommutator := function( g, h, com )
    local tmp,r,l,i;
    tmp := [g,h];
    r := tmp[com[1]];

    l := Length( com );
    for i in [2..l] do
        r := Comm( r, tmp[com[i]] );
    od;
    return r;
end;


BCH_WeightOfCommutator := function( com, wx, wy )
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

BCH_CheckWeightOfCommutator := function( com, wx, wy, class )
    local w;
    w := BCH_WeightOfCommutator( com, wx, wy );
    if w > class then
        return false;
    else
        return true;
    fi;
end;

# Very SIMPLE Implementation !!
#
# compute x*y, where "*" is the BCH operation 
# wx ..... weight of x
# wy ..... weight of y
# class ...... nilpotency class of the lie algebra,
#              so brackets of Length class + 1 are always 0.
#
# info ..... strings which determines how lie bracket is evaluated
#
BCH_Star_Simple := function( recBCH, x, y, wx, wy, class, info  )
    local i,r,bchSers,com,a,term,max,min,bound;
    bchSers := recBCH.bchSers;
    
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
    # max + min* (bound-1 ) <= class
    bound := Int( (class-max)/min + 1 );

    # up to bound  compute the commutators and add them.
    # Note that the list contains commutators of length i at position i-1.
    for i in [1..bound-1] do
        for term in bchSers[i] do
            com := term[2];
            # check if weight of commutator is not to big
            if BCH_CheckWeightOfCommutator( com, wx, wy, class ) then
                a := BCH_EvaluateLieBracket( x, y, com, info );
                r := r + term[1]*a; 
            fi;
        od;
    od;
    
    return r;
end;


BCH_WeightOfLieAlgElm := function( recLieAlg, elm )
    local basisL,coeff,i,e;
    basisL := Basis( recLieAlg.L );
    coeff := Coefficients( basisL, elm );
    for i in [1..Length(coeff)] do
        e := coeff[i];
        if e <> 0 then
            return recLieAlg.recTGroup.weights[i];
        fi;
    od;
    return recLieAlg.recTGroup.class;
end;

BCH_AbstractLog_Simple_ByExponent := function( recLieAlg, recBCH,  exp )
    local  gensL,r,i,e,x,y,wx,wy,class;

    # setup
    gensL := GeneratorsOfAlgebra( recLieAlg.L );
    class := recLieAlg.recTGroup.class;

    r := Zero( recLieAlg.L );
    # go backwards through exponents. 
    # Which direction is more efficient ?
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
                wx := recLieAlg.recTGroup.weights[i];
                wy := BCH_WeightOfLieAlgElm( recLieAlg, y );
                r := BCH_Star_Simple( recBCH,x,y,wx,wy,class, "lieAlgebra" );
            fi;
        fi;
    od;
    return r;
end;

BCH_AbstractLogCoeff_Simple_ByExponent := function( recLieAlg, recBCH, exp )
    local r;
    r := BCH_AbstractLog_Simple_ByExponent( recLieAlg, recBCH,  exp );
    return Coefficients( Basis( recLieAlg.L ), r);
end;

BCH_AbstractLog_Simple_ByElm := function( recLieAlg, recBCH,  g )
    local exp,hl,l,expNN;
    
    # get exponets with respect to gens of T group NN
    # (g might an elment of PcpGroup which contains NN as a normal sugroup,
    # then Exponents(g) would be longer hl(NN))
    exp := Exponents( g );
    hl := HirschLength( recLieAlg.recTGroup.NN );
    l := Length( exp );
    expNN := exp{[l-hl+1..l]};
    return BCH_AbstractLog_Simple_ByExponent( recLieAlg, recBCH,  expNN );
end;

BCH_AbstractLogCoeff_Simple_ByElm := function( recLieAlg, recBCH,  g )
    local r;
    r := BCH_AbstractLog_Simple_ByElm( recLieAlg, recBCH,  g );
    return Coefficients( Basis( recLieAlg.L ), r);
end;

BCH_EvaluateLieBracketsInTermsOfLogarithmsSers 
                            := function( recBCH, recLieAlg, g,h,wg,wh )
    local bchLBITOL,r,max,min,bound,class,i,term,com,a,log_a;
    bchLBITOL := recBCH.bchLBITOL;
  
    r := Zero( recLieAlg.L );

    # compute upper bound for the Length of commutators, which 
    # can be involved
    class := recLieAlg.recTGroup.class;
    max := Maximum( wg,wh );
    min := Minimum( wg,wh );
    # max + min* (bound-1 ) <= class
    bound := Int( (class-max)/min + 1 );

    # up to bound  compute the commutators and add them.
    # Note that the list contains commutators of length i at position i-1.
    for i in [1..bound-1] do
        for term in bchLBITOL[i] do
            com := term[2];
            # check if weight of commutator is not to big
            if BCH_CheckWeightOfCommutator( com, wg, wh, class ) then
                # evaluate commutator in the group
                a := BCH_EvaluateGroupCommutator( g, h, com );
                # map to the Lie algebra
                log_a := BCH_AbstractLog_Simple_ByElm( recLieAlg, recBCH, a );
                r := r + term[1]*log_a; 
            fi;
        od;
    od;
    
    return r;

end;

# output is a list, as required for SetEntrySCTable
BCH_LieAlgElm2CoeffGenList := function( L, x )
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

#recBCH := BCH_FullBCHInformation( 5 );

#N := PcpGroupByMatGroup( PolExamples(4 ) );
#recTGroup := BCH_TGroupRec( N );  
#recLieAlg := BCH_SetUpLieAlgebraRecordByMalcevbasis( recTGroup );

# F := FreeGroup( 3 );
# N := NilpotentQuotient( F, 5 );
# recTGroup := BCH_TGroupRec( N );
# recLieAlg := BCH_SetUpLieAlgebraRecordByMalcevbasis( recTGroup );

BCH_ComputeStructureConstants := function( args )
    local factors, index_x, index_y, indices,indicesNoCenter,l,g,h,wg,wh,
          lie_elm,ll,gens,T,out,recLieAlg,recBCH;

    # setup
    recLieAlg := args[2];
    recBCH := args[1];
    indices := recLieAlg.recTGroup.indices;
    l := Length( indices );
    indicesNoCenter := indices{[1..l-1]};
    gens := GeneratorsOfGroup( recLieAlg.recTGroup.NN );
    T := recLieAlg.scTable;

    # go through blocks of generators backwards, 
    # where a block corresponds to the generators of the factors
    # of the upper central series
    for factors in Reversed( indicesNoCenter ) do
        for index_y in Reversed( factors ) do
            for index_x in  [1..index_y-1]   do
                g := gens[index_x];
                h := gens[index_y];
                wg := recLieAlg.recTGroup.weights[index_x];
                wh := recLieAlg.recTGroup.weights[index_y];
                # compute [log(g),log(h)]
                lie_elm := BCH_EvaluateLieBracketsInTermsOfLogarithmsSers( 
                                            recBCH, recLieAlg, g,h, wg, wh );
                # set entry in structure constant table 
                ll := BCH_LieAlgElm2CoeffGenList( recLieAlg.L, lie_elm );
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
