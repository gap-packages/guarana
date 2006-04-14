# Code for symbolic log and exp.
#      


# n ..... number of variables
# st ..... string for example  "x"
BCH_RationalVariableList := function( n, st )
     return List(
                  List( [1..n], i->Concatenation( st, String(i) ) ),
                  x->Indeterminate( Rationals, x : new ) );
end;

BCH_GenericElement := function( n, vars )
    return [[1..n], vars{[1..n]}];
end;

# domain is list of indezes we want
# 
BCH_GenericElementByDomain := function( vars, domain )
    return [ domain, vars{domain}];
end;

BCH_GenericBasisElement := function( i, vars )
    return [[i],[vars[i]]];    
end;


# list_x .... a list describing the lie algebra element x
#             alpha_i log(n_i) is represented as [ [i],[alpha_i] ]
#             alpha_i log(n_i) + alpha_j log(n_j) is represented as
#             [ [i,j],[alpha_i,alpha_j] ] and so on 
# 
BCH_Sum_Symbolic := function( list_x, list_y )
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


# list_x .... a list describing the lie algebra element x
#             alpha_i log(n_i) is represented as [ [i],[alpha_i] ]
#             alpha_i log(n_i) + alpha_j log(n_j) is represented as
#             [ [i,j],[alpha_i,alpha_j] ]
# 
BCH_EvaluateLieBracket_Symbolic := function( list_x, list_y, scTable )
    local index_x,index_y,coeff_x,coeff_y,prod,list_x1,list_x2,
          list_y1,list_y2,sum1,sum2,l;
    
    # catch trivial case 
    if Length( list_x[1] )=0 then
        return [[],[]];
    fi;
    if Length( list_y[1] )=0 then
        return [[],[]];
    fi;

    if Length( list_x[1]  ) = 1 and Length( list_y[1] )=1 then 
        index_x := list_x[1][1];
        index_y := list_y[1][1];
        coeff_x := list_x[2][1];
        coeff_y := list_y[2][1];
       
        prod := POL_CopyVectorList( scTable[index_x][index_y] );
        prod[2] := coeff_x*coeff_y*prod[2];

        return prod;
    elif Length( list_x[1] ) > 1 then 
        l := Length( list_x[1] );
        # split x
        list_x1 := [ list_x[1]{[1]} , list_x[2]{[1]} ];
        list_x2 := [ list_x[1]{[2..l]} , list_x[2]{[2..l]} ];
        
        sum1 := BCH_EvaluateLieBracket_Symbolic( list_x1, list_y, scTable );
        sum2 := BCH_EvaluateLieBracket_Symbolic( list_x2, list_y, scTable );

        return BCH_Sum_Symbolic( sum1, sum2 );

    elif Length( list_y[1] ) > 1 then
        l := Length( list_y[1] );
        # split y
        list_y1 := [ list_y[1]{[1]} , list_y[2]{[1]} ];
        list_y2 := [ list_y[1]{[2..l]} , list_y[2]{[2..l]} ];
        
        sum1 := BCH_EvaluateLieBracket_Symbolic( list_x, list_y1, scTable );
        sum2 := BCH_EvaluateLieBracket_Symbolic( list_x, list_y2, scTable );

        return BCH_Sum_Symbolic( sum1, sum2 );
    fi; 
    
    Error( "wrong input" );
end;

# longer lie bracktes
# x corresponds to 1
# y corresponds to 2 
BCH_EvaluateLongLieBracket_Symbolic := function( list_x, list_y, com, scTable )
    local r,l,tmp,i;
    tmp := [list_x,list_y];
    r := tmp[com[1]];

    l := Length( com );
    for i in [2..l] do
            r := BCH_EvaluateLieBracket_Symbolic( r, tmp[com[i]], scTable );
    od;
    return r;
end;


# list_x .... a list describing the lie algebra element x
#             alpha_i log(n_i) is represented as [ [i],[alpha_i] ]
#             alpha_i log(n_i) + alpha_j log(n_j) is represented as
#             [ [i,j],[alpha_i,alpha_j] ], etc...
# 
BCH_ComputeStarPolys 
                := function( recBCH, list_x, list_y, wx, wy, class, scTable )
    local i,r,bchSers,com,a,term,max,min,bound;
    bchSers := recBCH.bchSers;
    
    # start with terms which are not given by Lie brackets
    r := BCH_Sum_Symbolic( list_x, list_y );

    # trivial check 
    if Length( list_x[1] ) = 0  or Length( list_y[1] ) =  0 then
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
                a := BCH_EvaluateLongLieBracket_Symbolic( 
                                            list_x, list_y, com, scTable );
                r := BCH_Sum_Symbolic( r,[a[1],term[1]*a[2]] );
                #r := r + term[1]*a; 
            fi;
        od;
    od;
    
    return r;
end;

# For Log and Exp I just need 
# log x_i * sum_{j=i}^l beta_j log x_j
# So these are the polynomials I should save. It is NOT necessary to
# save x_i against a generic element of L.
# compute polynomials which can used to speed up star operation
#
# Example:
# exams_F2c := BCH_Get_FNG_TGroupRecords( 2, 9 );;
# recLieAlgs_bch_F2c := List( [2..Length( exams_F2c )], x-> BCH_LieAlgebraByTGroupRec( recBCH9,exams_F2c[x] ));;
# recLieAlg := recLieAlgs_bch_F2c[4];
#  BCH_AddStarPolynomialsToRecLieAlg( recLieAlg, recBCH9 );
#
BCH_AddStarPolynomialsToRecLieAlg := function( recLieAlg, recBCH )
    local i,n,vars_x,vars_y,x_i,elm_y,wx,wy,recStarPols,c,star_pols;

    # get variable for polynomials
    n := HirschLength( recLieAlg.recTGroup.NN );
    vars_x := BCH_RationalVariableList( n, "x" ); 
    vars_y := BCH_RationalVariableList( n, "y" );  

    # compute polynomials
    c := recLieAlg.recTGroup.class;
    star_pols := [];
    for i in [1..n] do 
        x_i := BCH_GenericBasisElement( i, vars_x );
        elm_y := BCH_GenericElementByDomain( vars_y, [i..n] );
        wx := recLieAlg.recTGroup.weights[i];
        wy := recLieAlg.recTGroup.weights[i];
        star_pols[i] := BCH_ComputeStarPolys( 
                 recBCH, x_i, elm_y,  wx, wy, c, recLieAlg.scTable );
    od;

    recStarPols := rec( vars_x := vars_x,
                        vars_y := vars_y,
                        pols := star_pols );
    recLieAlg.recStarPols := recStarPols; 
    return 0;                  
end;

# IN: x ........  [ [i], [x_i] ] 
#                 where x_i in general will be a variable.
#                 It can be also a polynomial, or some other ring element
#     y .......   [ [i,...,n], [ y_i,...,y_n]] 
#                 
# Star symbolically computed, for the input
# x_i * \sum_{j=i}^n \alpha_j y_y.
# For Log and Exp these are the kind of star operation which will 
# be needed. 
#
BCH_Star_Symbolic_SingleVersusGeneric := function( recLieAlg, x, y  )
    local n,vars_y,vars_x,indets,vals,pols,result,i,res,index_x,range_result;

    # simple test 
    if not x[1][1] = y[1][1] then Error( "wrong start index of y\n" ); fi;

    # setup
    n := HirschLength( recLieAlg.recTGroup.NN );
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


# For a given Lie algebra L(N) of dimension l, compute polynomials
# p_1,...,p_l such that 
# 
# Log( g_1^e_1 ... g_l^e_l ) = p_1( e )Log g_1 + ... + p_l( e ) Log g_l 
# 
# where (g_1,...,g_l) is a Malcev basis for N
# 
# Note: We could speed this up via using recStarPols
# 
BCH_ComputeSymbolicLogPolynomials := function( recLieAlg, recBCH )
    local n, vars_e,c,log_pols,tail,x_i,w_x_i,w_tail,i;

    # get variable for polynomials
    n := HirschLength( recLieAlg.recTGroup.NN );
    vars_e := BCH_RationalVariableList( n, "e" ); 

    # compute recursively polynomials as follows
    # log( g_1^e_1...g_n^e_n ) = log( g_1^e_1 )*(log( g_2^e_2...g_n^e_n ))
    #                          = e_1 log g_1 *  tail  
    c := recLieAlg.recTGroup.class;
    log_pols := [];
    tail := BCH_GenericBasisElement( n, vars_e );
    for i in Reversed( [1..(n-1)] ) do 
        x_i := BCH_GenericBasisElement( i, vars_e );
        w_x_i := recLieAlg.recTGroup.weights[i];
        w_tail := recLieAlg.recTGroup.weights[i+1];
        tail := BCH_ComputeStarPolys( recBCH, x_i, tail,  w_x_i, 
                                      w_tail, c, recLieAlg.scTable );
    od;

    return rec( pols := tail, vars_e := vars_e );
end;

# converts  [[i+1..n],[x_{i+1}..x_n]] to 
#           [[i..n],  [0,x_{i+1}..x_n]]
# 
BCH_Convert1 := function( x )
    local x_new,i,n;

    x_new := [];
    i := x[1][1] -1;
    n := x[1][Length(x[1])];
    x_new[1] := [i..n];
    x_new[2] := Concatenation( [ 0*x[2][1] ], x[2] );
    return x_new;
end;


BCH_ComputeSymbolicLogPolynomialsByStarPols := function( recLieAlg, recBCH )
    local n,vars_e,c,log_pols,tail,x_i,i;

    # get variable for polynomials
    n := HirschLength( recLieAlg.recTGroup.NN );
    vars_e := BCH_RationalVariableList( n, "e" ); 

    # compute recursively polynomials as follows
    # log( g_1^e_1...g_n^e_n ) = log( g_1^e_1 )*(log( g_2^e_2...g_n^e_n ))
    #                          = e_1 log g_1 *  tail  
    c := recLieAlg.recTGroup.class;
    tail := BCH_GenericBasisElement( n, vars_e );
    for i in Reversed( [1..(n-1)] ) do 
        x_i := BCH_GenericBasisElement( i, vars_e );
        tail := BCH_Convert1( tail );
        tail := BCH_Star_Symbolic_SingleVersusGeneric( recLieAlg, x_i, tail  );
    od;

    return rec( pols := tail, vars_e := vars_e );

end;

BCH_AddLogPolynomialsToLieAlgRecord := function( recLieAlg, recBCH )
    local pols, recLogPols;
    recLogPols := BCH_ComputeSymbolicLogPolynomialsByStarPols( 
                                                      recLieAlg, recBCH );
    recLieAlg.recLogPols := recLogPols;
    return 0;
end;




# IN: exp_n .... exponent vector of element n in \hat{N}
#    
# OUT: coeffcients of Log n, computed with the polynomials,
#      describing Log
#
# Example:
# exams_F3c := BCH_Get_FNG_TGroupRecords( 3, 5 );;
# recLieAlgs_bch_F3c := List( [2..Length( exams_F3c )], x-> BCH_LieAlgebraByTGroupRec( recBCH9,exams_F3c[x] ));;
#  BCH_AddStarPolynomialsToRecLieAlg( recLieAlgs_bch_F3c[5], recBCH9 );
#  BCH_AddLogPolynomialsToLieAlgRecord( recLieAlgs_bch_F3c[5], recBCH9 );
#  exp_n := BCH_Random_IntegralExpVector( recLieAlgs_bch_F3c[5], 2^12 );  
#  BCH_Logarithm_Symbolic( recLieAlgs_bch_F3c[5], exp_n );
#
#  compare 
#  BCH_AbstractLog_Simple_ByExponent( recLieAlgs_bch_F3c[5], recBCH9,exp_n);
#
# Example 2:
# 
# exams_unitr_2 := BCH_Get_Unitriangular_TGroupRecords( 10, 2 );
# recLieAlgs_bch_unitr_2 := List( [2..Length( exams_unitr_2)-3], x-> BCH_LieAlgebraByTGroupRec( recBCH9,exams_unitr_2[x] ));; 
#  BCH_AddStarPolynomialsToRecLieAlg( recLieAlgs_bch_unitr_2[5], recBCH9 );
#  BCH_AddLogPolynomialsToLieAlgRecord( recLieAlgs_bch_unitr_2[5], recBCH9 );
#  exp_n := BCH_Random_IntegralExpVector( recLieAlgs_bch_unitr_2[5], 2^10 );  
#  BCH_Logarithm_Symbolic( [recLieAlgs_bch_unitr_2[5], exp_n] );
# 
BCH_Logarithm_Symbolic := function( args )
    local n,indets,pols,coeffs,i,coeff,recLieAlg,exp_n;

    # setup    
    recLieAlg := args[1];
    exp_n := args[2];
    n := HirschLength( recLieAlg.recTGroup.NN ); 
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


BCH_ComputeSymbolicExpPolynomialsByStarPols := function( recLieAlg, recBCH )
    local n,vars_a,c,exp_pols,tail,a_bar,divider;

    # get variable for polynomials
    n := HirschLength( recLieAlg.recTGroup.NN );
    vars_a := BCH_RationalVariableList( n, "a" ); 

    

    # compute recursively polynomials as follows
    # exp( a_1 log g_1 + ... + a_n log g_n ) 
    # = exp( a_1 log g_1 ) * exp( -a_1 log_1 ) *
    #                               exp( a_1 log g_1 + ... + a_n log g_n )
    # = g_1^a_1 * exp( -a_1 log g_1 ) * exp( a_1 log g_1 + ... + a_n log g_n )
    # and then recurse
    c := recLieAlg.recTGroup.class;
    exp_pols := [];
    tail := BCH_GenericElement( n , vars_a  );
    for i in [1..n]  do 
        # get divider 
        a_bar := tail[2][1];
        divider := [ [i], [ (-1)* a_bar ] ];
        tail:= BCH_Star_Symbolic_SingleVersusGeneric(recLieAlg,divider,tail);
        Remove( tail[1], 1 );
        Remove( tail[2], 1 );
        Add( exp_pols, a_bar );
    od;
  
    return rec( pols := [[1..n],exp_pols], vars_a := vars_a );

end;

BCH_AddExpPolynomialsToLieAlgRecord := function( recLieAlg, recBCH )
    local pols, recExpPols;
    recExpPols := BCH_ComputeSymbolicExpPolynomialsByStarPols( 
                                                      recLieAlg, recBCH );
    recLieAlg.recExpPols := recExpPols;
    return 0;
end;


#  coeffs_x := BCH_Random_IntegralExpVector( recLieAlg, 2^12 );  
#  BCH_Exponential_Symbolic( recLieAlg, coeffs_x );
#  
#  compare 
#  BCH_Abstract_Exponential_ByVector( recBCH9, recLieAlg, coeffs_x );
BCH_Exponential_Symbolic := function( args )
    local n,indets,pols,exp,i,e;
    
    # setup
    recLieAlg := args[1];
    coeffs_x := args[2];
    n := HirschLength( recLieAlg.recTGroup.NN ); 
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

# Example usage:
# exams_unitr_2 := BCH_Get_Unitriangular_TGroupRecords( 10, 2 ); 
# recLieAlgs_bch_unitr_2 := List( [2..Length( exams_unitr_2 )], x-> BCH_LieAlgebraByTGroupRec( recBCH9,exams_unitr_2[x] ));;
# BCH_AddStarLogAndExpPols( [recLieAlgs_bch_unitr_2[5], recBCH9] );
# 
# exams_unitr_3 := BCH_Get_Unitriangular_TGroupRecords( 8, 3 ); 
# recLieAlgs_bch_unitr_3 := List( [2..Length( exams_unitr_3 )], x-> BCH_LieAlgebraByTGroupRec( recBCH9,exams_unitr_3[x] ));;
# BCH_AddStarLogAndExpPols( [recLieAlgs_bch_unitr_3[5], recBCH9] );
BCH_AddStarLogAndExpPols := function( args )
    local recLieAlg,recBCH;
    recLieAlg := args[1];
    recBCH := args[2];
    BCH_AddStarPolynomialsToRecLieAlg( recLieAlg, recBCH ); 
    BCH_AddLogPolynomialsToLieAlgRecord( recLieAlg, recBCH9 ); 
    BCH_AddExpPolynomialsToLieAlgRecord( recLieAlg, recBCH9 ); 
    return 0;
end;


# Next steps:
# - give Exp and Log an option, with which you can chosse the star operation
# - give collection an option
# - test the new collection and compare it with the old method


# test functions
BCH_Random_IntegralExpVector := function( recLieAlg, range )
    local n,ll,vec;
    n := HirschLength( recLieAlg.recTGroup.NN );
    ll := [ - range .. range ];
    vec := List( [ 1 .. n ], function ( x )
            return RandomList( ll );
        end );
    return vec;
end;

