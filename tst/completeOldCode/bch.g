#############################################################################
#
#                                                                 7.10.2005
# Code for the computation of coefficients of the BCH-series
# in k variables. That is, we compute an expression for z in
# 
# exp(y_1)* exp(y_2) *...* exp(y_k) = exp( z ) 
# 
# First we compute z as a polynomial in the x_i. For this purpose,
# we use the paper of Reinsch "A simple expression for the terms in the
# Baker-Campbell-Hausdorff series"
#
# These terms have now to be converted in expression containing only
# Lie brackets. It is not clear to me, how you can do this for k variables,
# but this is also not really needed. 
# What is is more important for us is the case
#
# Long Commutator in exp( y_1 ),exp(y_2)  = exp( z )
#
# For this we replace the y_i by y_1,-y_1,y_2,-y_2.
#
# The resulting polynomial in y_1,y_2 can be transformed into an 
# expression using Lie brackets in y_1, y_2 with the trick of Vaughan-Lee
#
# The aim of all this, is the computation of structure constants and
# an efficient log and exp. We use
#
# [logx,logy]= log( [x,y] ) + higher commutators in x,y           (*)
#
# Maybe it is possible to compute a series for this, i.e. we would
# get [logx,logy] = Q-linear combination of logs of group commutators.
# Approach: Compute this expression for a group of fixed nilpotency 
# class. Then we can probably show that a group of nilpotency class
# c+1 contains the old series plus some new terms. 
# We can compute this by replacing in (*) the higher commutators in x,y
# chi(x,y)_L by chi(x,y)_G + higher commutators in x,y, which we 
# compute with the extended BCH-series.
#  
# Maybe we can also find a formula for log applied to 
# g_1^e_1 * ... * g_k^e_k where k varies between 1 and 10. 
# What we need for this is exactly the Lie bracket expression for 
# z in exp(z) = exp(y_1) * ... exp( y_n )
#
# TODO: 
# - BCH series for any commutator in x,y
#   (also possible with Reinsch by using F,F^-1,G,G^-1)
# - BCH in k variable can be used for efficient log.
#   Can we deal with elements of the factors of the upper central
#   series simultaneously ? 
#   If you have a group of class 7 with 60 generators, it would 
#   be an overkill to use the BCH formula in 60 Variables!
# - Should we use basic commutators for all series ?

#############################################################################
#
# IN: v ..... upper or lower trinagular matrix with 0s on the diagonal
#
# OUT: exp(x)
#
BCH_Exponential := function( A, v )
    local n,exp, vPow, fac, i;
    n := Length( v[1] );
    exp := One( A )* IdentityMat(n) + v;
    vPow := v;
    fac := 1;
    for i in [2..n-1] do
        vPow := vPow*v;
        fac := fac*i;
        exp := exp + (fac^-1)*vPow;
    od;
    return exp;
end;

#############################################################################
#
# IN: x ..... upper or lower trinagular matrix with ones on the diagonal
#
# OUT: log(x)
#
BCH_Logarithm := function( x )
    local n,mat,matPow,log,i;
    n := Length( x[1] );
    mat := x - x^0;
    matPow := mat;
    log := mat;
    for i in [2..n-1] do
        matPow := matPow*mat;
        log := log + (-1)^(i-1)*(i)^-1*matPow;
    od;
    return log;
end;

#############################################################################
#
# matrix F as in paper of Reinsch 
#
BCH_Compute_F := function( n )
    local F,i,j;
    F := IdentityMat( n+1 );
    for i in [1..n+1] do
        for j in [i+1..n+1] do
            F[i][j] := 1/Factorial( j-i );
        od;
    od;
    return F;
end;

# compute F^pow where pow = +- 1
BCH_Compute_F_Power := function( n, pow )
    local mat,i,F_pow;
    mat := 0*IdentityMat( n+1 );
    for i in [1..n] do
        mat[i][i+1] := 1;
    od;
    F_pow := BCH_Exponential( Rationals, pow*mat ); 
    return F_pow;
end;

BCH_Compute_G := function( n, PolRing )
    local G,i,j,prod,I;
    I := IndeterminatesOfPolynomialRing( PolRing );
    G := IdentityMat( n+1 );
    for i in [1..n+1] do
        for j in [i+1..n+1] do
            prod := Product( I{[i..j-1]} );
            G[i][j] := 1/Factorial( j-i )* prod;
        od;
    od;
    return G;
end;

# compute G^pow where pow = +- 1
BCH_Compute_G_Power := function( n, PolRing, pow )
    local mat,i,G_pow,I;
    mat := 0*IdentityMat( n+1 );
    I := IndeterminatesOfPolynomialRing( PolRing );
    for i in [1..n] do
        mat[i][i+1] := I[i];
    od;
    G_pow := BCH_Exponential( Rationals, pow*mat ); 
    return G_pow;
end;



BCH_Compute_H := function( n, PolRing )
    local H,i,j,prod,I;
    I := IndeterminatesOfPolynomialRing( PolRing );
    H := IdentityMat( n+1 );
    for i in [1..n+1] do
        for j in [i+1..n+1] do
            prod := Product( I{[(i+n)..(j+n-1)]} );
            H[i][j] := 1/Factorial( j-i )* prod;
        od;
    od;
    return H;
end;

#############################################################################
#
# IN: mat ............ unipotent matrix
#     pow ............ rational number
# 
# OUT: mat^pow
BCH_RationalPowerOfUnipotentMat := function( mat, pow )
    local log;
    log := BCH_Logarithm( mat );
    return BCH_Exponential( Rationals, pow*log ); 
end;


#############################################################################
#
# Extension of the code of Reinsch. We do computations with k variables
#
# IN:
# n ...........  n-th term of the extended BCH series should be computed
# k ...........  k variables, i.e. exp(y_1)*...*exp(y_k)= exp(z)
# A ...........  Free associative algebra with one containing the y_i
#
# For simplification of the programming we replace the matrix F, which 
# originally does not contain varialbles with one similar to G, 
# i.e. we get a list of k matrices called H, whose i-th matrix
# contains only y_i
#
#A := FreeAssociativeAlgebraWithOne( Rationals, k, "x" );
#
BCH_Compute_HMatrices := function( n, k, A )
    local gens,i,j,mat,H;
    gens := GeneratorsOfAlgebraWithOne( A );
    H := [];
    for i in [1..k] do
        mat := Zero( A ) * IdentityMat( n+1 );
        for j in [1..n] do
            mat[j][j+1] := gens[i];
        od; 
        H[i] := BCH_Exponential( A, mat );
    od;
    return H;
end;


# typical input string = "1/180*x_1*x_2*x_5"
BCH_SigmaExpression2XY := function( n, string, A )
    local ll,k,gens,monomial,start,ll2,i,rat;

    ll := SplitString( string, "*" );
    k := Length( ll );

    #check if start with rational number
    rat := 1;
    start := 1;
    if ll[1][1] <> 'x' then;
        rat := Rat( ll[1] );
        start := 2;
    fi;

    ll2 := List( [start..k], x-> Int( ll[x]{[3]}) );

    # compute monomial
    gens := GeneratorsOfAlgebra( A ); 
    monomial := gens[1];
    for i in [1..n] do
        if i in ll2 then
            monomial := monomial*gens[3];
        else
            monomial := monomial*gens[2];
        fi;
    od;
    return rat*monomial;
end;

#############################################################################
#
# Operator T corresponding to Reinsch
# IN
# string ..........  polynomial in sigmaexpressions 
#                    typical input string = "1/180*x_1*x_2*x_5"
# OUT
# result of operator T.
# 
# Example:
# gap> BCH_SigmaExpression2List( 7, "1/180*x_1*x_2*x_5" );
# [ 1/180, [ 2, 2, 1, 1, 2, 1, 1 ] ]
#
BCH_SigmaExpression2List := function( n, string )
    local ll,k,start,ll2,i,rat,result;

    ll := SplitString( string, "*" );
    k := Length( ll );

    #check if start with rational number or minus sign
    rat := 1;
    start := 1;
    if ll[1][1] <> 'x' then
        if ll[1][1] = '-' and ll[1][2] = 'x' then
            rat := -1;
            ll[1] := ll[1]{[2..Length(ll[1])]};
        else
            rat := Rat( ll[1] );
            start := 2;
        fi;
    fi;

    # extract indices
    ll2 := List( [start..k], 
            x-> Int( ll[x]{[3..Length(ll[x])]}) );

    # compute monomial
    result := [];
    for i in [1..n] do
        if i in ll2 then
            # the 1 stands for y
            Add( result, 2);
        else
            Add( result, 1);
        fi;
    od;
    return [rat,result];
    #return    [result, rat ];
end;

#############################################################################
#
# IN:
# string .............. represents a MONOMIAL of a Free Algebra with One
#                       for example "(4/5)*y.2*y.1"
# OUT
# A list describing this element, in the  example [ 4/5, [ 2, 1 ] ]
#
BCH_Monomial2List := function( string )
    local ll,b,a,rat,inds,index,i,power,j,ss;
    # get rid of new lines commands
    ll := ReplacedString( string, "\n", "" );
    ll := SplitString( ll, "*" );
    b := Length( ll );

    # get coefficient
    a := Length( ll[1] );
    rat := Rat( ll[1]{[2..a-1]} );

    # get list which describes monomial
    inds := [];
    for i in [2..b] do
        # deal with powers
        power := 1;
        ss := SplitString( ll[i], "^" );
        if Length( ss ) > 1 then
            power := Int( ss[2] );
        fi;
        a := Length( ss[1] );
        index := Int( ss[1]{[3..a]} );
        for j in [1..power] do
            Add( inds, index );
        od;
    od;
    return [ rat, inds ];
end;


#############################################################################
#
# A := FreeAlgebraWithOne( Rationals, 2 );
# typical input string = "-1/720*x_1*x_2+1/180*x_3*x_5";
# Conversion corresponding to Reinsch
# Output lies in A
#
BCH_SumOfSigmaExpresions2XY := function( n, string, A )
    local ll1,ll2,ll3,k, gens,pol,i;
    ll1 :=  ReplacedString( string, "-", "+-" );
    ll2 := SplitString( ll1, "+" );
    ll3 := Filtered( ll2, x-> x <> "" );
    k := Length( ll3 );

    gens := GeneratorsOfAlgebra( A );
    pol := gens[1]*0;
    for i in [1..k] do
        pol := pol + BCH_SigmaExpression2XY( n, ll3[i], A );
    od;
    return pol;
end;

#############################################################################
#
# typical input string = "-1/720*x_1*x_2+1/180*x_3*x_5";
# Conversion corresponding to Reinsch
#
BCH_SumOfSigmaExpresions2List := function( n, string )
    local ll1,ll2,ll3,k,i,result;
    ll1 :=  ReplacedString( string, "-", "+-" );
    ll2 := SplitString( ll1, "+" );
    ll3 := Filtered( ll2, x-> x <> "" );
    k := Length( ll3 );

    result := [];
    for i in [1..k] do
        Add( result , BCH_SigmaExpression2List( n, ll3[i] ));
    od;
    return result;
end;

#############################################################################
#
# IN:
# string .............. represents a POLYNOMIAL of a Free Algebra with One
#                       for example "(1/2)*y.1*y.2+(-1/2)*y.2*y.1"
# OUT
# A list describing this element, in the  example 
# [ [ 1/2, [ 1, 2 ] ], [ -1/2, [ 2, 1 ] ] ]
#
BCH_Polynomial2List := function( string )
    local ll,y,result,i;
    ll := SplitString( string, "+" );
    y := Length( ll );
    result := [];
    for i in [1..y] do
        Add( result , BCH_Monomial2List( ll[i] ));
    od;
    return result;
end;

BCH_WeightInY := function( list )
    local weight,i;
    weight := 0;
    for i in list do
        if i = 2 then 
            weight := weight + 1;
        fi;
    od;
    return weight;
end;

#############################################################################
#
# 
# IN:
# ll ............ a list (for example [ -1/720, [ 1, 2 ] ] ) which
#                 specifies a monomial in x,y over Q
#                 Convention: 1 stands for x
#                             2 stands for y
# OUT:
# a list which specifies an elment of the Lie algebra, written 
# with Lie brackets, which is the image under the map of Vaughan-Lee
#
BCH_MonomialList2LieBracketList := function( ll )
    local k,result,i,l,sum,weight;
    k := Length( ll );
    result := [];
    for i in [1..k] do
        # take only terms starting with yx
        if ll[i][2][1]=2 and ll[i][2][2]=1 then
            l := StructuralCopy( ll[i] );
            # weight in y
            weight := BCH_WeightInY( l[2] );
            #sum := Sum( l[2] );
            l[1] := l[1]/weight;
            Add( result, l );
        fi;
    od;
    return result;
end;

#############################################################################
#
# 
# IN:
# ll ............ a list (for example [ -1/720, [ 1, 2 ] ] ) which
#                 specifies a monomial in y_1,...,y_k over Q
#                 Convention:  
#                 1 stands for  y_1
#                 2 stands for  y_2
#                 and so on
# OUT:
# a list which specifies an elment of the Lie algebra, written 
# with Lie brackets, which is the image under the map of Jacobson:
# d(m) = 1/i [m] 
# where "i" is the weigth of the monomial.
#
# TODO Is this function really correct for more than 2 variables ? 
#
BCH_MonimalInKVariablesList2LieBracketList := function( ll )
    local k,result,i,l,weight;
    k := Length( ll );
    result := [];
    for i in [1..k] do
        l := StructuralCopy( ll[i] );
        weight := Length( l[2] );
        l[1] := l[1]/weight;
        Add( result, l );
    od;
    return result;
end;


#############################################################################
#
# Calculation of the coefficients of the BCH series of weigt n.
# 
# Eample:
# 
# gap> BCH_ComputeCoefficientsAndBrackets( 2 );
# [ [ -1/2, [ 1, 0 ] ] ]
# gap> BCH_ComputeCoefficientsAndBrackets( 3 );
# [ [ -1/12, [ 1, 0, 1 ] ], [ 1/12, [ 1, 0, 0 ] ] ]
#
BCH_ComputeCoefficientsAndBrackets := function( n )
     local P,F,G,elm,string,ll,lie ;
     P := PolynomialRing( Rationals, n );
     F := BCH_Compute_F( n );;
     G := BCH_Compute_G( n, P );;
     elm := BCH_Logarithm( F*G )[1][n+1];
     string := String( elm );
     ll := BCH_SumOfSigmaExpresions2List( n, string );
     lie := BCH_MonomialList2LieBracketList( ll ); 
     return lie;
end;

# compute bch series up to weight c
BCH_Series := function( c )
    local sers, terms,i;
    sers := [];
    for i in [2..c] do
        terms := BCH_ComputeCoefficientsAndBrackets( i );
        sers := Concatenation( sers, terms );
    od;
    return sers;
end;

# compute bch series up to weight c
BCH_SeriesOrdered := function( c )
    local sers, terms,i;
    sers := [];
    for i in [2..c] do
        terms := BCH_ComputeCoefficientsAndBrackets( i );
        Add( sers, terms );
    od;
    return sers;
end;



BCH_ComputeMonomials := function( n )
     local P,F,G,elm,string,A,ll,lie ;
     P := PolynomialRing( Rationals, n );
     F := BCH_Compute_F( n );;
     G := BCH_Compute_G( n, P );;
     elm := BCH_Logarithm( F*G )[1][n+1];
     string := String( elm );
     #A := FreeAlgebraWithOne( Rationals, 2 );
     ll := BCH_SumOfSigmaExpresions2List( n, string );
     return ll;
end;

#############################################################################
#
# IN:
# n ......... weight, must be bigger than 1.
# word ...... list specifying word in e^y_1,e^y_2
#             For example [-1,-2,1,2] means y_1^-1 * y_2^-1* y_1 * y_2
#            
# Compute z where e^z = word( e^x, e^y ).
# Note that word can be any product of e^x,e^y,e^-x,e^-y
#
BCH_ComputeCoefficientsAndBracketWord := function( n, word )
     local P,F,G,elm,string,ll,lie,mats,i,prod,l ;
     P := PolynomialRing( Rationals, n );
     mats := [];
     # the following could be done better for longer words
     # by saving the F,F^-1,G,G^-1
     l := Length( word );
     for i in [1..l] do
         if AbsoluteValue( word[i] ) = 1 then 
            Add( mats, BCH_Compute_F_Power( n, SignInt( word[i] ) ));
         elif AbsoluteValue( word[i] ) = 2 then 
            Add( mats, BCH_Compute_G_Power( n, P, SignInt( word[i] ) ));
         else 
            Error( "Wrong Input\n" );
         fi;
     od;
     prod := Product( mats );
     elm := BCH_Logarithm( prod )[1][n+1];
     string := String( elm );
     ll := BCH_SumOfSigmaExpresions2List( n, string );
     lie := BCH_MonomialList2LieBracketList( ll ); 
     return lie;
end;

BCH_InverseWord := function( word )
    local res;
    res := Reversed( word );
    res := -1*res;
    return res;
end;

BCH_Comm2Word := function( comm )
    local l,res,word,wordInv;
    l := Length( comm );
    if l = 2 then 
         res := [ -comm[1],-comm[2],comm[1],comm[2] ];
    else  
         word := BCH_Comm2Word( comm{[1..l-1]} );
         wordInv := BCH_InverseWord( word );
         res := Concatenation( wordInv, [-comm[l]], word, [comm[l]] );
     fi;
     return res;
end;

#############################################################################
#
# IN:
# n ......... weight, must be bigger than 1.
# kappa ..... list specifyinga commutator in e^y_1,e^y_2
#             For example [1,2] means [y_1,y_2]
#            
# Compute z where e^z = kappa( e^x, e^y ).
#
BCH_ComputeCoefficientsAndBracketsByCommutator := function( n, kappa )
    local word;
    # check if 
    if Length( kappa ) = n then
        return [ [1,kappa] ];
    else
        word := BCH_Comm2Word( kappa );
        return BCH_ComputeCoefficientsAndBracketWord( n, word );
    fi;
end;


# BCH in k variables 
BCH_ComputeMonomialsExtendedBCHSeries := function( n, k )
     local A,H,mat,elm,string,ll;
     A := FreeAssociativeAlgebraWithOne( Rationals, k, "y" );     
     H := BCH_Compute_HMatrices( n, k, A );
     mat := Product( H );
     elm := BCH_Logarithm( mat )[1][n+1];
     string := StringPrint( elm );
     ll := BCH_Polynomial2List( string );
     return ll;
end;

# BCH in k variables
# e^z = e^y_1* ... * e^y_k
#
BCH_ComputeCoefficientsAndBracketsExtendedBCH :=  function( n, k )
     local A,H,mat,elm,string,ll,lie;
     A := FreeAssociativeAlgebraWithOne( Rationals, k, "y" );     
     H := BCH_Compute_HMatrices( n, k, A );
     mat := Product( H );
     elm := BCH_Logarithm( mat )[1][n+1];
     string := StringPrint( elm );
     ll := BCH_Polynomial2List( string );
     lie := BCH_MonimalInKVariablesList2LieBracketList( ll );
     return lie;
end;


#############################################################################
#
# Compute all terms of the BCH series up to a given weight wCom
# of all commutators up to weight wSers. 
# That is, we compute for all group commutators
# kappa of weight up to w2, all terms of weight up to  w1 in the expression
# z where 
# kappa( e^x, e^y )  = e^z
# 
# 
BCH_ComputeCommutatorSeries := function( wSers, wCom )
    local coms,i,a,lie,j,n,com,series,terms;

    if wSers < wCom then
        Error( "The BCH series of a commutator of Length wCom starts with\n",
               "Lie brackets of Length wCom.\n",
               "Therefore wSers should be at least as big as wCom." );
    fi;
 
    # produce a list coms, containing all possible commutators up 
    # to weight wCom
    coms := [];
    coms[1] := -1;
    coms[2] := [ [1,2], [2,1] ];
    for i in [3..wCom] do
        coms[i] := [];
        for a in coms[i-1] do 
            Add( coms[i], Concatenation( a, [1] ) );
            Add( coms[i], Concatenation( a, [2] ) );
        od;
    od;
    
    # for each commutator compute its BCH series up to weight wSers
    lie := [];
    lie[1] := -1;
    for i in [2..wCom] do
        lie[i] := [];
        for j in [1..Length(coms[i])] do 
            com := coms[i][j];
            series := [];
            # the BCH of com starts with weight length of com
            for n in [Length(com)..wSers] do
                terms := BCH_ComputeCoefficientsAndBracketsByCommutator( 
                                                                n, com );
                series := Concatenation( series, terms );
            od;
            Add( lie[i], series );
        od;
    od;
    
    return rec( coms := coms, lie := lie );
end;

#############################################################################
#
#
# kappa( logx, logy ) = log( kappa(x,y) ) +  tail
# factor*kappa( logx, logy ) = factor*log( kappa(x,y) ) + factor* tail
#
BCH_GetTail := function( kappa, sersOfComs, factor )
    local w, pos,sers,i,l;
    w := Length( kappa );
    pos := Position( sersOfComs.coms[w], kappa );
    
    sers := StructuralCopy( sersOfComs.lie[w][pos] );
    l := Length( sers );
    # get rid of starting term kappa( logx, logy )
    sers := sers{[2..l]};

    # change sign because we have to bring to the other side
    # and multiply it with factor
    for i in [1..l-1] do
        sers[i][1] := -factor*sers[i][1];
    od;
    return sers;
end;


# add lie brackets that are similar.
# for example [ [1/4, [ 2, 1, 2, 1 ] ], [1/4, [ 2, 1, 2, 1 ] ]]
# becomes [  [1/2, [ 2, 1, 2, 1 ] ] ]
BCH_SimplifySers := function( sers )
    local res,l,term,contained,term2,pos,i,res_without_zeros;
    res := [];
    l := Length( sers );
    if l = 0 then
        return [];
    fi;
    Add( res, sers[1] );
    for i in [2..l] do
        term := sers[i];
        # test if term[2] is already contained in res
        contained := false;
        for pos in [1..Length(res)] do
            term2 := res[pos];
            if term2[2] = term[2] then
                contained := true;
                break;
            fi;
        od;
        if contained then
            # summarize
            res[pos][1] := res[pos][1] + term[1];           
        else            
            # else add term
            Add( res, term );
        fi;
    od;

    # delete entries with zero coeff
    res_without_zeros := [];
    for pos in [1..Length( res )] do    
        # if zero coeff then add
        if res[pos][1] <> 0 then
            Add( res_without_zeros, res[pos] );
        fi;
    od;
       
    return res_without_zeros;
end;


#############################################################################
#
# Aim:
#[logx,logy] = Q-linear combination of logs of group commutators.
#
#
# c .............. class of nilpotent groups for which the result
#                  of the computation should be used. 
#                  That is, c is the maximal length of a group commutator
#                  which is non-trivial.
#
# Approach: Compute this expression for a group of fixed nilpotency 
# class. Then we can probably show that a group of nilpotency class
# c+1 contains the old series plus some new terms. 
# We can compute this by replacing in (*) the higher commutators in x,y
# chi(x,y)_L by chi(x,y)_G + higher commutators in x,y, which we 
# compute with the extended BCH-series.
#
#
BCH_LieBracketInTermsOfLogs := function( c )
    local sersOfComs,logs,term,lieTerms,tail,i,logs_ordered,le;

    # get bch series of all commutators up to weight c
    sersOfComs := BCH_ComputeCommutatorSeries( c, c );

    # we know the first term of the logs expressions
    logs := [ [1,[1,2]] ];
    
    # get remaining terms in logx, logy
    lieTerms := BCH_GetTail( [1,2], sersOfComs, 1 );

    # replace iteratively all remaining Lie bracktes
    i := 1;
    while i <= Length( lieTerms ) do
        term := lieTerms[i];
        # replace it by factor*log( kappa(x,y) ) + factor*tail
        Add( logs, term );
        tail := BCH_GetTail( term[2], sersOfComs, term[1] );
        Append( lieTerms, tail );
        i := i+1;
    od;
    logs := BCH_SimplifySers( logs );

    # order result corresponding to weight
    logs_ordered := [];
    for i in [2..c] do
        logs_ordered[i-1] := [];
    od;
    for term in logs do
        le := Length( term[2] );
        Add( logs_ordered[le-1], term );
    od;
  
    return logs_ordered;
end;

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

BCH_FullBCHInformation := function( w  )
    local bchSers, bchLBITOL;
    bchSers := BCH_SeriesOrdered( w );
    bchLBITOL :=  BCH_LieBracketInTermsOfLogs( w );

    return rec( bchSers := bchSers, 
                bchLBITOL := bchLBITOL );
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
