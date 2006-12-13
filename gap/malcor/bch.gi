#############################################################################
##
#W bch.gi               Guarana package                        Bjoern Assmann
##
## Code for the computation of coefficients of the BCH-series
## in k variables. That is, we compute an expression for z in
## 
## exp(y_1)* exp(y_2) *...* exp(y_k) = exp( z ) 
## 
## First we compute z as a polynomial in the x_i. For this purpose,
## we use the paper of Reinsch 'A simple expression for the terms in the
## Baker-Campbell-Hausdorff series'
## Then we transform it to an expression only involving Liebrackets  
## with the Dynkin bracket operator in  Michael Vaughan-Lee's book.
##
## In addition, this file contains code for the computation of 
## identities that are related to BCH formula. 
## For example we compute an identity that expresses Lie brackets
## in terms of a linear combination of logarithms of group 
## commutators.  
##
## TODO/Questions: 
## - BCH in k variable can be used for efficient log.
##   Can we deal with elements of the factors of the upper central
##   series simultaneously ? 
##   If you have a group of class 7 with 60 generators, it would 
##   be an overkill to use the BCH formula in 60 Variables!
## - Should we use basic commutators for all series ?

##############################################################################
##
#F GUARANA.Exponential( A, v )                                                
##
## IN: v ..... upper or lower trinagular matrix with 0s on the diagonal
##
## OUT: exp(x) where exp is given by the usual power series. 
##
GUARANA.Exponential := function( A, v )
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
##
#F GUARANA.Logarithm( x )
##
## IN: x ..... upper or lower trinagular matrix with ones on the diagonal
##
## OUT: log(x) where log is given by the usual power series. 
##
GUARANA.Logarithm := function( x )
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
# matrix F,G and H as in paper of Reinsch 
#
GUARANA.Compute_F := function( n )
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
GUARANA.Compute_F_Power := function( n, pow )
    local mat,i,F_pow;
    mat := 0*IdentityMat( n+1 );
    for i in [1..n] do
        mat[i][i+1] := 1;
    od;
    F_pow := GUARANA.Exponential( Rationals, pow*mat ); 
    return F_pow;
end;

GUARANA.Compute_G := function( n, PolRing )
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
GUARANA.Compute_G_Power := function( n, PolRing, pow )
    local mat,i,G_pow,I;
    mat := 0*IdentityMat( n+1 );
    I := IndeterminatesOfPolynomialRing( PolRing );
    for i in [1..n] do
        mat[i][i+1] := I[i];
    od;
    G_pow := GUARANA.Exponential( Rationals, pow*mat ); 
    return G_pow;
end;

GUARANA.Compute_H := function( n, PolRing )
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
GUARANA.RationalPowerOfUnipotentMat := function( mat, pow )
    local log;
    log := GUARANA.Logarithm( mat );
    return GUARANA.Exponential( Rationals, pow*log ); 
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
# originally does not contain any variables, with one similar to G, 
# i.e. we get a list of k matrices called H, whose i-th matrix
# contains only y_i
#
# A := FreeAssociativeAlgebraWithOne( Rationals, k, "x" );
#
GUARANA.Compute_HMatrices := function( n, k, A )
    local gens,i,j,mat,H;
    gens := GeneratorsOfAlgebraWithOne( A );
    H := [];
    for i in [1..k] do
        mat := Zero( A ) * IdentityMat( n+1 );
        for j in [1..n] do
            mat[j][j+1] := gens[i];
        od; 
        H[i] := GUARANA.Exponential( A, mat );
    od;
    return H;
end;

# typical input string = "1/180*x_1*x_2*x_5"
GUARANA.SigmaExpression2XY := function( n, string, A )
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
# gap> GUARANA.SigmaExpression2List( 7, "1/180*x_1*x_2*x_5" );
# [ 1/180, [ 2, 2, 1, 1, 2, 1, 1 ] ]
#
GUARANA.SigmaExpression2List := function( n, string )
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
GUARANA.Monomial2List := function( string )
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
GUARANA.SumOfSigmaExpresions2XY := function( n, string, A )
    local ll1,ll2,ll3,k, gens,pol,i;
    ll1 :=  ReplacedString( string, "-", "+-" );
    ll2 := SplitString( ll1, "+" );
    ll3 := Filtered( ll2, x-> x <> "" );
    k := Length( ll3 );

    gens := GeneratorsOfAlgebra( A );
    pol := gens[1]*0;
    for i in [1..k] do
        pol := pol + GUARANA.SigmaExpression2XY( n, ll3[i], A );
    od;
    return pol;
end;

#############################################################################
#
# typical input string = "-1/720*x_1*x_2+1/180*x_3*x_5";
# Conversion corresponding to Reinsch
#
GUARANA.SumOfSigmaExpresions2List := function( n, string )
    local ll1,ll2,ll3,k,i,result;
    ll1 :=  ReplacedString( string, "-", "+-" );
    ll2 := SplitString( ll1, "+" );
    ll3 := Filtered( ll2, x-> x <> "" );
    k := Length( ll3 );

    result := [];
    for i in [1..k] do
        Add( result , GUARANA.SigmaExpression2List( n, ll3[i] ));
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
GUARANA.Polynomial2List := function( string )
    local ll,y,result,i;
    ll := SplitString( string, "+" );
    y := Length( ll );
    result := [];
    for i in [1..y] do
        Add( result , GUARANA.Monomial2List( ll[i] ));
    od;
    return result;
end;

GUARANA.WeightInY := function( list )
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
GUARANA.MonomialList2LieBracketList := function( ll )
    local k,result,i,l,sum,weight;
    k := Length( ll );
    result := [];
    for i in [1..k] do
        # take only terms starting with yx
        if ll[i][2][1]=2 and ll[i][2][2]=1 then
            l := StructuralCopy( ll[i] );
            # weight in y
            weight := GUARANA.WeightInY( l[2] );
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
GUARANA.MonimalInKVariablesList2LieBracketList := function( ll )
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
# gap> GUARANA.ComputeCoefficientsAndBrackets( 2 );
# [ [ -1/2, [ 1, 0 ] ] ]
# gap> GUARANA.ComputeCoefficientsAndBrackets( 3 );
# [ [ -1/12, [ 1, 0, 1 ] ], [ 1/12, [ 1, 0, 0 ] ] ]
#
GUARANA.ComputeCoefficientsAndBrackets := function( n )
    local indeces, indets, P, F, G, elm, string, ll, lie, i;
     # create polynomial ring with indeterminates called x_...
     indeces := [1..n];
     indets := [];
     for i in indeces do 
	 Add( indets, Concatenation( "x_", String( i ) ));
     od;
     P := PolynomialRing( Rationals, indets );
     F := GUARANA.Compute_F( n );;
     G := GUARANA.Compute_G( n, P );;
     elm := GUARANA.Logarithm( F*G )[1][n+1];
     string := String( elm );
     ll := GUARANA.SumOfSigmaExpresions2List( n, string );
     lie := GUARANA.MonomialList2LieBracketList( ll ); 
     return lie;
end;

# compute bch series up to weight c
GUARANA.BchSeries := function( c )
    local sers, terms,i;
    sers := [];
    for i in [2..c] do
        terms := GUARANA.ComputeCoefficientsAndBrackets( i );
        sers := Concatenation( sers, terms );
    od;
    return sers;
end;

# compute bch series up to weight c
GUARANA.BchSeriesOrdered := function( c )
    local sers, terms,i;
    sers := [];
    for i in [2..c] do
        terms := GUARANA.ComputeCoefficientsAndBrackets( i );
        Add( sers, terms );
    od;
    return sers;
end;


GUARANA.ComputeMonomials := function( n )
     local P,F,G,elm,string,A,ll,lie ;
     P := PolynomialRing( Rationals, n );
     F := GUARANA.Compute_F( n );;
     G := GUARANA.Compute_G( n, P );;
     elm := GUARANA.Logarithm( F*G )[1][n+1];
     string := String( elm );
     #A := FreeAlgebraWithOne( Rationals, 2 );
     ll := GUARANA.SumOfSigmaExpresions2List( n, string );
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
GUARANA.ComputeCoefficientsAndBracketWord := function( n, word )
     local P,F,G,elm,string,ll,lie,mats,i,prod,l,indeces,indets;

     # create polynomial ring with indeterminates called x_...
     indeces := [1..n];
     indets := [];
     for i in indeces do 
	 Add( indets, Concatenation( "x_", String( i ) ));
     od;
     P := PolynomialRing( Rationals, indets );

     mats := [];
     # the following could be done better for longer words
     # by saving the F,F^-1,G,G^-1
     l := Length( word );
     for i in [1..l] do
         if AbsoluteValue( word[i] ) = 1 then 
            Add( mats, GUARANA.Compute_F_Power( n, SignInt( word[i] ) ));
         elif AbsoluteValue( word[i] ) = 2 then 
            Add( mats, GUARANA.Compute_G_Power( n, P, SignInt( word[i] ) ));
         else 
            Error( "Wrong Input\n" );
         fi;
     od;
     prod := Product( mats );
     elm := GUARANA.Logarithm( prod )[1][n+1];
     string := String( elm );
     ll := GUARANA.SumOfSigmaExpresions2List( n, string );
     lie := GUARANA.MonomialList2LieBracketList( ll ); 
     return lie;
end;

GUARANA.InverseWord := function( word )
    local res;
    res := Reversed( word );
    res := -1*res;
    return res;
end;

GUARANA.Comm2Word := function( comm )
    local l,res,word,wordInv;
    l := Length( comm );
    if l = 2 then 
         res := [ -comm[1],-comm[2],comm[1],comm[2] ];
    else  
         word := GUARANA.Comm2Word( comm{[1..l-1]} );
         wordInv := GUARANA.InverseWord( word );
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
GUARANA.ComputeCoefficientsAndBracketsByCommutator := function( n, kappa )
    local word;
    # check if 
    if Length( kappa ) = n then
        return [ [1,kappa] ];
    else
        word := GUARANA.Comm2Word( kappa );
        return GUARANA.ComputeCoefficientsAndBracketWord( n, word );
    fi;
end;


# BCH in k variables 
GUARANA.ComputeMonomialsExtendedBCHSeries := function( n, k )
     local A,H,mat,elm,string,ll;
     A := FreeAssociativeAlgebraWithOne( Rationals, k, "y" );     
     H := GUARANA.Compute_HMatrices( n, k, A );
     mat := Product( H );
     elm := GUARANA.Logarithm( mat )[1][n+1];
     string := StringPrint( elm );
     ll := GUARANA.Polynomial2List( string );
     return ll;
end;

# BCH in k variables
# e^z = e^y_1* ... * e^y_k
#
GUARANA.ComputeCoefficientsAndBracketsExtendedBCH :=  function( n, k )
     local A,H,mat,elm,string,ll,lie;
     A := FreeAssociativeAlgebraWithOne( Rationals, k, "y" );     
     H := GUARANA.Compute_HMatrices( n, k, A );
     mat := Product( H );
     elm := GUARANA.Logarithm( mat )[1][n+1];
     string := StringPrint( elm );
     ll := GUARANA.Polynomial2List( string );
     lie := GUARANA.MonimalInKVariablesList2LieBracketList( ll );
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
GUARANA.ComputeCommutatorSeries := function( wSers, wCom )
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
                terms := GUARANA.ComputeCoefficientsAndBracketsByCommutator( 
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
GUARANA.GetTail := function( kappa, sersOfComs, factor )
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
GUARANA.SimplifySers := function( sers )
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
# [logx,logy] = Q-linear combination of logs of group commutators.
#
#
# c .............. class of nilpotent groups for which the result
#                  of the computation should be used. 
#                  That is, c is the maximal length of a group commutator
#                  which is non-trivial.
#
# Approach: Compute this expression for a group of fixed nilpotency 
# class. A group of nilpotency class
# c+1 contains the old series plus some new terms. 
# We can compute this by replacing in (*) the higher commutators in x,y
# chi(x,y)_L by chi(x,y)_G + higher commutators in x,y, which we 
# compute with the extended BCH-series.
#
# Example:
# gap> GUARANA.LieBracketInTermsOfLogs( 3 );
# [ [ [ -1, [ 2, 1 ] ] ], [ [ 1/2, [ 2, 1, 2 ] ], [ 1/2, [ 2, 1, 1 ] ] ] ]
#
#
GUARANA.LieBracketInTermsOfLogs := function( c )
    local sersOfComs,logs,term,lieTerms,tail,i,logs_ordered,le;

    # get bch series of all commutators up to weight c
    sersOfComs := GUARANA.ComputeCommutatorSeries( c, c );

    # we know the first term of the logs expressions
    logs := [ [-1, [2,1]] ];
    
    # get remaining terms in logx, logy
    lieTerms := GUARANA.GetTail( [2,1], sersOfComs, -1 );

    # replace iteratively all remaining Lie bracktes
    i := 1;
    while i <= Length( lieTerms ) do
        term := lieTerms[i];
        # replace it by factor*log( kappa(x,y) ) + factor*tail
        Add( logs, term );
        tail := GUARANA.GetTail( term[2], sersOfComs, term[1] );
        Append( lieTerms, tail );
        i := i+1;
    od;
    logs := GUARANA.SimplifySers( logs );

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

############################################################################
##
#F ComputeBchIdentities
##
## IN 
## w ................... maximal length of Lie bracket
## 
## OUT
## 
## A record containg
## length ............ maximal length of Lie bracket. 
## bchSers ........... a list containing all terms of the Bch formula
##		       up to length w
## bchLBITOL....... .. a list containing all terms up to length w of 
##                     the identity that expresses a lie bracket 
##                     in terms of logs of group commutators. 
## 
## Example: 
## gap> GUARANA.ComputeBchIndentities( 3 );
## rec(
##   length := 3,
##   bchSers := [ [ [ -1/2, [ 2, 1 ] ] ], [ [ -1/12, [ 2, 1, 2 ] ], 
##                  [ 1/12, [ 2,1, 1 ] ] ] ],
##   bchLBITOL := [ [ [ 1, [ 1, 2 ] ] ], [ [ 1/2, [ 2, 1, 2 ] ],
##	            [ 1/2, [ 2, 1, 1 ] ] ] ] )
##				   
## 
GUARANA.ComputeBchIndentities := function( w )
    local bchSers, bchLBITOL;
    bchSers := GUARANA.BchSeriesOrdered( w );
    bchLBITOL :=  GUARANA.LieBracketInTermsOfLogs( w );

    return rec( length := w,
                bchSers := bchSers,
                bchLBITOL := bchLBITOL );
end;

