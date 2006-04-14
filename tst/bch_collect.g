#
# Code for collection in polycyclic groups using the BCH - Formula
#

#############################################################################
#
# IN:  G.......... pc group
#      N.......... normal subgroup of G, 
#                  N is a T-group
# 
# OUT  record which contains
#      GG   ...........  group isomorphic to G and which has 
#                        a presentation which goes through the 
#                        upper central series of N and G/N
#      NN  ............. group isomorphic to N which has a presntation
#                        which goes through the upper central series of N
#      G,N ............. as input
#      sersN ........... upper central series of N
#
BCH_NewUnderlyingCollectorParentGrp := function( G, N )
    local sersN, sersG, GG,NN;
  
    # compute series which goes through the factors of the upper central
    # series and G/N
    sersN := UpperCentralSeries( N );
    sersG := [ G ];
    Append( sersG, sersN );

    # compute pcp group on generators which goes through sersG
    GG := PcpGroupBySeries( sersG ); 

    # get image of N in GG
    NN := PreImage( GG!.bijection, N );  

    return rec( GG := GG, G := G, NN := NN, N :=N,  sersN := sersN );
end;

#############################################################################
#
# IN: R .........   R contains output as computed by 
#                   BCH_NewUnderlyingCollectorParentGrp
#
# OUT record which contains the datastructure which is needed
#     to compute the Lie algebra corresponding to R.N via the 
#     BCH method.
#
BCH_TGroupRec_ByUpperCentralSeries := function( R )
    local sers,i, NN, l,leFactors,s,indices,class,
          wei, weights,a, largestAbelian;

    sers := R.sersN;
    NN := R.NN;
 
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

    return rec( N := R.N , NN := NN,
                sers := sers, indices := indices,
                class := class, weights := weights,
                largestAbelian := largestAbelian );  
end;

#############################################################################
#
# IN: g ........ element of parent group GG (new presentation of G )
#     NN ....... new presentation of T-group N
#     liealg ... liealgebra record of N
# OUT: Automorphism of L(NN) corresponding to the conjugation action
#      of g
#
BCH_LieAutMatrix := function( g, NN, recLieAlg, recBCH )
    local AutMat,pcp,exps,n,i,coeff;
    
    AutMat := [];
    pcp := Pcp( NN );
    exps := LinearActionOnPcp( [g], pcp )[1];
    n := Length( exps );
    for i in [1..n] do
        # compute coefficients of log n_i^g
        coeff := BCH_AbstractLogCoeff_Simple_ByExponent( 
                                            recLieAlg, recBCH, exps[i] ); 
        Add( AutMat, coeff );
    od;
    return AutMat;
   
end;


#############################################################################
#
# IN: G ..... pc-group
#     N ..... T-group which is normal in G
#
# OUT: record containing
#      all datastructures which are needed to do fast
#      conjugation.
#
BCH_FastConjugationRec := function( G, N, recBCH )
    local recNewParent, recTGroup,recLieAlg,pcpGG,n,l,m,factorGens,lieAuts;

    # change underlying generating set of G
    recNewParent := BCH_NewUnderlyingCollectorParentGrp( G, N );
    
    # set up lie algebra and compute structure constants 
    recTGroup := BCH_TGroupRec_ByUpperCentralSeries( recNewParent );
    recLieAlg := BCH_LieAlgebraByTGroupRec( recBCH, recTGroup ); 


    pcpGG := Pcp( recNewParent.GG );
    n := HirschLength( N );
    l := Length( pcpGG );
    m := l-n;
    factorGens := pcpGG{[1..m]};

    # compute corresponding autom. of L(NN)
    lieAuts := List( pcpGG, x -> BCH_LieAutMatrix( 
                            x, recNewParent.NN, recLieAlg, recBCH ));

    return rec( recNewParent := recNewParent,
                recTGroup := recTGroup,
                recLieAlg := recLieAlg,
                factorGens := factorGens,
                lieAuts := lieAuts,
                recBCH := recBCH,
                info := "bch" );
end;



#############################################################################
#
# IN: FCR ....... fast multiplication record
#     n ......... element of T-group N
#     c ......... element of parent group of N
#                 c = w( factorGens )
# OUT: n^c computed via the BCH lie algebra
#
BCH_FastConjugation := function( FCR, n, c, recBCH )
    local exp,logn,lieAut,l,expl;

    exp := Exponents( c );

    # comput log n
    logn := BCH_AbstractLogCoeff_Simple_ByElm( FCR.recLieAlg, recBCH, n ); 

    #get corresponding automorphism of L(N)
    lieAut := MappedVector( exp, FCR.lieAuts );

    #apply autom. to log(n)
    l := logn*lieAut;

    #compute exponents of exponential  of l 
    expl :=  BCH_Abstract_Exponential_ByVector ( recBCH, FCR.recLieAlg, l );

    return expl;
end;


#############################################################################
#
# IN: FCR ............. fast conjugation record
#     log1_coeff....... element of lie algebra represented as coefficient 
#                       vector
#     exp_n2.......     exponent of element n2 in T-group
#     recBCH
#
# OUT:
# normal form of  exp( log1) n_2 
# computed in an efficient way, by sucessively dividing 
# off and applying conjugation
#
# Use that in the group we have
# n_1^x_1...n_l^x_l n_1^y_1 = n_1^x_1+y_1 (n_2^x_2 ... n_l^x_l)^n_1^y_1
#
# Let log1 = sum alpha_i log n_i = log( n^x ) 
# Then for the computation of exp( log1*log2 ) we do the following:
# Divide alpa_1 log n_1 = x_1 log n_1 up and apply to the rest 
# the automorphism corresponding to the conjugation with n_1^y_1.
# Save x_1+y_1 as the first exponent of exp( log1*log2 ).
# Continue to divide up and to apply the conjugation automorphisms. 
# 
BCH_ExpOfStarProduct := function( FCR, log1_coeff, exp_n2, recBCH )
    local indices,basis,class,l,log1,log2,tail,coeffs,largestAbelian,
          exp_r,i,j,divider,w_divider,w_tail,aut,le,exp_r_2ndPart,factor;

    # setup 
    indices := FCR.recLieAlg.recTGroup.indices;
    basis := Basis( FCR.recLieAlg.L );
    class := FCR.recLieAlg.recTGroup.class;
    l := Length( FCR.factorGens );

     # get element of lie algebra
    log1 := LinearCombination( basis, log1_coeff ); 

    # divide up and apply conjugation trick untill remaining elment lies in an
    # abelian group.
    tail := log1;
    coeffs := StructuralCopy( log1_coeff );
    largestAbelian := FCR.recLieAlg.recTGroup.largestAbelian;
    exp_r := [];
    #for i in [1..largestAbelian] do
    for i in [1..largestAbelian-1] do
        factor := indices[i];
        for j in factor do 
            # get element to divide of
            divider := -coeffs[j]*basis[j];

            # save exponents of divider plus y_j
            Add( exp_r, coeffs[j]+exp_n2[j] );

            # divide off
            w_divider := BCH_WeightOfLieAlgElm ( FCR.recLieAlg, divider );
            w_tail := BCH_WeightOfLieAlgElm ( FCR.recLieAlg, tail );
            tail := BCH_Star_Simple( recBCH, divider, tail, w_divider, 
                                     w_tail, class, "lieAlgebra"  );
        
            # compute matrix corresponding to automorph. and apply it
            ##BCH_RationalPowerOfUnipotentMat is very slow !
            ##aut := BCH_RationalPowerOfUnipotentMat( 
            ##                        FCR.lieAuts[l+j], exp_n2[j] );
            coeffs := Coefficients( basis, tail );
            #aut := FCR.lieAuts[l+j]^exp_n2[j];
            #coeffs := coeffs*aut;
            coeffs := coeffs*FCR.lieAuts[l+j]^exp_n2[j];
            tail := LinearCombination( basis, coeffs );
        od;
    od;

    # test intermediate result
    le := Length( exp_r );
    if not coeffs{[1..le]} = 0 *  coeffs{[1..le]} then
        Error( "Failure in  BCH_ExpOfStarProduct\n" );
    fi;

    # get the remaining coefficients 
    exp_r_2ndPart := coeffs{[le+1..Length(coeffs)]} + 
                                    exp_n2{[le+1..Length(coeffs)]};
    
    return Concatenation( exp_r, exp_r_2ndPart );
end;

#############################################################################
#
# IN: FCR ....... fast multiplication record of a group which is
#                 the semidiret product of an abelian and a T-Group
#                 so G = A semi N, for group elements g = a(g)n(g)
#
#                 Attention: we assume just power relations of the form
#                 g_i^r_i = 1
#
#     g1,g2 ..... random element of G
#
#                 c = w( factorGens )
# OUT: exponent vector of g1*g2
#       n^c computed via the matrix lie algebra
#
BCH_FastMultiplicationAbelianSemiTGroup := function( FCR, recBCH, g1, g2 )
    local l,hl,exp1,exp2,exp1_n,exp2_n,exp2_a,rels,i,log,logMat,mat1,mat2,
          mat, exp_n,exp,exp_a;

    # get length of abelian part and Hirsch length of N
    l := Length( FCR.factorGens );
    hl := HirschLength( FCR.recNewParent.NN );

    # get exponent vectors
    exp1 := Exponents( g1 );
    exp2 := Exponents( g2 );
    exp1_n := exp1{[l+1..l+hl]};
    exp2_n := exp2{[l+1..l+hl]};
    exp2_a := exp2{[1..l]};

    # the computation now follows the scheme:
    # g1g2 = a(g1)n(g1) a(g2)n(g2)
    #      = a(g1)a(g2)  n(g1)^a(g2) n(g2)

    # get exponent vector of a(g1*g2) and reduce mod relative orders
    exp_a := exp1{[1..l]} + exp2{[1..l]};
    rels := RelativeOrdersOfPcp( Pcp( FCR.recNewParent.GG ) );
    for i in [1..l] do
        if rels[i] <> 0 then
            exp_a[i] := exp_a[i] mod rels[i];
        fi;
    od;

    # map n(g1) to lie algebra
    log := BCH_AbstractLogCoeff_Simple_ByExponent( 
                                        FCR.recLieAlg, recBCH, exp1_n ); 

    # apply a(g2) to it; thus we get log( n(g1)^a(g2) )
    for i in [1..l] do
        log := log*FCR.lieAuts[i]^exp2_a[i];
    od;

    # compute the normal form of    exp(log)n(g_2)
    exp_n := BCH_ExpOfStarProduct( FCR, log, exp2_n, recBCH );

    exp :=  Concatenation( exp_a, exp_n );

    return exp;
end;

BCH_FastMultiplicationAbelianSemiTGroup_RuntimeVersion := function( args )
    local FCR, recBCH, g1, g2, exp;
    FCR := args[1];
    recBCH := FCR.recBCH;
    g1 := args[2];
    g2 := args[3];
    
    exp := BCH_FastMultiplicationAbelianSemiTGroup( FCR, recBCH, g1, g2 );
    return exp;
end;

#############################################################################
#
# Test functions
#


BCH_Test_ExpOfStarProduct := function( FCR,  recBCH )
    local l,hl,n1,n2,exp1,exp2,exp1_n,exp2_n,log,exp_n,n,exp_n_usual,
          log1,log2,w_log1,w_log2,pro,exp_n_bch,class;

    # get length of abelian part and Hirsch length of N
    l := Length( FCR.factorGens );
    hl := HirschLength( FCR.recNewParent.NN );
    class := FCR.recLieAlg.recTGroup.class;

    # get two random elements in N
    n1 := Random( FCR.recNewParent.NN );
    n2 := Random( FCR.recNewParent.NN );
    exp1 := Exponents( n1 );
    exp2 := Exponents( n2 );
    exp1_n := exp1{[l+1..l+hl]};
    exp2_n := exp2{[l+1..l+hl]};


    # map first one to Lie algebra, and apply ExpOfStarProduct
    log := BCH_AbstractLogCoeff_Simple_ByExponent( 
                               FCR.recLieAlg, recBCH, exp1_n );
    exp_n := BCH_ExpOfStarProduct( FCR, log, exp2_n, recBCH );
    

    # compute product in usual way
    n := n1*n2;
    exp_n_usual := Exponents(n){[l+1..l+hl]};


    # compute product with bch star
    log1 := BCH_AbstractLog_Simple_ByExponent( FCR.recLieAlg, recBCH, exp1_n );
    log2 := BCH_AbstractLog_Simple_ByExponent( FCR.recLieAlg, recBCH, exp2_n );
    w_log1 := BCH_WeightOfLieAlgElm ( FCR.recLieAlg, log1 );
    w_log2 := BCH_WeightOfLieAlgElm ( FCR.recLieAlg, log2 );
    pro  := BCH_Star_Simple( recBCH, log1, log2, w_log1, w_log2, class, 
                                                            "lieAlgebra"  );
    exp_n_bch := BCH_Abstract_Exponential_ByElm( recBCH, FCR.recLieAlg, pro );
           

    # compare
    return [ exp_n = exp_n_usual, exp_n_bch = exp_n_usual] ;

end;
