#############################################################################
##
#W setup.gi               GUARANA package                     Bjoern Assmann
##
## Methods for setting up the data structures that are needed for 
## Mal'cev collection.
##
#H  @(#)$Id$
##
#Y 2006
##
##

#############################################################################
##
#F GUARANA.PcpGroupByPcs( G, pcs )
##
## IN
## G ........................ polycyclically presented group
## pcs ...................... polycyclic sequence of G
##
## OUT
## Polycyclically presented group isomorphic to G, with 
## a polycyclic presentation with respect to pcs. 
##
GUARANA.PcpGroupByPcs := function (G, pcs )
    local sers, l, H, i;
    # built corresponding series
    sers := [];
    l := Length( pcs );
    for i in [1..l+1] do
	Add( sers, Subgroup( G, pcs{[i..l]} ));
    od;

    # get pcp 
    H := PcpGroupBySeries( sers );

    return H;
end;
    
#############################################################################
##
#F GUARANA.InitialSetupCollecRecord( args )
##
## IN
## G ................... polycyclic group whose polycyclic presentation
##                       is given with respect to a nice polycyclic 
##                       sequence, that is a sequence that fullfills
##                       all requirements that are needed for the 
##                       Malcev collection.
##         pcs = (f_1,...f_r,c_{r+1},...,c_{r+s},n_{r+s+1},...,n_{r+s+t})
##         (n_{r+s+1},...,n_{r+s+t}) is a Mal'cev basis of a normal
##                       T group N. 
##                       C =< c_{r+1},...,c_{r+s} > is a T-group and an 
##                       almost nilpotent supplement to N in G. 
##                       CN/N is free abelian.
##                       G/CN is finite.
## indeces ............. list containing three lists l1,l2,l3
##                       l1 = [1...r]
##                       l2 = [r+1,...,r+s]
##                       l3 = [r+s+1,...,r+s+t]
## N ................... normal T-group of N. 
## NN .................. pcp group isomorphic to N. 
##                       Note that in the gap sense NN is not a subgroup 
##                       of G. N and NN must have the same underlying 
##                       Mal'cev basis. 
## C ..................  T-group as described above. 
## CC .................  pcp group isomorphic to C. 
##                       The polycyclic presentation of CC must be 
##                       given with respect to Mal'cev basis 
##                       starting with (c_{r+s},...,c_{r+s}).
##                       CC!.bijection must be a bijection between 
##                       CC and C. 
##
## Example
## ll := GUARANA.SomePolyMalcevExams( 3 );
## R := GUARANA.InitialSetupCollecRecord( ll );
##
## TODO
## Maybe add some tests to check whether the input is correct.
##
GUARANA.InitialSetupCollecRecord := function( args )
    local G, indeces, N, NN, C, CC,lengths;
    G := args[1];
    indeces := args[2];
    N := args[3];
    NN := args[4];
    C := args[5];
    CC := args[6];
    lengths := [ Length(indeces[1]), Length( indeces[2] ),
	         Length( indeces[3] ) ];
    return rec( G := G, indeces := indeces, lengths := lengths,
                N := N, NN := NN, 
                C := C, CC := CC,
		collCN_method := "star" );
end;

#############################################################################
##
#F GUARANA.AddTGroupRecs( malcevRec )
#F GUARANA.AddLieAlgRecs( malcevRec )
##
## IN
## malcevRec ............... record for malcev collection
##
## EFFECT
## T-group records of C,N, 
## respectively Lie algebra records of L(N),L(C) are 
## added to malcevRec. 
##
GUARANA.AddTGroupRecs := function( malcevRec )
    malcevRec.recNN := GUARANA.TGroupRec( [malcevRec.NN ]);
    malcevRec.recCC := GUARANA.TGroupRec( [malcevRec.CC] );
end;

GUARANA.AddLieAlgRecs := function( malcevRec )
    malcevRec.recL_NN := GUARANA.LieAlgebraByTGroupRec( [malcevRec.recNN,
                                                        "gen"] );
    malcevRec.recL_CC := GUARANA.LieAlgebraByTGroupRec( [malcevRec.recCC,
                                                        "gen"] );
end;

#############################################################################
##
#F GUARANA.ComputeLieAutMatrix( malcevRec, g )
##
## IN 
## g ....................... element of parent group G.
## malcevRec ............... malcev record of G.
##
## OUT
## matrix of the Lie algebra automorphism of L(N) corresponding to  
## the group autmorphism of N, that maps n to n^g
##
GUARANA.ComputeLieAutMatrix := function( malcevRec, g )
    local autMat, pcpN, exps, n, coeff, i;
    autMat := [];
    pcpN := Pcp( malcevRec.N );
    exps := LinearActionOnPcp( [g], pcpN )[1];
    n := Length( exps );
    for i in [1..n] do
        # compute coefficients of log n_i^g
        coeff := GUARANA.AbstractLogCoeff_Simple_ByExponent( 
		 malcevRec.recL_NN,  exps[i] ); 
        Add( autMat, coeff );
    od;
    return autMat;
end;
 
#############################################################################
##
#F GUARANA.AddLieAuts( malcevRec )
##
## IN 
## malcevRec ............... malcev record of G
##
## EFFECT
## all lie automorphism matrices of the polycyclic sequence of G
## are added to malcevRec. 
##
GUARANA.AddLieAuts := function( malcevRec )
    local lieAuts, pcpG, g, autMat, i;
    lieAuts := [];
    pcpG := Pcp( malcevRec.G );

    # compute lie auts corresponding to conjugation action
    for i in [1..Length(pcpG)] do 
	g := pcpG[i];
        autMat := GUARANA.ComputeLieAutMatrix( malcevRec, g );
	Add( lieAuts, autMat );
    od;
    malcevRec.lieAuts := lieAuts;
end;

#############################################################################
##
#F GUARANA.MapFromCCtoNN( malcevRec, cc )
##
## IN
## malcevRec ...................... malcev record of G
## cc ............................. group element of CC corresponding
##
## OUT
## the corresponding elment in NN if c in C \cap N. 
## fail otherwise.
## 
GUARANA.MapFromCCtoNN := function( malcevRec, cc )
    local isoCC_C, c, isoNN_N, nn;
    
    # map from CC to C \cap N
    isoCC_C := malcevRec.CC!.bijection;
    c := Image( isoCC_C, cc );

    # map from C \cap N to NN 
    if c in malcevRec.N then
	isoNN_N := malcevRec.NN!.bijection;
	nn := PreImage( isoNN_N, c );
	return nn;
    else
	return fail;
    fi;
end;

#############################################################################
##
#F GUARANA.AddImgsOf_LCcapN_in_LN( malcevRec )
##
## IN
## malcevRec ....................... malcev record
## 
## EFFECT
## Coefficient vectors with respect to the basis of L(N) 
## of the part of the Malcev basis
## of C that lies in N  are stored.
##
GUARANA.AddImgsOf_LCcapN_in_LN := function( malcevRec )
    local hl_fa, hl_c, basis, elms_NN, nn, imgs_LCapN_in_LN, 
          expVec_nn, log_nn, i, args;

    # Get malcev basis of  C \cap N < C.
    # Note that this malcev basis corresponds to the basis of 
    # L(C \cap N ) < L(C )
    hl_fa := Length( malcevRec.indeces[2] ); # hirschlength of CN/N
    hl_c  := HirschLength( malcevRec.CC ); 
    basis := Pcp( malcevRec.CC ){[hl_fa+1..hl_c]};

    # get its image in NN
    elms_NN := [];
    for i in [1..Length(basis)] do
	nn := GUARANA.MapFromCCtoNN( malcevRec, basis[i] );
	Add( elms_NN, nn );
    od;

    # get its image in L(N)
    imgs_LCapN_in_LN := [];
    for i in [1..Length(elms_NN)] do 
	nn := elms_NN[i];
	expVec_nn := Exponents( nn );
	args := [malcevRec.recL_NN, expVec_nn, "vecByVec" ];
	log_nn := GUARANA.AbstractLog( args );
        Add( imgs_LCapN_in_LN, log_nn );
    od;

    malcevRec.imgs_LCapN_in_LN := imgs_LCapN_in_LN;
end;

#############################################################################
##
#F GUARANA.MapFromLCcapNtoLN( malcevRec, l_cc )
##
## IN
## malcevRec .................... malcev record
## l_cc ......................... lie algebra element of L(C)
##                                (in fact L(CC)) given by 
##                                coefficient vector.
##
## OUT
## Coefficient vector of the corresponding in element in L(NN)
## if l_cc in L(C\cap N ). Ohterwise fail is returned.
##
GUARANA.MapFromLCcapNtoLN := function( malcevRec, l_cc )
    local hl_fa, coeffs_1, hl_C, coeffs_2, imgs;
    
    # test wheter l_cc in L(C\cap N ) by checking whether the first 
    # hl(CN/N) coordinates of l_cc are equal to zero.
    hl_fa := Length( malcevRec.indeces[2] );# hirsch lenght of CN/N
    coeffs_1 := l_cc{[1..hl_fa]};
    if coeffs_1 <> coeffs_1*0 then
        return fail;	
    fi;

    # map it to L(N)
    hl_C := Length( l_cc );
    coeffs_2 := l_cc{[hl_fa+1..hl_C]};
    imgs := malcevRec.imgs_LCapN_in_LN;
    if Length( coeffs_2 ) <> Length( imgs ) then
	Error( " " );
    fi;
    return LinearCombination( imgs, coeffs_2 );
end;

#############################################################################
##
#F GUARANA.G_CN_WeightVector( rels )
##
## IN
## rels ....................... relative orders of G/CN
##
## OUT
## A weight vector that can be used to associate a number
## to the exponent vector of a group element of G/CN
## and vice versa.
##
GUARANA.G_CN_WeightVector := function( rels )
    local vec, n, i;

    # catch trivial case
    if Length( rels ) = 0 then
	return [];
    fi;
    
    vec := [1];
    n := Length( rels );
    for i in [1..n-1] do
	Add( vec, vec[i]*rels[n-i+1] );
    od;
    return Reversed( vec );
end;

#############################################################################
##
#F GUARANA.G_CN_ExpVectorToNumber( exp, w_vec )
#F GUARANA.G_CN_NumberToExpVector( num, w_vec )
##
## IN
## exp .................... exponent vector of an element of G/CN
## n ...................... number of an element of G/CN
## w_vec .................. a weight vector 
##
## OUT
## The number corresponding to exp, or the exponent vector corresponding to n
##
GUARANA.G_CN_ExpVectorToNumber := function( exp, w_vec )
    local num;
    num := ScalarProduct( exp, w_vec );
    return num+1;
end;
    
GUARANA.G_CN_NumberToExpVector := function( num, w_vec )
    local exp, n, div, i, num2;
   
    num2 := num -1;
    exp := [];
    n := Length( w_vec );
    for i in [1..n] do
	div := Int( num2/w_vec[i] );
	Add( exp, div );
	num2 := num2 - div*w_vec[i];
    od;
    return exp;
end;

#############################################################################
##
#F GUARANA.G_CN_InitialSetup( malcevRec )
#F GUARANA.G_CN_AddMultTable( malcevRec )
#F GUARANA.G_CN_Setup( malcevRec )
##
## IN
## malcevRec ......................... malcevRecord
## 
## EFFECT
## The multiplication table of G/CN is added to malcevRec 
##
GUARANA.G_CN_InitialSetup := function( malcevRec )
    local rels_G, rels, w_vec, order;

    # compute weight vector 
    rels_G := RelativeOrdersOfPcp( Pcp(malcevRec.G ) );
    rels := rels_G{ malcevRec.indeces[1] };
    w_vec := GUARANA.G_CN_WeightVector( rels );

    # compute order of the finite quotient G/CN 
    order := Product( rels );

    # add info record
    malcevRec.G_CN := rec( w_vec := w_vec, 
                           order := order,
		           rels := rels );
end;

GUARANA.G_CN_AddMultTable := function( malcevRec )
    local prods, pcpG, order, w_vec, exp_g, exp_h, g, h, k, i, j;

    prods := [];
    pcpG := Pcp( malcevRec.G );
    
    # catch trivial case 
    if Length( malcevRec.G_CN.rels ) = 0 then 
	return prods;
    fi;
    
    # go through all possible products and store their normal form
    order := malcevRec.G_CN.order;
    w_vec := malcevRec.G_CN.w_vec;
    for i in [1..order] do
	Add( prods, [] );
	for j in [1..order] do 
	    exp_g := GUARANA.G_CN_NumberToExpVector( i, w_vec );
	    exp_h := GUARANA.G_CN_NumberToExpVector( j, w_vec );
	    g := GUARANA.GrpElmByExpsAndPcs( pcpG, exp_g );
	    h := GUARANA.GrpElmByExpsAndPcs( pcpG, exp_h );
	    k := g*h;
	    Add( prods[i], Exponents( k ) );
	od;
    od;
    
    malcevRec.G_CN.multTable := prods;
end;

GUARANA.G_CN_Setup := function( malcevRec )
    GUARANA.G_CN_InitialSetup( malcevRec );
    GUARANA.G_CN_AddMultTable( malcevRec );
end;

#############################################################################
##
#F GUARANA.G_CN_LookUpProduct( malcevRec, exp1, exp2 )
##
## IN
## malcevRec ....................... malcev record
## exp1, exp2 ...................... exponent vector of two elments g1,g2
##                                   of G/CN
##
## OUT 
## Exponent vector ( with respect to the pcs of G ) of g1*g2
##
GUARANA.G_CN_LookUpProduct := function( malcevRec, exp1, exp2 )
    local w_vec, num1, num2;
    w_vec := malcevRec.G_CN.w_vec;
    num1 := GUARANA.G_CN_ExpVectorToNumber( exp1, w_vec );
    num2 := GUARANA.G_CN_ExpVectorToNumber( exp2, w_vec );
    return malcevRec.G_CN.multTable[num1][num2];
end;

#############################################################################
##
#F GUARANA.AddCompleteMalcevInfo( malcevRec )
##
## IN
## malcevRec ................. malcev record as given by 
##                             InitialSetupCollecRecord.
## EFFECT
## All information/data structures that are needed to do the 
## Malcev collection are added.
## 
## Example
## ll := GUARANA.SomePolyMalcevExams( 3 );
## R := GUARANA.InitialSetupCollecRecord( ll );
## GUARANA.AddCompleteMalcevInfo( R );
##
GUARANA.AddCompleteMalcevInfo := function( malcevRec )
    GUARANA.AddTGroupRecs( malcevRec );
    GUARANA.AddLieAlgRecs( malcevRec );
    GUARANA.AddLieAuts( malcevRec );
    GUARANA.AddImgsOf_LCcapN_in_LN( malcevRec );
    GUARANA.G_CN_Setup( malcevRec );
end;

#############################################################################
##
#E
