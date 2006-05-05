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
GUARANA.InitialSetupCollecRecord := function( args )
    local G, indeces, N, NN, C, CC;
    G := args[1];
    indeces := args[2];
    N := args[3];
    NN := args[4];
    C := args[5];
    CC := args[6];
    return rec( G := G, indeces := indeces, N := N, NN := NN, 
                C := C, CC := CC );
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

# TODO
# setup functions
# -function for mapping elms from C\capN < C to 
#  and more importantly for the elms of the Lie algebra. 
# -multiplication table of finite part on top. 
# -function for the full setup. 
# 

## Example
## ll := GUARANA.SomePolyMalcevExams( 3 );
## R := GUARANA.InitialSetupCollecRecord( ll );
## GUARANA.AddCompleteMalcevInfo( R );
##
GUARANA.AddCompleteMalcevInfo := function( malcevRec )
    GUARANA.AddTGroupRecs( malcevRec );
    GUARANA.AddLieAlgRecs( malcevRec );
    GUARANA.AddLieAuts( malcevRec );
end;
# 
# collection functions
# -idea represent internally elments only as vectors
#  (maybe splitted in three parts ? )
# -computations with powers of autom of N. 
#  Two possibilities, n given by group exponent vector or by 
#  Lie algebra coeff vector
# -computations with consecutive powers of automorphisms. 
#  again n given in two ways. 
# -Powering in C. 
# -Inversion in CN
# -Powering in CN
# -Collection in G
#
# construct more examples
# 

#############################################################################
##
#E
