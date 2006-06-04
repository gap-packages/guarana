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
## R := GUARANA.InitialSetupMalcevCollector( ll );
##
## TODO
## Maybe add some tests to check whether the input is correct.
##
GUARANA.InitialSetupMalcevCollector := function( args )
    local G, indeces, N, NN, C, CC, lengths, cn_fam, cn_elms_type, obj;

    G := args[1];
    indeces := args[2];
    N := args[3];
    NN := args[4];
    C := args[5];
    CC := args[6];
    lengths := [ Length(indeces[1]), Length( indeces[2] ),
	         Length( indeces[3] ) ];

    # create family and types for elements of CN and ...
    cn_fam := NewFamily( "MalcevCNFamily",
                          IsMalcevCNElement,
                          IsMalcevCNElement );
    cn_elms_type := NewType( cn_fam, IsMalcevCNElementRep );

    obj := rec( G := G, indeces := indeces, lengths := lengths,
                N := N, NN := NN, 
                C := C, CC := CC,
                mo_NN := "unknown" , mo_CC := "unknown",
                cn_fam := cn_fam,
                cn_elms_type := cn_elms_type );

    obj :=  Objectify( NewType( MalcevCollectorFamily, 
                               IsMalcevCollectorRep and
			       IsMutable ),
		      obj );
    return obj;
end;

#############################################################################
##
## Funtions to view and print a Malcev collector
##
InstallMethod( PrintObj, "for Malcev collector", [ IsMalcevCollectorRep ],
function( malColl )
    Print( "<<Malcev collector>>\n",
           "  C : ", malColl!.mo_CC, "\n",
           "  N : ", malColl!.mo_NN  );
end );

#F GUARANA.MO_AddTGroupRecs( malcevCollector )
#F GUARANA.AddMalcevObjects( malcevCollector )
##
##
## EFFECT
## T-group records of C,N, 
## respectively Malcev objects of L(N),L(C) are 
## added to malCol. 
##
GUARANA.MO_AddTGroupRecs := function( malcevCollector )
    malcevCollector!.recNN := GUARANA.TGroupRec( [malcevCollector!.NN ]);
    malcevCollector!.recCC := GUARANA.TGroupRec( [malcevCollector!.CC] );
end;

GUARANA.AddMalcevObjects := function( malCol )
    malCol!.mo_NN := MalcevObjectConstruction( malCol!.recNN ); 
    malCol!.mo_CC := MalcevObjectConstruction( malCol!.recCC ); 
end;

#############################################################################
##
#F GUARANA.MO_ComputeLieAutMatrix( malCol, g )
##
## IN 
## g ....................... element of parent group G.
## malCol .................. Malcev collector
##
## OUT
## matrix of the Lie algebra automorphism of L(N) corresponding to  
## the group autmorphism of N, that maps n to n^g
##
GUARANA.MO_ComputeLieAutMatrix := function( malCol, g )
    local pcpN, autMat, exps, n, n_i_g, log_n_i_g, coeff, i;

    # catch trivial case
    if HirschLength( malCol!.N ) = 0 then 
        return EmptyMatrix( 0 );
    fi;

    pcpN := Pcp( malCol!.N );
    autMat := [];
    exps := LinearActionOnPcp( [g], pcpN )[1];
    n := Length( exps );
    for i in [1..n] do
        # compute coefficients of log n_i^g
        n_i_g := MalcevGrpElementByExponents( malCol!.mo_NN, exps[i] );
        log_n_i_g := Log( n_i_g );
        coeff := Coefficients( log_n_i_g ); 
        Add( autMat, coeff );
    od;
    return autMat;
end;
 
#############################################################################
##
#F GUARANA.Add_C_LieAuts( malCol )
##
## all lie automorphism matrices corresponding to the basis elms of C
## are added to the Malcev collector. 
##
GUARANA.Add_C_LieAuts := function( malCol )
    local lieAuts, pcpC, c, autMat, i;

    lieAuts := [];
    pcpC := Pcp( malCol!.C );

    # compute lie auts corresponding to conjugation action
    for i in [1..Length(pcpC)] do 
	    c := pcpC[i];
        autMat := GUARANA.MO_ComputeLieAutMatrix( malCol, c );
	Add( lieAuts, autMat );
    od;

    malCol!.C_lieAuts := lieAuts;
end;

#############################################################################
##
#F GUARANA.MO_MapFromCCtoNN( malCol, cc )
##
## IN
## malCol ...................... .. Malcev collector of G
## cc ............................. group element of CC corresponding
##
## OUT
## the corresponding elment in NN if c in C \cap N. 
## fail otherwise.
## 
GUARANA.MO_MapFromCCtoNN := function( malCol, cc )
    local isoCC_C, c, isoNN_N, nn;
    
    # map from CC to C \cap N
    isoCC_C := malCol!.CC!.bijection;
    c := Image( isoCC_C, cc );

    # map from C \cap N to NN 
    if c in malCol!.N then
	    isoNN_N := malCol!.NN!.bijection;
	    nn := PreImage( isoNN_N, c );
	    return nn;
    else
	    return fail;
    fi;
end;

#############################################################################
##
#F GUARANA.MO_AddImgsOf_LCcapN_in_LN( malCol )
##
## Coefficient vectors with respect to the basis of L(N) 
## of the part of the Malcev basis
## of C that lies in N  are stored.
##
## Example
## ll := GUARANA.SomePolyMalcevExams( 9 );
## R := GUARANA.InitialSetupMalcevCollector( ll );
## GUARANA.AddTGroupRecs( R );
## GUARANA.AddMalcevObjects( R );
##
GUARANA.MO_AddImgsOf_LCcapN_in_LN := function( malCol )
    local hl_fa, hl_c, basis, elms_NN, nn, imgs_LCapN_in_LN, nn2, log_nn, i;

    # Get malcev basis of  C \cap N < C.
    # Note that this malcev basis corresponds to the basis of 
    # L(C \cap N ) < L(C )
    hl_fa := Length( malCol!.indeces[2] ); # hirschlength of CN/N
    hl_c  := HirschLength( malCol!.CC ); 
    basis := Pcp( malCol!.CC ){[hl_fa+1..hl_c]};

    # get its image in NN
    elms_NN := [];
    for i in [1..Length(basis)] do
	    nn := GUARANA.MO_MapFromCCtoNN( malCol, basis[i] );
	    Add( elms_NN, nn );
    od;

    # get its image in L(N)
    imgs_LCapN_in_LN := [];
    for i in [1..Length(elms_NN)] do 
	    nn := elms_NN[i];
        nn2 := MalcevGrpElementByExponents( malCol!.mo_NN, Exponents(nn) );
        log_nn := Log( nn2 ); 
        Add( imgs_LCapN_in_LN, log_nn );
    od;

    malCol!.imgs_LCapN_in_LN := imgs_LCapN_in_LN;
end;

#############################################################################
##
#F GUARANA.MO_MapFromLCcapNtoLN( malCol, l_cc )
##
## IN
## malCol... .................... malcev collector
## l_cc ......................... lie algebra element of L(C)
##                                (in fact L(CC)) given 
##                                as Malcev lie element.
##
## OUT
## Coefficient vector of the corresponding in element in L(NN)
## if l_cc in L(C\cap N ). Ohterwise fail is returned.
##
##
GUARANA.MapFromLCcapNtoLN := function( malCol, l_cc )
    local hl_fa, coeffs, coeffs_1, hl_C, coeffs_2, imgs;

    # catch trivial case
    if l_cc = 0*l_cc then 
        return MalcevLieElementByWord( malCol!.mo_NN, [[],[]] ); 
    fi;
    
    # test wheter l_cc in L(C\cap N ) by checking whether the first 
    # hl(CN/N) coordinates of l_cc are equal to zero.
    hl_fa := Length( malCol!.indeces[2] );# hirsch lenght of CN/N
    coeffs := Coefficients( l_cc );
    coeffs_1 := coeffs{[1..hl_fa]};
    if coeffs_1 <> coeffs_1*0 then
        return fail;	
    fi;

    # map it to L(N)
    hl_C := Length( coeffs );
    coeffs_2 := coeffs{[hl_fa+1..hl_C]};
    imgs := malCol!.imgs_LCapN_in_LN;
    if Length( coeffs_2 ) <> Length( imgs ) then
	    Error( " " );
    fi;
    return LinearCombination( imgs, coeffs_2 );
end;

#############################################################################
##
#F GUARANA.MO_AddCompleteMalcevInfo( malCol )
##
## IN
## malCol .................... Malcev collector as given by 
##                             InitialSetupMalcevCollector.
## EFFECT
## All information/data structures that are needed to do the 
## Malcev collection are added.
## 
## Example
## ll := GUARANA.SomePolyMalcevExams( 3 );
## R := GUARANA.InitialSetupCollecRecord( ll );
## GUARANA.MO_AddCompleteMalcevInfo( R );
##
GUARANA.MO_AddCompleteMalcevInfo := function( malCol )

    GUARANA.MO_AddTGroupRecs( malCol );
    GUARANA.AddMalcevObjects( malCol );
    GUARANA.Add_C_LieAuts( malCol );
    GUARANA.MO_AddImgsOf_LCcapN_in_LN( malCol );

end;

#############################################################################
##
#F GUARANA.SetupMalcevCollector( args )
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
## OUT
## Malcev record that contains all data needed  for Mal'cev collection.
##
## Example
## ll := GUARANA.SomePolyMalcevExams( 3 );
## R := GUARANA.SetupMalcevCollector( ll );
##
##
GUARANA.SetupMalcevCollector := function( args )
    local malCol;
    malCol := GUARANA.InitialSetupMalcevCollector( args );
    GUARANA.MO_AddCompleteMalcevInfo( malCol );
    return malCol;
end;

InstallGlobalFunction( MalcevCollectorConstruction,
function( args )
    return GUARANA.SetupMalcevCollector( args );
end);

# TODO
#
# functions for creating an CN - element by an cn vector.
# 
# Creat object GUA_CN_Elm
#                  malcev_elm_c
#                  malcev_elm_n
#              Methods
#              GUA_Exponents_CN
#              GUA_Exponents_C
#              GUA_Exponents_N
#              GUA_Coefficients_LN
#              GUA_Coefficients_LC
#
#        object MalcevElement
#                   grp_elm
#                   lie_elm
#                   (not both of them have to be known).
#        need a constructor
#        a type,...
#
#              * for IsMalcevElm IsMalcevElm
#              (mulitplication can happen in group or in Lie algebra).
#              * for IsMalcevElm IsMatrix   
# Create object GUA_G_Elm
#                  exp_f
#                  GUA_CN_Elm


#############################################################################
##
#E
