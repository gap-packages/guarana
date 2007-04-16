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
    local G, indeces, N, NN, C, CC, lengths, cn_fam, cn_elms_type, obj,
          g_fam, g_elms_type, rels;

    G := args[1];
    indeces := args[2];
    N := args[3];
    NN := args[4];
    C := args[5];
    CC := args[6];
    lengths := [ Length(indeces[1]), Length( indeces[2] ),
	         Length( indeces[3] ) ];
    rels := RelativeOrdersOfPcp( Pcp(G ) );

    # create family and types for elements of CN and ...
    cn_fam := NewFamily( "MalcevCNFamily",
                          IsMalcevCNElement,
                          IsMalcevCNElement );
    cn_elms_type := NewType( cn_fam, IsMalcevCNElementRep );
    g_fam := NewFamily( "MalcevGFamily",
                          IsMalcevGElement,
                          IsMalcevGElement );
    g_elms_type := NewType( g_fam, IsMalcevGElementRep );

    obj := rec( G := G, indeces := indeces, lengths := lengths,
                N := N, NN := NN, 
                C := C, CC := CC,
                rels := rels,
                mo_NN := "unknown" , mo_CC := "unknown",
                cn_fam := cn_fam,
                cn_elms_type := cn_elms_type, 
                g_fam := g_fam, 
                g_elms_type := g_elms_type,
                mult_method := "standard" 
                # if we set this to symbolic then we have symbolic 
                # du Sautoy collection, (only prototype implementation).
                );

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
           "  F : ", malColl!.rels{malColl!.indeces[1]} , "\n", 
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
#F GUARANA.Add_G_LieAuts( malCol )
##
## all lie automorphism matrices corresponding to the pcs elms of G
## are added to the Malcev collector. 
##
GUARANA.Add_G_LieAuts := function( malCol )
    local lieAuts, pcpG, g, autMat, i;

    lieAuts := [];
    pcpG := Pcp( malCol!.G );

    # compute lie auts corresponding to conjugation action
    for i in [1..Length(pcpG)] do 
	    g := pcpG[i];
        autMat := GUARANA.MO_ComputeLieAutMatrix( malCol, g );
	Add( lieAuts, autMat );
    od;

    malCol!.G_lieAuts := lieAuts;
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

## IN elm ..................... MalcevGenElement of C
##                              that lies in the intersection of C and N
##
## OUT 
## The
GUARANA.MapFrom_MOC_2_MON := function( malCol, elm )
    local x, x_img;

    x := LieElement( elm );
    x_img := GUARANA.MapFromLCcapNtoLN( malCol, x );
    if x_img = fail then
        return fail;
    else
        return MalcevGenElementByLieElement( x_img );
    fi;
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
#F GUARANA.MO_G_CN_InitialSetup( malCol )
#F GUARANA.MO_G_CN_AddMultTable( malCol )
#F GUARANA.MO_G_CN_Setup( malCol )
##
## EFFECT
## The multiplication table of G/CN is added to the malcev collector
##
GUARANA.MO_G_CN_InitialSetup := function( malCol )
    local rels_G, rels, w_vec, order;

    # compute weight vector 
    rels_G := RelativeOrdersOfPcp( Pcp(malCol!.G ) );
    rels := rels_G{ malCol!.indeces[1] };
    w_vec := GUARANA.G_CN_WeightVector( rels );

    # compute order of the finite quotient G/CN 
    order := Product( rels );

    # add info record
    malCol!.G_CN := rec( w_vec := w_vec, 
                         order := order,
		                 rels := rels );
end;

GUARANA.MO_G_CN_AddMultTable := function( malCol )
    local prods, collG, order, w_vec, exp_g, exp_h, g, h, k, i, j, exps_k;

    prods := [];
    collG := Collector( malCol!.G );
    
    # catch trivial case 
    if Length( malCol!.G_CN.rels ) = 0 then 
	return prods;
    fi;
    
    # go through all possible products and store their normal form
    order := malCol!.G_CN.order;
    w_vec := malCol!.G_CN.w_vec;
    for i in [1..order] do
	    Add( prods, [] );
	    for j in [1..order] do 
	        exp_g := GUARANA.G_CN_NumberToExpVector( i, w_vec );
	        exp_h := GUARANA.G_CN_NumberToExpVector( j, w_vec );
	        g := GUARANA.GrpElmByExpsAndCollNC( collG, exp_g );
	        h := GUARANA.GrpElmByExpsAndCollNC( collG, exp_h );
	        k := g*h;
            # only save exp vector, since the generation of 
            # malcevGelements is a little bit expensive
	        Add( prods[i], Exponents( k ) );

                #exps_k := Exponents( k );
                #Add( prods[i], MalcevGElementByExponentsNC( malCol, exps_k ));
	    od;
    od;
    
    malCol!.G_CN.multTable := prods;
end;

GUARANA.MO_G_CN_Setup := function( malCol )
    GUARANA.MO_G_CN_InitialSetup( malCol );
    GUARANA.MO_G_CN_AddMultTable( malCol );
end;

#############################################################################
##
#F GUARANA.MO_G_CN_LookUpProduct( malCol, exp1, exp2 )
##
## IN
## malCol ....................... malcev collector
## exp1, exp2 ...................... exponent vector of two elments g1,g2
##                                   of G/CN
##
## OUT 
## MalcevGElement  g1*g2
##
GUARANA.MO_G_CN_LookUpProduct := function( malCol, exp1, exp2 )
    local w_vec, num1, num2,exp_vec;
    w_vec := malCol!.G_CN.w_vec;
    # catch case G/CN = 1
    if Length( w_vec ) = 0 then 
        return GUARANA.G_Identity( malCol );
    fi;

    # get exp vector
    num1 := GUARANA.G_CN_ExpVectorToNumber( exp1, w_vec );
    num2 := GUARANA.G_CN_ExpVectorToNumber( exp2, w_vec );
    exp_vec := malCol!.G_CN.multTable[num1][num2];

    return  MalcevGElementByExponentsNC( malCol, exp_vec );
end;

GUARANA.GetCNElmentByGExpVector := function( malCol, exps )
    local indeces, exps_cn, cn_elm;
    indeces := malCol!.indeces;
    if exps{indeces[1]} <> 0*exps{indeces[1]} then
        Error( "Exponent vector does not correspond to CN element \n" );
    else
        exps_cn := exps{ Concatenation(indeces[2],indeces[3]) };
        cn_elm := MalcevCNElementByExponents( malCol, exps_cn );
        return cn_elm;
    fi;
end;

#############################################################################
##
#F GUARANA.MO_CconjF_AddInfo( malCol )
##
## EFFECT
## Let (f_1,...,f_r,c_1,...,c_s,n_1,...,n_t ) be a pcs suitable
## for Malcev collection. This function stores
## all normal forms of the type 
## c_i^f ( where f is product of the f_j )
## in the Malcev collector.
##
GUARANA.MO_CconjF_AddInfo := function( malCol )
    local CconjF, pcpG, indeces, order, w_vec, ind_C, cconjF, 
          exp_f, c_i, f, k, exps_k, cn_elm, i, j; 
    CconjF := [];
    pcpG := Pcp( malCol!.G );
    indeces := malCol!.indeces;

    # catch trivial case 
    if Length( malCol!.G_CN.rels ) = 0 then 
        malCol!.CconjF := CconjF;
	    return 0;
    fi;
    
    order := malCol!.G_CN.order;
    w_vec := malCol!.G_CN.w_vec;
    ind_C := malCol!.indeces[2];
   
    # go through all c_i
    for i in ind_C do 
	    cconjF := [];
	    exp_f := GUARANA.G_CN_NumberToExpVector( i, w_vec );
	    # go through all elms of finite part
	    for j in [1..order] do 
	        exp_f := GUARANA.G_CN_NumberToExpVector( j, w_vec );
	        c_i := pcpG[i]; 
	        f := GUARANA.GrpElmByExpsAndPcs( pcpG, exp_f );
	        k := c_i^f;
	        # Add( cconjF, Exponents( k ) );
            exps_k := Exponents( k );
            # get corresponding Malcev CN element
            cn_elm := GUARANA.GetCNElmentByGExpVector( malCol, exps_k );
            # get also lie information about cn
            Coefficients( cn_elm!.c );
            Coefficients( cn_elm!.n );

            Add( cconjF, cn_elm );
	    od;
	    Add( CconjF, cconjF );
    od;

    malCol!.CconjF := CconjF;
end;

#############################################################################
##
#F GUARANA.MO_CconjF_LookUp( malCol, i, exp_f )
##
## IN
## malCol .......................    Malcev collector
## i ............................... index of c_i, which is a part 
##                                   of the pcs of CN/N.
##                                   i determines the position of c_i in 
##                                   the pcs of CN/N.
## exp_f ........................... short exponent vector (with respect 
##                                   to G/CN) of an element of the finite 
##                                   part.
## 
## OUT
## (c_i)^f as stored in the Malcev record.
##
GUARANA.MO_CconjF_LookUp := function( malCol, i, exp_f )
    local w_vec, num_f;

    # check input 
    if not i+malCol!.lengths[1] in malCol!.indeces[2] then
        Error( " \n" );
    fi;
    # catch trivial case G/CN =1 
    #TODO

    # get number that corresponds to f
    w_vec := malCol!.G_CN.w_vec;
    num_f := GUARANA.G_CN_ExpVectorToNumber( exp_f, w_vec ); 

    return malCol!.CconjF[i][num_f]; 
end;

#############################################################################
##
#F GUARANA.MO_F_AddInverseInfo( malCol )
##
## IN
## malCol ....................... Malcev record 
##
## EFFECT
## Store f^-1 for all representatives f of the finite part G/CN
##
GUARANA.MO_F_AddInverseInfo := function( malCol )
    local inversesF, pcpG, order, w_vec, exp_f, f, f_inv, exps, elm, j;

    inversesF := [];
    pcpG := Pcp( malCol!.G );

    # catch trivial case 
    if Length( malCol!.G_CN.rels ) = 0 then 
        malCol!.inversesF := inversesF;
	    return 0;
    fi;
    
    order := malCol!.G_CN.order;
    w_vec := malCol!.G_CN.w_vec;

    # go through all elms of finite part
    for j in [1..order] do 
	    exp_f := GUARANA.G_CN_NumberToExpVector( j, w_vec );
	    f := GUARANA.GrpElmByExpsAndPcs( pcpG, exp_f );
	    f_inv := f^-1;   
	    #Add( inversesF, Exponents( f_inv ) );
        exps := Exponents( f_inv );
        elm := MalcevGElementByExponents( malCol, exps );
        Add( inversesF, elm );
    od;
    malCol!.inversesF := inversesF;
    return 0;
end;

#############################################################################
##
#F GUARANA.MO_F_LookupInverse( malCol, exp_f )
##
## IN
## malCol ........................... malcev Record
## exp_f ........................... short exponent vector (with respect 
##                                   to G/CN) of an element of the finite 
##                                   part.
## 
## OUT
## f^-1  as stored in the Malcev record.
##
GUARANA.MO_F_LookupInverse := function( malCol, exp_f )
    local w_vec, num_f;

    # catch trivial case G/CN = 1
    if Length( exp_f )=0 then 
        return GUARANA.G_Identity( malCol );
    fi;

    # get number that corresponds to f
    w_vec := malCol!.G_CN.w_vec;
    num_f := GUARANA.G_CN_ExpVectorToNumber( exp_f, w_vec ); 

    return malCol!.inversesF[num_f];
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
    GUARANA.Add_G_LieAuts( malCol );
    GUARANA.MO_AddImgsOf_LCcapN_in_LN( malCol );
    GUARANA.MO_G_CN_Setup( malCol );
    GUARANA.MO_CconjF_AddInfo( malCol );
    GUARANA.MO_F_AddInverseInfo( malCol );
end;

GUARANA.SetAllMethodsToSymbolic := function( malCol )
    local mo_CC, mo_NN;

    mo_CC := malCol!.mo_CC;
    mo_NN := malCol!.mo_NN;
    
    SetLogAndExpMethod( mo_CC, "pols" ); 
    #SetStarMethod( mo_CC, "pols" );

    SetLogAndExpMethod( mo_NN, "pols" ); 
    #SetStarMethod( mo_NN, "pols" );
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
    GUARANA.SetAllMethodsToSymbolic( malCol );

    # add dt polynomials and set them as defaul method
    # TODO: Should we add them before we setup the malcev correspondence ?
    AddDTPolynomials( malCol!.mo_NN );
    AddDTPolynomials( malCol!.mo_CC );
    SetMultiplicationMethod( malCol!.mo_NN , GUARANA.MultMethodIsCollection );
    SetMultiplicationMethod( malCol!.mo_CC, GUARANA.MultMethodIsCollection );

    return malCol;
end;

InstallGlobalFunction( MalcevCollectorConstruction,
function( args )
    return GUARANA.SetupMalcevCollector( args );
end);

#############################################################################
##
#E
