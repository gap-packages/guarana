#############################################################################
##
#W almsup.gi              GUARANA package                     Bjoern Assmann
##
## Methods for computing nilpotent almost supplements
##
#H  @(#)$Id$
##
#Y 2006

#############################################################################
#
# In: A abelian group, A < G, A given as Pcp group, A normal in G
#     U < G
#     G parent group
#
# Out: C_A(U)
#
GUARANA.Supple_CentralizerAbelianGroup := function( A, U, G )
    local C,g,pcp,rel,mat,fix,i,elms;
    # setup
    C := A;
    elms := GeneratorsOfGroup( U );

    # find iterated centralizers under the action of elms
    for g in  elms  do
        # get pcp and its relation matrix
        # This step applies only if we have finite relative orders.
        pcp := Pcp( C );
        if Length( pcp ) = 0 then return C; fi;
        rel := ExponentRelationMatrix( pcp );

        # compute action by g on pcp
        mat := LinearActionOnPcp( [g], pcp )[1];
        mat := mat - mat^0;
        Append( mat, rel );

        # compute fixed point space
        fix := PcpNullspaceIntMat( mat, Length( mat ) );
        for i in [1..Length(fix)] do
            fix[i] := MappedVector( fix[i]{[1..Length(pcp)]}, pcp );
        od;
 
        # get subgroup which centralizes all elements of elms, which 
        # were picked up so far.
        C := Subgroup( G, fix );
    od;
    return C;
end;

#############################################################################
#
# In: K < U < G, K normal in U
#
# Out: C_K(U) 
#       which is the intersection of Z(U) and K
#
GUARANA.Supple_CentralizerSubgroup := function( K, U )
    local A; 

    # centre of K is normal in U, since Z(K) is charac. in K and
    # K normal in U (thus conj. by u is an autom.)
    A := Centre( K );
    return GUARANA.Supple_CentralizerAbelianGroup( A, U, U );
end;

GUARANA.Supple_NaturalHomomorphismByPcpNC := function ( pcp )
    local  G, F, N, gens, imgs, hom;
    G := GroupOfPcp( pcp );
    N := SubgroupByIgs( G, DenominatorOfPcp( pcp ) );
    F := PcpGroupByPcp( pcp );
    UseFactorRelation( G, N, F );
    gens := ShallowCopy( GeneratorsOfPcp( pcp ) );
    imgs := ShallowCopy( Igs( F ) );
    Append( gens, DenominatorOfPcp( pcp ) );
    Append( imgs, List( DenominatorOfPcp( pcp ), function ( x )
            return One( F );
        end ) );
    hom := GroupHomomorphismByImagesNC( G, F, gens, imgs );
    SetKernelOfMultiplicativeGeneralMapping( hom, N );
    return hom;
end;


GUARANA.Supple_NaturalHomomorphismNC := function ( G, N )
    if Size( N ) = 1  then
        return IdentityMapping( G );
    fi;
    return GUARANA.Supple_NaturalHomomorphismByPcpNC( Pcp( G, N ) );
end;

#############################################################################
#
# In: K < U, K normal in U 
#
# Out: upper U-central series of K, which may terminate a proper subgroup
#      of K
# 
GUARANA.Supple_UpperU_CentralSeries := function( K, U )
    local C, upp, N, nat, U_img, K_img, C_img;

    C := TrivialSubgroup( U );
    upp := [C];
    N := GUARANA.Supple_CentralizerSubgroup( K, U );
    while IndexNC( N, C ) > 1 do
        C := N;
        Add( upp, C );
        nat := GUARANA.Supple_NaturalHomomorphismNC( U, C );
        U_img := Image( nat  );
        K_img := Image( nat, K );
        C_img := GUARANA.Supple_CentralizerSubgroup( K_img, U_img );
        N := PreImage( nat, C_img );
    od;
    return Reversed( upp );
   
end;

#############################################################################
# the following code would be suitable for a more general approach.
# It has to be revised

# A is normal subgroup of the parent group G
# gens .... elements of parent group G
# pcp  .... pcp of A/L
# nat  .... hom from A to A/L
GUARANA.Supple_LinearActionOnFactorPcp := function ( gens, pcp, nat )
    local result, x, act, y,y1,y2,y3;

    result := [];
    for x in gens do    
        act := [];
        for y in AsList( pcp ) do 
            y1 := PreImagesRepresentative( nat, y );
            y2 := y1^x;
            y3 := Image( nat, y2 );
            Add( act, ExponentsByPcp( pcp, y3 ));
        od;   
        Add( result, act );
    od;
    return result;
end;

GUARANA.Supple_UHyperCentre := function( U, A )
    local pcpU, gensU,L,nat,AL,g,pcp,rel,mat,fix,i,C,m,l,Matt,Mat;

    pcpU  := Pcp( U );
    gensU := AsList( pcpU );
    L     := TrivialSubgroup( A );

    while true do
        nat := NaturalHomomorphism( A, L );
        AL := Image( nat );

        # compute C_A/L( U ), i.e. the elements of A/L which are fixed under
        # the conjugation action of U

          # get pcp and its relation matrix, only necessary if A has factors
          # with finite orders
          pcp := Pcp( AL );
          rel := ExponentRelationMatrix( pcp );
          if Length( pcp ) = 0 then return A; fi;

          # compute action by gensU on pcp
          Mat := [];
          for g in gensU do
              mat := GUARANA.Supple_LinearActionOnFactorPcp( [g], pcp, nat )[1];
              mat := mat - mat^0;
              Append( mat, rel );
              Add( Mat, mat );
          od;

          # write everything in one big matrix (g_1 g_2 ... )
          m := Length( mat );
          Matt := [];
          for i in [1..m] do
              l := List( Mat, x-> x[i] );
              Add( Matt, Concatenation( l ) );
          od;
          
          # compute fixed point space of U in AL
          fix := PcpNullspaceIntMat( Matt, Length( Matt ) );
          for i in [1..Length(fix)] do
              fix[i] := MappedVector( fix[i]{[1..Length(pcp)]}, pcp );
          od;
          C := Subgroup( AL, fix );
          
        # test if L is already big enough
        if Size( C ) = 1 then
            return L;
        fi; 
        L := PreImage( nat, C );
    od;          
end;

#############################################################################
##
#F GUARANA.Supple_RefinedLowerCentralSeries( <G> )
##
## IN
## G ... nilpotent polycyclic group
##
## OUT
## efa series of G which refines the lower central series by
## finite - by - (torsion-free) factors.
##
GUARANA.Supple_RefinedLowerCentralSeries := function( G )
    local ser, ref, i, A, B, pcp, gens, rels, n, free, fini, U, s, t, f;

    ser := LowerCentralSeries( G );
    ref := [G];
    for i in [1..Length( ser ) - 1] do
        # refine abelian factor A/B
        A := ser[i];
        B := ser[i+1];
        pcp := Pcp( A, B, "snf" );
        gens := GeneratorsOfPcp( pcp );
        rels := RelativeOrdersOfPcp( pcp );
        n    := Length( gens );

        # take the free part for the top factor
        free := Filtered( [1..n], x -> rels[x] = 0 );
        fini := Filtered( [1..n], x -> rels[x] > 0 );
        if Length( free ) > 0 then
            f := AddToIgs( Igs(B), gens{fini} );
            U := SubgroupByIgs( G, f );
            Add( ref, U );
        else
            U := A;
        fi;

        # the torsion subgroup
        if Length( fini ) > 0 then
            s := Factors( Lcm( rels{fini} ) );
            f := gens{fini};
            for t in s do
                f := List( f, x -> x ^ t );
                f := AddToIgs( Igs(B), f );
                U := SubgroupByIgs( G, f );
                Add( ref, U );
            od;
        fi;
    od;
    return ref;
end;

#############################################################################
##
#F GUARANA.Supple_NilpAlmSup( G, N, H )
##
## IN
## G ... polycyclic group
## N ... normal nilpotent subgroup of G
## H ... subgroup with, H/N abelian and [G:H]< \infty
##
## OUT
## nilpotent almost supplement for N in G
## 
## TODO: Ueberlege ob ich immer mit den Homomorphismen rechnen muss.
## z.B. beim U-hyper centre
##
GUARANA.Supple_NilpAlmSup := function( G, N, H )
    local sers,U,l,i,nat,A,sers2,L,nat_L,N_i_img,U_img,CR,s,K_img;

    # determine H-invariant linear abelian series of N with central factors 
    sers := GUARANA.Supple_RefinedLowerCentralSeries( N );

    U := H;
    l := Length( sers ) - 1;

    for i in [1..l] do
        # hom U -> U/N_i+1
        nat := GUARANA.Supple_NaturalHomomorphismNC( U, sers[i+1] );
        
        # compute ascending U central series of A = N_i/N_i+1
        A := Image( nat, sers[i] );
        U_img := Image( nat );
        sers2 := GUARANA.Supple_UpperU_CentralSeries( A, U_img );
          
        # proceed only if series does not reach A, since
        # if A = U-hpercentre, U/N_i+1 is nilpotent
        if sers2[1] <> A  then
            # get  U-hypercentre L/N_i+1 in N_i/N_i+1
            L := PreImage( nat, sers2[1]);

            # hom U -> U/L
            nat_L := GUARANA.Supple_NaturalHomomorphismNC( U, L );
            # N_i/L and U/L
            N_i_img := Image( nat_L, sers[i] );
            U_img := Image( nat_L ); 
            if Order( A ) < infinity then
                # compute complement K/L to N_i/L in U/L
                CR := CRRecordBySubgroup( U_img, N_i_img );
                s :=  OneCocyclesEX( CR ).transl;
                K_img := ComplementCR( CR, s );
            else 
                # compute almost complement K/L to N_i/L in U/L
                K_img := GUARANA.Supple_AlmostComplement( U_img , N_i_img );
            fi;
            U := PreImage( nat_L, K_img );
        fi;
    od;
    return U;
end;

GUARANA.Supple_ComputeHAndN := function( G )
    local N, nat, G_img, G_img_fit, H_img;
    N := FittingSubgroup( G );
    nat := GUARANA.Supple_NaturalHomomorphismNC( G, N );
    G_img := Image( nat );
    G_img_fit := FittingSubgroup( G_img );
    H_img := Centre( G_img_fit );
    return rec( G := G, H := PreImage( nat, H_img ), N := N );
end;
