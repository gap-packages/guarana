# The code in this file was used to construct the examples used 
# in the conference paper "Algorithmic use of the Mal'cev correspondence"
#
# 
#

LoadPackage( "nq" );


#############################################################################
#
# Code for the computation of polycyclic presentations of triangular
# matrix groups over a maximal order of a number field
#
#############################################################################

FindNonZeroEntry := function ( mat )
    # only for unipotent matrices with one non-zero entry
    local n,i,j;
    n := Length( mat );
    for i in [1..n-1] do
        for j in [i+1..n] do 
            if mat[i][j] <> 0 then
                return [[i,j],mat[i][j]];
            fi;
        od;
    od;
    return fail;
end;

DetermineRelationsD := function( col, U_O, dim )
    local o,r,i,gens;

    gens := GeneratorsOfGroup( U_O );
    # order of fundamental unit
    o := Order( gens[1] );
    r := Length( gens );

    # set power relation of torsion unit
    for i in [0..dim-1] do
        SetRelativeOrder( col, 1+i*r, o );
    od;
    return;
end;

PositionsUnipotentGens := function( pos, n, dim, noD  )
    # pos = [i,j] specify a matrix entry
    # n is the dimension of the maximal order
    # dim is the dimension of the matrices
    local i,j,weight,sum,ll,k;
    i := pos[1];
    j := pos[2];

    # number of subdiagonal
    weight := j-i;

    # number of "n-boxes" before the wanted subdiagonal
    sum := Sum( [(dim-1 -weight +2) ..(dim-1)] );
  
    ll := [];
    for k in [1..n] do 
        Add( ll, sum*n + (i-1)*n + k );
    od;
    return ll + noD;
end;

UnipoMat2Word := function( mat, O, noD, noU )
    # only for unipotent matrices with one non-zero entry
    local n,dim, ent, coeff,pos,ll,word,i;

    dim := Length( mat );
    # dimension of maximal order
    n := Length( O );

    ent :=  FindNonZeroEntry( mat );
    if ent = fail then return fail; fi;
    coeff := Coefficients( O, ent[2] );     

    pos := ent[1];
    # get the positions of the corresponding generators in the U(n,O)
    ll := PositionsUnipotentGens( pos, n, dim, noD  );
    
    # construct the word
    word := [];
    for i in [1..n] do  
        if coeff[i] <> 0 then
            Add( word, ll[i] );
            Add( word, coeff[i] );
        fi;
    od;

    return word;
end; 

ConstructD := function( dim, U_O , F)
    local ll, gens, r, i,j, mat;

    ll := [];
    gens := GeneratorsOfGroup( U_O );
    r := Length( gens ); 
 
    for i in [1..dim] do
        for j in [1..r] do
            mat :=  IdentityMat( dim, F );
            mat[i][i] := gens[j];
            Add( ll, mat );
        od;
    od;

    return ll;

end;

ConstructU := function( dim, O, F )
    local ll,w,i,mat,o;
    ll := [];
    for w in [1..(dim-1)] do
        for i in [1..(dim-w)] do
            for o in O do
                mat :=  IdentityMat( dim, F );
                mat[i][i+w] := o;
                Add( ll, mat );
            od;
        od;
    od;
    return ll;
end;

DetermineRelationsDonU := function( gensTr, noD,noU,col,O )
    local i,d,j,u,mat,word;
    for i in [1..noD] do
        d := gensTr[i];
        for j in [noD+1..noD+noU] do
            u := gensTr[j];
            mat := u^d;
            word := UnipoMat2Word( mat, O, noD, noU );
            #Print( "i ", i, " j ", j , " word ", word, "\n");
            SetConjugate( col, j, i, word );

            # also inverse
            mat := u^(d^-1);
            word := UnipoMat2Word( mat, O, noD, noU );
            #Print( "i ", i, " j ", j , " word ", word, "\n");
            SetConjugate( col, j, i, word );
        od;
    od;
end;

DetermineRelationsU := function( gensTr, noD,noU,col,O )
    local i,u1,j,u2,mat,word, word2;
    for i in [noD+1..(noD+noU)] do
        u1 := gensTr[i];
        for j in [i+1..noD+noU] do
            u2 := gensTr[j];
            mat := Comm( u2, u1 );
            #Display( mat );
            word := UnipoMat2Word( mat, O, noD, noU );
            #Print( "i ", i, " j ", j , " word ", word, "\n");
            if not word = fail then 
                word2 := [j,1];
                Append( word2, word );
                SetConjugate( col, j, i, word2 );
            fi;

            # also inverse
            mat := Comm( u2,(u1^-1));
            #Display( mat );
            word := UnipoMat2Word( mat, O, noD, noU );
            #Print( "i ", i, " j ", j , " word ", word, "\n");
            if not word = fail then
                word2 := [j,1];
                Append( word2, word ); 
                SetConjugate( col, j, i, word2 );
            fi;
        od;
    od;
end;



PresentTriang := function( dim, pol )
    local F,O,U_O,isoU_O,n,r,noD,noU,col,D,U,gensTr;

    # setup
    F :=  FieldByPolynomial( pol );
    O := MaximalOrderBasis( F );
    U_O := UnitGroup( F );
    #isoU_O := IsomorphismPcpGroup( U_O );

    # dimension of maximal order
    n := Length( O );
    # number fundamental units plus 1 (torsion unit)
    r := Length( GeneratorsOfGroup( U_O ) );

    # number of gens of D and U
    noD := r*dim;
    noU := n*((dim-1)*dim)/2;
     
    # get collector
    col := FromTheLeftCollector( noD + noU );

    # determine relations wihtin the elements of D
    DetermineRelationsD( col, U_O, dim );

    D := ConstructD( dim, U_O , F);
    U := ConstructU( dim, O, F );

    gensTr := [];
    Append( gensTr, D );
    Append( gensTr, U );

    # determine the relations of D action on U
    DetermineRelationsDonU( gensTr, noD,noU,col,O );
    # determine conjugation relations in U
    DetermineRelationsU( gensTr, noD,noU,col,O );
 
    UpdatePolycyclicCollector( col );
    return rec( Tr := PcpGroupByCollector( col ), noD := noD );

end;


#############################################################################
#
ExamplesOfSomeTGroups := function()
    local G, FF, F,n;
    n := 10;
    G := List( [1..n], x-> ExamplesOfSomePcpGroups(x) );
    FF := List( G, FittingSubgroup );
    F := List( [1..n], x-> PcpGroupByPcp( Pcp( FF[x] ) ) );
    List( F, IsNilpotent );
    return rec( G := G, FF := FF, F := F );

end;



#############################################################################
#
# IN F_nc ..... free nilpotent group on n gens and class c
#
SC_GetSomeAutomorphsimOfF_nc := function( F_nc, n ) 
    local F,A,gensA,gensN,k,auts,ll,gensF,imgs,g,beta,i,aut;

    # get automorphsim induced from free group
    F := FreeGroup( n );
    gensF := GeneratorsOfGroup( F );
    A := AutomorphismGroup( F );
    gensA := GeneratorsOfGroup( A );
    gensN := Pcp( F_nc );
    k := Length( gensA );
    auts := [];
    for i in [1..k] do
        aut := gensA[i];
        ll := List( gensF, x->x^aut );
        imgs := List( ll, x -> MappedWord( x, gensF, gensN{[1..n]} ) );
        Add( auts, GroupHomomorphismByImagesNC(F_nc,F_nc,gensN{[1..n]},imgs ));
    od;

    # get  automorphism of Bryant
    g := gensN[1]*Comm( Comm( gensN[1], gensN[2] ), gensN[1] );
    beta := GroupHomomorphismByImagesNC( F_nc,F_nc,gensN{[1,2]}, [g,gensN[2]]);
    Add( auts, beta );

    return auts;
end;


#############################################################################
# IN: N ..... T-group given by a finite presentation
#     auts... List Automorphisms of N, where <auts> is abelian
# 
# OUT: polycyclic presentation of the semidirect product H \rhd N
#      where H is a free abelian group with Length(auts) generators, which
#      acts like auts on N
#
SC_TGroupByAbelianGroup := function( N, auts )
    local n,m,l,coll,pcpN,exp1,i,j,con,exp,exp2,genList,g;
    
    #get collector
    n := HirschLength( N );
    m := Length( auts );
    l := m+n;
    coll := FromTheLeftCollector( l );
    
    # define conjugation relations within N,   g_i^g_j
    pcpN := Pcp( N );
    exp1 := List( [1..m], x-> 0 );
    for i in [2..n] do 
        for j in [1..(i-1)] do
            # compute g_i^g_j
            con := pcpN[i]^pcpN[j];
            exp2 := Exponents( con ); 
            exp := Concatenation( exp1, exp2 );
            genList:=POL_Exp2GenList(exp);
            SetConjugate( coll, i+m, j+m, genList);

            # compute g_i^(g_j^-1)
            con := pcpN[i]^(pcpN[j]^-1);
            exp2 := Exponents( con ); 
            exp := Concatenation( exp1, exp2 );
            genList:=POL_Exp2GenList(exp);
            SetConjugate( coll, i+m, -(j+m), genList);
        od;
    od;

    #define relations coming from the action of H on N
    for j in [1..m] do
        for i in [1..n] do
            #action of g_j on g_{i+m}
            g := pcpN[i]^(auts[j]);
            exp2 := Exponents( g );
            exp := Concatenation( exp1, exp2 );
            genList:=POL_Exp2GenList(exp);
            SetConjugate( coll, i+m, j, genList);

            #action of g_j^-1 onf g_{i+m}
            g := pcpN[i]^(auts[j]^-1);
            exp2 := Exponents( g );
            exp := Concatenation( exp1, exp2 );
            genList:=POL_Exp2GenList(exp);
            SetConjugate( coll, i+m, -j, genList);
        od;
    od;

    UpdatePolycyclicCollector( coll);
    return PcpGroupByCollector( coll );
    
end;



SC_Exams_Help1 := function( R )
    local G,m,ll,N;
    G := R.Tr;
    m := Length( Pcp(G) );
    ll := [R.noD+1..m];
    N := Subgroup( G, Pcp(G){ll} );     
    return rec( G := G, N := N );
end;

#TODO: look at examples coming from Polenta.

#############################################################################
#
# Examples of some nilpotent-by-abelian groups
# 
#
SC_Exams := function( n )
    local x,pol,R,G,m,ll,N,n1,n2,N2,F,F_nc,hl,auts,T,a;

    x := Indeterminate( Rationals );

    # first 8 groups as in paper
    if n=1 then 
        pol := x^2-3;
        R := PresentTriang( 3, pol );
        return SC_Exams_Help1( R );
    fi;
    if n=2 then 
        pol := x^2-3;
        R := PresentTriang( 4, pol );
        return SC_Exams_Help1( R );
    fi;
    if n=3 then
        pol := x^2-3;
        R := PresentTriang( 5, pol );
        return SC_Exams_Help1( R );
    fi;
    if n=4 then
        pol :=  x^3 - x^2 + 4; 
        R := PresentTriang( 3, pol );
        return SC_Exams_Help1( R );
    fi;
    if n=5 then
        pol :=  x^3 - x^2 + 4; 
        R := PresentTriang( 4, pol );
        return SC_Exams_Help1( R );
    fi;
    #if n=6 then
    #    pol :=  x^3 - x^2 + 4; 
    #    R := PresentTriang( 5, pol );
    #    return SC_Exams_Help1( R );
    #fi;
    if n=6 then
        m := 2;
        F := FreeGroup( m );
        F_nc := NilpotentQuotient( F, 4 );
        hl := HirschLength( F_nc );
        auts := SC_GetSomeAutomorphsimOfF_nc( F_nc, m );
        a := auts[1]*auts[2]*auts[3]^3;
        G := SC_TGroupByAbelianGroup( F_nc, [a] );
        N := Subgroup( G, Pcp(G){[2..hl+1]} );
        T := rec( G := G, N := N );
        return T;
    fi;
    if n=7 then
        m := 2;
        F := FreeGroup( m );
        F_nc := NilpotentQuotient( F, 5 );
        hl := HirschLength( F_nc );
        auts := SC_GetSomeAutomorphsimOfF_nc( F_nc, m );
        a := auts[1]*auts[2]*auts[3]^3;
        G := SC_TGroupByAbelianGroup( F_nc, [a] );
        N := Subgroup( G, Pcp(G){[2..hl+1]} );
        T := rec( G := G, N := N );
        return T;
    fi;
    if n=8 then
        m := 3;
        F := FreeGroup( m );
        F_nc := NilpotentQuotient( F, 4 );
        hl := HirschLength( F_nc );
        auts := SC_GetSomeAutomorphsimOfF_nc( F_nc, m );
        a := auts[1]*auts[3]^3;
        G := SC_TGroupByAbelianGroup( F_nc, [a] );
        N := Subgroup( G, Pcp(G){[2..hl+1]} );
        T := rec( G := G, N := N );
        return T;
    fi;
    ############################################################33

    if n=11 then
        pol :=  x^3 - x^2 + 4; 
        R := PresentTriang( 4, pol );
        return SC_Exams_Help1( R );
    fi;
    if n=12 then 
        pol :=  x^3 - x^2 + 4; 
        R := PresentTriang( 3, pol );
        G := R.Tr;
        m := Length( Pcp(G) );
        ll := [R.noD+1..m];
        N := Subgroup( G, Pcp(G){ll} );     
        n1 := Pcp(N)[2]*Pcp(N)[4];
        n2 := Pcp(N)[5]*Pcp(N)[1]^-1;
        N2 := NormalClosure( G, Subgroup( G, [n1,n2]) );
        return rec( G := G, N := N2 );
    fi;
    if n=13 then
        pol := x^7-x^6+x^5-x^4-3*x^3+3*x^2-2*x+1;
        R := PresentTriang( 3, pol );
        return SC_Exams_Help1( R );
    fi;
    if n=14 then
        pol :=  x^3 - x^2 + 4;
        R := PresentTriang( 4, pol );
        return SC_Exams_Help1( R );
    fi;
    if n=15 then 
        pol :=  x^3 - x^2 + 4; 
        R := PresentTriang( 4, pol );
        G := R.Tr;
        m := Length( Pcp(G) );
        ll := [R.noD+1..m];
        N := Subgroup( G, Pcp(G){ll} );     
        n1 := Pcp(N)[2]*Pcp(N)[5]*Pcp(N)[13]^-1;
        n2 := Pcp(N)[3]*Pcp(N)[7]*Pcp(N)[11]^-1;
        N2 := NormalClosure( G, Subgroup( G, [n1,n2]) );
        return rec( G := G, N := N2 );
    fi;
    if n = 11 then
        pol := x^2 -3;
        R := PresentTriang( 5, pol );
        return SC_Exams_Help1( R );
    fi;
    if n = 21 then 
        m := 2;
        F := FreeGroup( m );
        F_nc := NilpotentQuotient( F, 5 );
        hl := HirschLength( F_nc );
        auts := SC_GetSomeAutomorphsimOfF_nc( F_nc, 2 );
        a := auts[1]*auts[2]*auts[3]^3;
        G := SC_TGroupByAbelianGroup( F_nc, [a] );
        N := Subgroup( G, Pcp(G){[2..hl+1]} );
        T := rec( G := G, N := N );
        return T;
    fi;
    if n = 29 then 
        # the calculation of the presentation needs 3000 sec.
        # Dimension = 34
        F := FreeGroup( 2 );
        return NilpotentQuotient( F, 6 );
    fi;
    if n in [-10..-1] then 
        G := ExamplesOfSomePcpGroups(-n);
        N := FittingSubgroup( G );
        return rec( G := G, N := N );
    fi;
    if n = -11 then
       G := ExamplesOfSomePcpGroups(13);
       N := FittingSubgroup( G );
       return rec( G := G, N := N );
    fi; 
end;




#g1 := SC_RandomPcpElement( T.G, 10 );
#g2 := SC_RandomPcpElement( T.G, 10 );
#g1 := SC_RandomPcpElement( FCR.GG, 10 );
#g2 := SC_RandomPcpElement( FCR.GG, 10 );


# T := SC_Exams( 2 );
#FCR := MAT_FastConjugationRec( T.G, T.N );

#n := Random( FCR.NN );
#n := SC_RandomPcpElement( FCR.NN, 10 );

#g := FCR.factorGens[1];
#c := g^100;

#n^c;
#MAT_FastConjugation( FCR, n, c );
#ExponentsByPcp( FCR.pcpNN, n^c );

SC_Exams_Help2 := function( dim, pol )
    local R,G,N,pcp,N2,phi,m,ll;
    R := PresentTriang( dim, pol );
    G := R.Tr;
    m := Length( Pcp(G) );
    ll := [R.noD+1..m];
    N := Subgroup( G, Pcp(G){ll} );
    pcp := Pcp( N );
    N2 := PcpGroupByPcp( pcp );
    IsNilpotent( N2 );
    return N2;
end;







#############################################################################
#
# IN:  n_fac ..... no. of gens of G/N
#      n_N   ..... no. of gens of N
#      aver ...... bound on absolut value of exponents 
#                  (measure how hard the group is )  
#
# OUT: Produces a random pc presentation of a group where the
#      tail of the generators generates a normal T-group
#      G > N > 1
SC_ProduceRandomPcGroup := function( n_fac, n_N, aver1, aver2 )
    local n,coll,l1,l2,i,j,exp1,exp2,exp3,exp4,exp,genList;
    
    # get collector
    n := n_fac + n_N;
    coll := FromTheLeftCollector( n );
        
    # produce list which is the domain of random exponents
    l1 := [-aver1..aver1];
    l2 := [-aver2..aver2];

    # define relations within nilpotent N   g_i^g_j
    exp1 := List( [1..n_fac], x-> 0 );
    for i in [2..n_N] do
        for j in [1..(i-1)] do
            exp2 := List( [1..(i-1)], x->0 );
            exp3 := [1];
            exp4 := List( [i+1..n_N], x->RandomList( l1 ) ); 
            exp := Concatenation( exp1, exp2, exp3, exp4 );
            genList:=POL_Exp2GenList(exp);
            SetConjugate( coll, i+n_fac, j+n_fac, genList);
        od;
    od;

    # define relations coming from the action of G/N on N
    #for j in [1..n_fac] do
    #    for i in [1..n_N] do
    #        #action of g_j on g_{i+n_fac}\
    #        exp2 := List( [1..n_N], x-> RandomList( l2 ) );
    #        exp := Concatenation( exp1, exp2 );
    #        genList:=POL_Exp2GenList(exp);
    #        SetConjugate( coll, i+n_fac, j, genList);
    #    od;
    #od;

    # define relations coming from the gens within G/N
     

    UpdatePolycyclicCollector( coll);
    return PcpGroupByCollector( coll );
end;


#############################################################################
# Code for the computation of automorphism of the free nilpotent group


DrumHerum := function()
        local n,F,gensF,A,gensA,N,NN,gensNN,aut,ll,imgs,psi,gensN,psi2,g,beta;

 n := 2;
 F := FreeGroup( n );
 gensF := GeneratorsOfGroup( F );
 A := AutomorphismGroup( F );
 gensA := GeneratorsOfGroup( A );
 N := NilpotentQuotient( F, 4 );
 gensN := GeneratorsOfGroup( N );
 
 #NN := Subgroup( N, gensN{[1..n]} );
 #gensNN := GeneratorsOfGroup( NN );
  
 #choose automorphism of F
 aut := gensA[2];
 ll := List( gensF, x->x^aut );
 imgs := List( ll, x -> MappedWord( x, gensF, gensN{[1,2]} ) );
 
 #psi := GroupHomomorphismByImagesNC( NN, NN, gensNN, imgs );
 # geht noch einfacher !
 #GroupHomomorphismByImagesNC( N, N, gensNN, imgs );
 #GroupHomomorphismByImagesNC( N, N, gensN{[1,2]}, imgs );

 
 psi2 := GroupHomomorphismByImagesNC( N, N, gensN{[1,2]}, imgs );


 g := gensN[1]*Comm( Comm( gensN[1], gensN[2] ), gensN[1] );
 beta := GroupHomomorphismByImagesNC( N, N, gensN{[1,2]}, [g, gensN[2]] );

end;





SC_OrbitStabilizerOnRightMod := function (G,erzs, x, length)
        local erzsF,n,k, bahn, stab, tran,tranF, g,gF,F,i, y, w, j;
                                                                              
        #erzs := ShallowCopy (GeneratorsOfGroup (G));
        n := Length( erzs );
        F := FreeGroup( n );
        erzsF :=  ShallowCopy (GeneratorsOfGroup( F ));
        Append (erzs, List (erzs, x -> x^-1));
        Append (erzsF, List (erzsF, x -> x^-1));
        bahn  := [x];
        stab  := [];
        tran  := [One (G)];
        tranF := [One(F)];
        i := 1;
        while i <= Length (bahn) and Length(bahn) <= length do
                y := bahn[i];
                #for g in erzs do
                for k in [1..Length(erzs)] do
                        g :=  erzs[k];
                        gF := erzsF[k];
                        w := y*g;
                        j := Position (bahn, w);
                        if j = fail then
                                Add (bahn, w);
                                Add (tran, tran[i]*g);
                                Add (tranF, tranF[i]*gF);
                                                                               
                        else
                                AddSet (stab, tran[i]*g/tran[j]);
                        fi;
                od;
                i := i + 1;
        od;
        return rec (bahn := bahn, stab := Subgroup (G, stab),
                    tran := tran, tranF := tranF);
end;
                                  
# F is a free group of rank k
# G is a free nilpotent group of rank
# x is an element of G which I want to express a word in the original
# generators                                         
SC_ExpressAsWord := function( F, G, x, length )
    local b,i,k,gens;
    k := Rank( F );
    gens := GeneratorsOfGroup( G ){[1..k]};
    b := SC_OrbitStabilizerOnRightMod(G, gens, One(G), length);
    i := Position( b.bahn, x );
    if i = fail then
        return fail;
    else
        return b.tranF[i];
    fi;
end;

SC_GetOrbitInfo := function( F, G, length )
    local b,i,k,gens;
    k := Rank( F );
    gens := GeneratorsOfGroup( G ){[1..k]};
    b := SC_OrbitStabilizerOnRightMod(G, gens, One(G), length);
    return b;
end;

SC_ExpressAsWord2 := function( F, G, x, orbit )
    local b,i,k,gens;
    b := orbit;
    i := Position( b.bahn, x );
    if i = fail then
        return fail;
    else
        return b.tranF[i];
    fi;
end;
