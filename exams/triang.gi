#############################################################################
##
#W triang.gi              GUARANA package                     Bjoern Assmann
##
## Code for constructing polycyclic presentations of triangular 
## matrix groups over the maximal
## order of a number field. 
##
#H  @(#)$Id$
##
#Y 2006
##
##
##

#############################################################################
##
##
GUARANA.Triang_FindNonZeroEntry := function ( mat )
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

#############################################################################
##
##
GUARANA.Triang_DetermineRelationsD := function( col, U_O, dim )
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

#############################################################################
##
##
GUARANA.Triang_PositionsUnipotentGens := function( pos, n, dim, noD  )
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

#############################################################################
##
##
GUARANA.Triang_UnipoMat2Word := function( mat, O, noD, noU )
    # only for unipotent matrices with one non-zero entry
    local n,dim, ent, coeff,pos,ll,word,i;

    dim := Length( mat );
    # dimension of maximal order
    n := Length( O );

    ent :=  GUARANA.Triang_FindNonZeroEntry( mat );
    if ent = fail then return fail; fi;
    coeff := Coefficients( O, ent[2] );     

    pos := ent[1];
    # get the positions of the corresponding generators in the U(n,O)
    ll := GUARANA.Triang_PositionsUnipotentGens( pos, n, dim, noD  );
    
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

#############################################################################
##
##
GUARANA.Triang_ConstructD := function( dim, U_O , F)
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

#############################################################################
##
##
GUARANA.Triang_ConstructU := function( dim, O, F )
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

#############################################################################
##
##
GUARANA.Triang_DetermineRelationsDonU := function( gensTr, noD,noU,col,O )
    local i,d,j,u,mat,word;
    for i in [1..noD] do
        d := gensTr[i];
        for j in [noD+1..noD+noU] do
            u := gensTr[j];
            mat := u^d;
            word := GUARANA.Triang_UnipoMat2Word( mat, O, noD, noU );
            #Print( "i ", i, " j ", j , " word ", word, "\n");
            SetConjugate( col, j, i, word );

            # also inverse
            mat := u^(d^-1);
            word := GUARANA.Triang_UnipoMat2Word( mat, O, noD, noU );
            #Print( "i ", i, " j ", j , " word ", word, "\n");
            SetConjugate( col, j, i, word );
        od;
    od;
end;

#############################################################################
##
##
GUARANA.Triang_DetermineRelationsU := function( gensTr, noD,noU,col,O )
    local i,u1,j,u2,mat,word, word2;
    for i in [noD+1..(noD+noU)] do
        u1 := gensTr[i];
        for j in [i+1..noD+noU] do
            u2 := gensTr[j];
            mat := Comm( u2, u1 );
            #Display( mat );
            word := GUARANA.Triang_UnipoMat2Word( mat, O, noD, noU );
            #Print( "i ", i, " j ", j , " word ", word, "\n");
            if not word = fail then 
                word2 := [j,1];
                Append( word2, word );
                SetConjugate( col, j, i, word2 );
            fi;

            # also inverse
            mat := Comm( u2,(u1^-1));
            #Display( mat );
            word := GUARANA.Triang_UnipoMat2Word( mat, O, noD, noU );
            #Print( "i ", i, " j ", j , " word ", word, "\n");
            if not word = fail then
                word2 := [j,1];
                Append( word2, word ); 
                SetConjugate( col, j, i, word2 );
            fi;
        od;
    od;
end;

#############################################################################
##
#F GUARANA.Triang_PresentTriang( dim, pol )
##
## IN
## dim ............. degree of matrix group.
## pol ............. polynomial over Q that is used to construct a 
##                   number field K.
## 
## OUT 
## A polycyclic presentation of the upper triangular group Tr(dim,K),
## i.e. the group containing all upper triangular matrices of degree
## "dim" over K. 
##
GUARANA.Triang_PresentTriang := function( dim, pol )
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

    # number of gens of D and U.
    # U denotes the upper unitriangular subgroup of Tr(dim,K)
    # and D denotes the subgroup of all diagonal matrices. 
    noD := r*dim;
    noU := n*((dim-1)*dim)/2;
     
    # get collector
    col := FromTheLeftCollector( noD + noU );

    # determine relations wihtin the elements of D
    GUARANA.Triang_DetermineRelationsD( col, U_O, dim );

    D := GUARANA.Triang_ConstructD( dim, U_O , F);
    U := GUARANA.Triang_ConstructU( dim, O, F );

    gensTr := [];
    Append( gensTr, D );
    Append( gensTr, U );

    # determine the relations of D acting on U
    GUARANA.Triang_DetermineRelationsDonU( gensTr, noD,noU,col,O );
    # determine conjugation relations in U
    GUARANA.Triang_DetermineRelationsU( gensTr, noD,noU,col,O );
 
    UpdatePolycyclicCollector( col );
    return rec( Tr := PcpGroupByCollector( col ), noD := noD );
end;


#############################################################################
##
#F GUARANA.Triang_UpperTriangAndUnitriang( R )
##
## IN
## R ............... record as poduced by Triang_PresentTriang.
## 
## OUT
## G a polycyclically presented group ismorphic to 
## the upper triangular group Tr(dim,K),
## i.e. the group containing all upper triangular matrices of degree
## "dim" over K. 
## N a polycyclically presented group, subgroup of G and ismorphic
## to Tr_2(dim,K).
##
## This function was previously called SC_Exams_Help1
##
GUARANA.Triang_UpperTriangAndUnitriang := function( R )
    local G,m,ll,N;
    G := R.Tr;
    m := Length( Pcp(G) );
    ll := [R.noD+1..m];
    N := Subgroup( G, Pcp(G){ll} );     
    return rec( G := G, N := N );
end;


