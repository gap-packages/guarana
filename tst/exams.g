
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
