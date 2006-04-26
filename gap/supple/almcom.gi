#############################################################################
##
#W almcom.gi              GUARANA package                     Bjoern Assmann
##
## Implementation of Algorithm Almost Complement, Eick habil p. 81
## given Group E polycyclic, M free abelian normal subgroup
##
#H  @(#)$Id$
##
#Y 2006
##
##
## TODO: -checke funktionen auf Extremfaelle

GUARANA.Supple_GetSys := function( CR )
    local sys,c,w;
    # add equations for relators
    sys := CRSystem( CR.dim, Length(CR.mats), CR.char );
    for c in CR.enumrels do
        w := CollectedRelatorCR( CR, c[1], c[2] );
        if IsBound( CR.extension) then
            AddEquationsCR( sys, w[1], w[2], false );
        else
            AddEquationsCR( sys, w[1], w[2], true );
        fi;
    od;
    return sys;
end;

# diese funktion ist noch nicht durchdacht bzgl der Extremfaelle,
# abgeschaut von KernelCR
GUARANA.Supple_GetBeta2AndSolution := function( CR )
    local mat, null, sol, vec,sys;

    # compute system (related to second cohomology map alpha_2)
    sys := GUARANA.Supple_GetSys( CR );

    # case that subgroup is full group, i.e. Length(CR.mats)=0
    if sys.len = 0 then return rec( mat := [], sol := [] ); fi;

    # get beta2
    mat := TransposedMat( sys.base );

    if Length( sys.base ) = 0  then
        sol := List( [1..sys.dim * sys.len], x -> 0 );
    else
        vec := Concatenation( CR.extension );
        sol := SolutionMat( mat, vec );
    fi;

    return rec( mat := mat, sol := sol );
end;

GUARANA.Supple_GetBeta2 := function( CR )
    local sys, c, w,beta2;
  
    # add equations for relators
    sys := CRSystem( CR.dim, Length(CR.mats), CR.char );
    for c in CR.enumrels do
        w := CollectedRelatorCR( CR, c[1], c[2] );
        if IsBound( CR.extension) then
            AddEquationsCR( sys, w[1], w[2], false );
        else
            AddEquationsCR( sys, w[1], w[2], true );
        fi;
    od;

    beta2 := TransposedMat( sys.base );
    return beta2;
end;

# operation on M_Q/M
GUARANA.Supple_OperOnMQFactorM := function( pp, act )
    local p,i,ii,mat,t;
    mat := act[1];
    t := act[2];
    p := Representative( pp );
    i := p * mat + t;
    ii := i + Integers^(Length(p));
    return ii;
end;

#############################################################################
##
#F GUARANA.Supple_ComputeU( CR, s )
##
## IN
## CR ... cohomology record for an extension of M by E/M=G
## s  ... elment of C^1(G,M_Q) which corresponds to the split extension
##        M_Q by G, s is given as vector of Length dim(M)*numberGens 
## 
## OUT
## Stab_G( 0 + M ), where (v+M)^g = v*g + g^deriv(s) + M
## 
## Comments
## 
GUARANA.Supple_ComputeU := function( CR, s )
     local dim,p,pp,n,delta,acts,act,rr,i;

     # get trivial element of M_Q/M 
     dim := CR.dim;
     p := List( [1..dim], x->0 );
     pp := p + Integers^dim;
  
     # get information how E/M acts on M_Q/M
     n := Length( CR.mats );
     # cut s in n pieces
     delta := CutVector( s, n ); 
     acts := [];
     for i in [1..n] do
         act := [CR.mats[i], delta[i]];
         Add( acts, act );
     od;

     rr := PcpOrbitStabilizer( pp, CR.factor, acts, GUARANA.Supple_OperOnMQFactorM );
    
     return rr;
end;

#############################################################################
##
#F GUARANA.Supple_ImageDerivation( der, w, CR )
##
## IN: 
## der ...  derivation on G = < g_1,...,g_n> given by a vector 
##          (s_1,...,s_n) where s_i = g_i^der.
## w  ...    word in g_1,...,g_n.
## CR ...   corresponing cohomology record.
##
## OUT
## vector which represents w(g)^der in M.
##
## Comments:
## This function can be certainly more effective
GUARANA.Supple_ImageDerivation := function( der, w, CR )
    local n,e,i,sum1,sum2,w1,w2;
    
    n := Length( w );
    #exponent and index of generator
    e := w[n][2];
    i := w[n][1];
    if n = 1 then
        # case word corresponds to power of generator       
        if e = 1 then
            return der[i];
        # case g^n
        elif e > 0 then
            sum1 := GUARANA.Supple_ImageDerivation( der, [[i,e-1]], CR )*CR.mats[i];
            sum2 := GUARANA.Supple_ImageDerivation( der, [[i,1]], CR );
            return sum1 + sum2;
        # case g^-n
        else 
            # (h^-1)^der = - (h^der)^(g^-1)
            return -GUARANA.Supple_ImageDerivation( der, [[i,-e]], CR )*(CR.mats[i]^-1);
        fi;
    else 
        w1 := w{[1..n-1]};
        w2 := w{[n]};
        sum1 := GUARANA.Supple_ImageDerivation( der, w1, CR )*(CR.mats[i]^e);
        sum2 := GUARANA.Supple_ImageDerivation( der, w2, CR );
        return sum1 + sum2;
    fi; 

end;

#############################################################################
##
#F GUARANA.Supple_AlmostComplement( E, M )
##
## IN
## E ... polycyclic group
## M ... free abelian normal subgroup
##
## OUT
## An almost complement for M in E. If this does not exist then 'fail'
## is returned.
##
## Comments:
##
GUARANA.Supple_AlmostComplement := function( E, M )
    local CR,s,U,der_s,gensU,i,n,gens,g,exp_img_u,img_u;

    # get cohomology data
    CR := CRRecordBySubgroup( E, M );
    
    # get s \in C^1(G,M_Q) whose image under beta_2 is the 2cocylce,
    # which corresponds to the extension M by E/M
    s := GUARANA.Supple_GetBeta2AndSolution( CR ).sol;
    if s = fail then
        Print( "Extension is not almost split\n" );
        return fail;
    fi;
    der_s := CutVector( s, Length( CR.factor ) );

    # compute U
    U := GUARANA.Supple_ComputeU( CR, s );

    # get generators of almost complement
    gensU := U.stab;
    n := Length( gensU );
    gens := [];
    for i in [1..n] do
        # compute image of u under derivation der_s
        exp_img_u := GUARANA.Supple_ImageDerivation( der_s, U.word[i], CR );
        img_u := MappedVector( exp_img_u, CR.normal ); 
        # get new generator of almost complement
        g := gensU[i]* img_u;
        Add( gens, g );
    od;
   
    return Subgroup( E, gens );

end;


