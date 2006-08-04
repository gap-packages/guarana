BinomialPolynomial := function( x, k )
    local   binom,  i;

    binom := 1;
    for i in [0..k-1] do
        binom := binom * (x-i) / (i+1);
    od;
    return binom;
end;

BinomialProduct := function( x, i, j )
    local   prod,  k;
    
    prod := 0 * x;
    for k in [0..j] do
        prod := prod + Binomial( j, k ) * Binomial( k+i, j ) * x^(k+i);
    od;
    return prod;
        
end;
MonomialDTMonomial := function( monom, n, x )
    local   m,  j;

    m := 1;
    for j in [1,3..Length(monom)-1] do
        if monom[j] = 0 then
            m := m * BinomialPolynomial( x[ n+1 ], monom[ j+1 ] );
        else
            m := m * BinomialPolynomial( x[ monom[j] ], monom[ j+1 ] );
        fi;
    od;
    return m;
end;

PolynomialDTPol := function( pols, i, n, x )
    local   exponents,  j,  monomials,  where,  monom,  k;

    exponents := [];
    exponents := [ x[n+1] + x[i] ];
    for j in [i+1..n] do Add( exponents,  x[j] ); od;

    if not IsInt(pols) then
        monomials := pols.evlist;
        where     := pols.evlistvec;
        
        for j in [1..Length(monomials)] do
            monom := monomials[j];
            monom := MonomialDTMonomial( monom{[5..Length(monom)]}, n, x );
            for k in [1,3..Length(where[j])-1] do
                exponents[ where[j][k]-i+1 ] :=
                  exponents[ where[j][k]-i+1 ] + where[j][k+1] * monom;
            od;
        od;
    fi;

    return exponents;
end;

PolynomialsDTCollector := function( collector )
    local   dtpols,  n,  x,  pols,  i;

    UpdatePolycyclicCollector( collector );
    dtpols := collector![ PC_DEEP_THOUGHT_POLS ];

    n := Length( dtpols );
    x := List( [1..n], i->Indeterminate( Rationals, 
                 Concatenation( "x",String(i) ) : new ) );
    Add( x, Indeterminate( Rationals, "y" : new ) );

    pols := [];
    for i in [1..n] do
        pols[i] := PolynomialDTPol( dtpols[i], i, n, x );
    od;

    return rec( indeterminates := x, polynomials := pols );
end;

##
##  Produce an exponent vector.
##
ExponentVector := function( n, c )

    return List( 
                 List( [1..n], i->Concatenation( c, String(i) ) ),
                 x->Indeterminate( Integers, x : new ) );           
end;

##
##  Multiply the exponent vectors u and z with the help of dtpols.
##
ProductDT := function( arg )
    local   u,           #  left factor
            z,           #  right factor,
            DT,          #  deep thought data
            n,           #  number of generators to be multiplied
            null,
            dtpols,      #  deep thought polynomials
            product,     #  result
            exps,        #  intermediate result during multiplication
            newexps,     #  intermediate result during multiplication
            i,  j;       #  loop variables
    
    u := arg[1];    z := arg[2];    DT := arg[3];
    if Length( arg ) = 4 then
        n := arg[4];
    else
        n := Length( DT.polynomials );
    fi;
    
    dtpols := DT.polynomials;

    null := u[1] * 0;

    # Collect the completed exponents in product.
    product := [];

    # Compute the exponents of the product.
    exps := ShallowCopy(u);
    for i in [1..n] do

        if z[i] <> null then
            newexps := ListWithIdenticalEntries( i-1, 0 );

            # Multiply with z_i from the right.
            Add( exps, z[i] );
            for j in [1..Length(dtpols[i])] do
                Add( newexps, 
                     Substitute( dtpols[i][j], DT.indeterminates, exps ) );
            od;

            Add( product, newexps[i] );
            newexps[i] := 0;
            exps := newexps;
        else
            Add( product, exps[i] );
            exps[i] := 0;
        fi;
            
    od;
    return product;
end;


##
##  Multiply an exponent vector with a generator exponent pair from the right
##  with respect to the given Deep Thought polynomials.
##
ProductExpVecGen := function( u, g, e, DT )
    local   dtpols,  n,  product,  replace,  p;
    
    dtpols := DT.polynomials;
    n := Length( dtpols );

    product := u{[1..g-1]};

    replace := Concatenation( ListWithIdenticalEntries( g-1, 0 ), u{[g..n]} );
    Add( replace, e );

    for p in dtpols[g] do
        Add( product, Substitute( p, DT.indeterminates, replace ) );
    od;
    return product;
end;


##
##  Solve the equation  u x = v with respect to the given Deep Thought
##  polynomials.
##
SolutionDT := function( u, v, DT )
    local   n,  x,  uu,  i,  null;

    # Number of exponents.
    n := Length(DT.polynomials);

    null := 0*u[1];

    # Solve.
    x := []; uu := u;
    for i in [1..n] do
        Add( x, v[i] - uu[i] );
        if x[i] <> null then
            uu := ProductExpVecGen( uu, i, x[i], DT );
        fi;
    od;
    return x;
end;
                
##
##  Invert an exponent vector.
##
InverseDT := function( u, DT )

    return SolutionDT( u, ListWithIdenticalEntries(
                   Length(DT.polynomials), 0 ), DT );
end;


##
##  Compute the commutator of two exponent vectors.
##
CommutatorDTSlow := function( u, v, DT )
    local uv, vu;

    Print( "Computing the product of u and v\n" );
    uv := ProductDT( u,v, DT );
    Print( "Computing the product of v and u\n" );
    vu := ProductDT( v,u, DT );
    Print( "Computing the inverse of vu\n" );
    vu := InverseDT( vu, DT );
    
    Print( "Computing the commutator of v and u\n" );
    return ProductDT( vu, uv, DT );
end;

##
##  Compute the commutator of two exponent vectors.
##
CommutatorDT := function( u, v, DT )
    local uv, vu;

    Print( "Computing the product of u and v\n" );
    uv := ProductDT( u,v, DT );
    Print( "Computing the product of v and u\n" );
    vu := ProductDT( v,u, DT );
    
    Print( "Computing the commutator of v and u\n" );
    return SolutionDT( vu, uv, DT );
end;

##
##  Compute the commutator [u,v] by solving the equation  vu x = uv with
##  respect to the given Deep Thought polynomials.
##
CommutatorSolutionDT := function( u, v, DT )
    local   n,  null,  x,  vl,  ul,  ur,  vr,  i;

    # Number of exponents.
    n := Length(DT.polynomials);

    null := 0*u[1];

    # Solve.
    x := []; 
    vl := ShallowCopy(v); ul := u;
    ur := ShallowCopy(u); vr := v;
    # vl * ul * x = vr * ur
    for i in [1..n] do
        Print( "#  collecting generator ", i, "\n" );
        Add( x, ur[i] + vr[i] - ul[i] - vl[i] );
        # move x past ul
        if x[i] <> null then
             ul := ProductExpVecGen( ul, i, x[i], DT );
        fi;
        # move ul past vl
        if ul[i] <> null then
             vl := ProductExpVecGen( vl, i, ul[i], DT );
        fi;
        # move vr past ur
        if vr[i] <> null then
             ur := ProductExpVecGen( ur, i, vr[i], DT );
        fi;
    od;
    return x;
end;
                


##
##  Compute the n-th Engel word of u and v.
##
EngelDT := function( u, v, n, DT )
    local   i;

    for i in [1..n] do
        u := CommutatorDT( u, v, DT );
    od;
    return u;
end;

ReduceMultivariatePol := function( p, q )
    local   F,  f,  zero,  mons,  cofs,  i,  k,  e,  prd,  c;

    F := FamilyObj( p );
    f := F!.zippedProduct;

    p := ExtRepOfObj( p );
    zero := p[1];
    p    := p[2];
    mons := p{[1,3..Length(p)-1]};
    cofs := p{[2,4..Length(p)]};

    for i in [1..Length(mons)] do
        cofs[i] := cofs[i] mod q;

        for k in [2,4..Length(mons[i])] do
            mons[i] := ShallowCopy(mons[i]);
            e := mons[i][k];
            if e <> 0 then 
                e := e mod (q-1);
                if e = 0 then e := q-1; fi;
            fi;
            mons[i][k] := e;
        od;
    od;

    # sort monomials
    SortParallel( mons, cofs, f[2] );

    # sum coeffs
    prd := [];
    i   := 1;
    while i <= Length(mons)  do
        c := cofs[i];
        while i < Length(mons) and mons[i] = mons[i+1]  do
            i := i+1;
            c := f[3]( c, cofs[i] );
        od;
        if c <> zero  then
            Add( prd, mons[i] );
            Add( prd, c );
        fi;
        i := i+1;
    od;

    return ObjByExtRep( F, [ zero, prd ] );
end;

CoefficientsMultivariatePol := function( p )

    p := ExtRepOfObj( p );
    p := p[2];
    return p{[2,4..Length(p)]};
end;

ConsistencyDT := function( P )
    local   n,  ev,  u,  v,  w,  uv,  uv_w,  vw,  u_vw;
    
    n := Length( P.indeterminates ) - 1;
    u := ExponentVector( n, "c" );
    v := ExponentVector( n, "b" );
    w := ExponentVector( n, "a" );
    
    Print( "Computing u v\n" );
    uv   := ProductDT( u,v,  P );
    Print( "Computing (u v) w\n" );
    uv_w := ProductDT( uv,w, P );
    Print( "Computing v w\n" );
    vw   := ProductDT( v,w,  P );
    Print( "Computing u (v w)\n" );
    u_vw := ProductDT( u,vw, P );
    
    return uv_w - u_vw;
end;


BinomialMonomial := function( fam, m )
    local   bm,  i,  y;
    
    bm := PolynomialByExtRep( fam, [[], 1] );
    for i in [1,3..Length(m)-1] do
        y := PolynomialByExtRep( fam, [ [m[i], 1], 1 ] );
        bm := bm * BinomialPolynomial( y, m[i+1] );
    od;
    return bm;
end;

AsBinomialPolynomial := function( q )
    local   fam,  null,  lm,  lc,  bm,  c,  asBinom;
    
    asBinom := [];
    fam := FamilyObj( q );
    null := 0 * q;
    while q <> null do
        lm := LeadingMonomial( q );
        lc := LeadingCoefficient( q );
        
        bm := BinomialMonomial( fam, lm );
        c := lc / LeadingCoefficient( bm );
        
        Add( asBinom, c );
        Add( asBinom, lm );
        
        q := q - c * bm;
    od;
    
    asBinom := Reversed( asBinom );
    return PolynomialByExtRep( fam, asBinom );
end;

AsOrdinaryPolynomial := function( q )
    local   p,  fam,  i;
    
    p := 0 * q;
    
    fam := FamilyObj( q );
    
    q := Reversed( ExtRepOfObj( q )[2] );
    for i in [1,3..Length(q)-1] do
        p := p + q[i] * BinomialMonomial( fam, q[i+1] );
    od;
    return p;
end;

FindLargestIntegralIndetermiante := function( p )
    local   fam,  ep,  terms,  gcd,  i,  m;
    
    fam := FamilyObj( p );
    p := AsBinomialPolynomial( p );
    ep := ExtRepOfObj( p );
    
    terms := ep[2];
    gcd  := Gcd( terms{[2,4..Length(terms)]} );
    for i in Reversed( [1,3..Length(terms)-1] ) do
        if Length( terms[i] ) = 2 and terms[i][2] = 1 and
           gcd mod terms[i+1] = 0 then
            m := ObjByExtRep( fam, [ ep[1], [terms[i], 1] ] );
            return [ m, terms[i+1] ];
        fi;
    od;
    return fail;
    
end;
  
FindSmallestIndeterminate := function( p )
    local   ep,  terme,  indet,  i,  m,  coeff,  fam;
    
    ep := ExtRepOfObj( p );
    terme := ep[2];
    indet := fail;
    for i in [1,3..Length( terme )-1] do
        m := terme[i];
        if Length( m ) = 2 and m[2] = 1 then
            indet := m; coeff := terme[i+1]; break;    
        fi;
    od;
    if indet = fail then return fail; fi;
    
    fam := FamilyObj( p );
    m   := ObjByExtRep( fam, [ ep[1], [ indet, 1 ] ] );
    
    return [ m, m - p / coeff ];
end;
                     
                     
                     
CentreDT := function( pols )
    local   n,  u,  v,  uv,  vu,  i,  w,  pair;
    
    n := Length( pols.indeterminates ) - 1;
    
    u := ExponentVector( n, "u" );
    v := ExponentVector( n, "v" );
    
    uv := ProductDT( u,v, pols );
    vu := ProductDT( v,u, pols );
    
    for i in [1..n] do
        w := uv[i] - vu[i];
        if w <> 0*w then
            pair := FindSmallestIndeterminate( w );
            if pair = fail then
                Error( "human interaction required: stage ", i, 
                       ", equation: ", w );
            else
                u[ pair[1] ] := pair[2];
            fi;
            uv := ProductDT( u,v, pols );
            vu := ProductDT( v,u, pols );
        fi;
    od;
    
    return u;
end;
            
CentralizerDT := function( v, pols )
    local   n,  u,  uv,  vu,  i,  w,  pair,  pos;
    
    n := Length( pols.indeterminates ) - 1;
    
    u := ExponentVector( n, "u" );
    uv := ProductDT( u,v, pols );
    vu := ProductDT( v,u, pols );
    for i in [1..n] do
        w := uv[i] - vu[i];
        if w <> 0*w then
            pair := FindLargestIntegralIndetermiante( w );
            if pair = fail then
                Error( "human interaction required: stage ", i, 
                       ", equation: ", w );
            else
                pos := Position(u,pair[1]);
                u[ pos ] := u[ pos ] - 1 / pair[2] * w;
            fi;
            uv := ProductDT( u,v, pols );
            vu := ProductDT( v,u, pols );
        fi;
    od;
    
    return u;
end;
            
CentralizerDTInteractive := function( v, pols )
    local   n,  u,  uv,  vu,  i,  w,  pair,  pos;
    
    n := Length( pols.indeterminates ) - 1;
    
    u := ExponentVector( n, "u" );
    uv := ProductDT( u,v, pols );
    vu := ProductDT( v,u, pols );
    for i in [1..n] do
        w := uv[i] - vu[i];
        if w <> 0*w then
            pair := FindLargestIntegralIndetermiante( w );
            pos  := Position(u,pair[1]);
            Error( "\nstage ", i, "\n",
                   "equation: ", w, "\n",
                   "pos: ", pos, "\n",
                   "pair: ", pair, "\n" );
#                u[ pos ] := u[ pos ] - 1 / pair[2] * w;
            uv := ProductDT( u,v, pols );
            vu := ProductDT( v,u, pols );
        fi;
    od;
    
    return u;
end;
            

TruncatePolynomials := function( pols, m )
    local   n,  trcted,  i;
    
    n := Length( pols.indeterminates ) - 1;
    
    if m > n then
        Error( "<n> must be at most the Hirsch length" );
    fi;
    
    trcted := rec();
    trcted.indeterminates      := pols.indeterminates{[1..m]};
    trcted.indeterminates[m+1] := pols.indeterminates[n+1];
    
    trcted.polynomials := [];
    for i in [1..m] do
        trcted.polynomials[i] := pols.polynomials[i]{[i..m]-i+1};
    od;
    return trcted;
end;

IsDivisiblePolynomial := function( p, m )
    local   c;

    p := AsBinomialPolynomial( p );
    p := ExtRepPolynomialRatFun( p );
    for c in p{[2,4..Length(p)]} do
        if c mod m <> 0 then return false; fi;
    od;
    return true;
end;

AbelianizedReducedDT := function( u, pols, DT )
    local   null,  exps,  g,  p,  pow,  h;
    
    null := u[1] * 0;
    exps := DT![PC_EXPONENTS];
    u := ShallowCopy( u );
    for g in [1..Length(u)] do
        if u[g] <> null then
            Print( "reducing exponent ", g, "\n" );
            if IsBound(exps[g]) and exps[g] <> 0 then
                if IsDivisiblePolynomial( u[g], exps[g] ) then
                    Print( "exponent ", g, " is divisible by ", 
                           exps[g], "\n" );
                    
                    if IsBound( DT![PC_POWERS][g] ) then
                        p := u[g] / exps[g];
                        pow := DT![PC_POWERS][g];
                        for h in [1,3..Length(pow)-1] do
                            u[ pow[h] ] := u[ pow[h] ] + pow[h+1] * p;
                        od;
                    fi;
                    u[g] := null;
                else
                    Error( "entry ", u[g], " (", g, ")", 
                           " in exponent vector is not divisible by ",
                           exps[g] );
                fi;
            else
                Error( "entry ", u[g], " (", g, ")", 
                       " in exponent vector is not zero" );
            fi;
        fi;
    od;
    return u;
end;
