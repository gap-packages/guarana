BinomialCoeff := function( m, n )
    local   bc,  k;

    bc := 1;
    for k in [1..n] do
        bc := ( bc * (m - k + 1) ) / k;
    od;
    return bc;
end;
    

PolyCoeff := function( arg )
    local   M,  x,  r,  C,  MM,  i;

    M := arg[1];
    if Length(arg) = 1 then
        x := Indeterminate( Rationals, "x" );
    elif IsString( arg[2] ) then
        x := Indeterminate( Rationals, arg[2] );
    else
        x := arg[2];
    fi;

    r := Length( M );
    M := M - IdentityMat( r );

    C := IdentityMat( r );
    MM := M;
    for i in [1..r-1] do 
        C := C + BinomialCoeff( x, i ) * MM; 
        MM := MM * M;
        if MM = MM * 0 then break; fi;
    od;
    return C;
end;

PolyCoeffAsBinomialPol := function( arg )
    local   M,  x,  r,  C,  i;

    M := arg[1];
    if Length(arg) = 1 then
        x := Indeterminate( Rationals, "x" );
    else
        x := Indeterminate( Rationals, arg[2] );
    fi;

    r := Length( M );
    M := M - IdentityMat( r );

    C := IdentityMat( r );
    for i in [1..r-1] do C := C + x^i * M^i; od;
    return C;
end;

SubstitutePolyCoeff := function( M, value )
    local   V,  i,  j;

    V := [];
    for i in [1..Length(M)] do
        V[i] := [];
        for j in [1..Length(M[i])] do
            V[i][j] := Value( M[i][j], value );
        od;
    od;
    return V;
end;

RandomUnipotentUpperTriangularMat := function( n, range )
    local   M,  i,  j;

    M := IdentityMat( n );
    for i in [1..n] do
        for j in [i+1..n] do
            M[i][j] := Random( range );
        od;
    od;
    return M;
end;

##
##  PolynomialCoefficientsMatrixRep( phi )
##  PolynomialCoefficientsMatrixRep( phi, indets )
##
PolynomialCoefficientsMatrixRep := function( arg )
    local   phi,  pcp,  d,  indeterminates,  P,  i,  M;

    phi := arg[1];
    pcp := Pcp( Source( phi ) );
    d   := Length( GeneratorsOfGroup( Range( phi ) )[1] );

    if Length( arg ) = 2 then
        indeterminates := arg[2];
    else
        indeterminates := 
          List( [1..Length(pcp)],
                i -> Indeterminate( Rationals,
                        Concatenation( "x", String(i) ) ) );
    fi;

    P := IdentityMat( d );
    for i in [1..Length(pcp)] do
        M := Image( phi, pcp[i] );
        P := P * PolyCoeff( M, indeterminates[i] );
    od;

    return rec( indeterminates := indeterminates,
                coefficients := PolynomialGaussRational( Flat(P) ) );
end;

PolynomialCoefficientsMatrices := function( matrices, indeterminates )
    local   P,  i;

    P := matrices[1]^0;
    for i in [1..Length(matrices)] do
        P := P * PolyCoeff( matrices[i], indeterminates[i] );
    od;

    return rec( indeterminates := indeterminates,
                coefficients := PolynomialGaussRational( Flat(P) ) );
end;

MatricesCoeffsGens := function( coeffs, genpols )
    local   d,  t,  matrices,  i,  matrix,  c,  v,  r;

    d := Length( genpols );
    t := coeffs.indeterminates[1];

    Print( "#  Compute matrices\n" );
    matrices := [];
    for i in [1..d] do
        matrix := [];
        for c in coeffs.coefficients do
            v := Value( c, coeffs.indeterminates{[i..d]}, genpols[i] );
            if IsCyclotomic( v ) then v := v + 0 * t; fi;
                
            r := ReducePolynomial( coeffs.coefficients, v );
            if r[2] * 0 <> r[2] then
                Error();
            fi;
            Add( matrix, r[1] );
        od;
        Add( matrices, TransposedMat(matrix) );
    od;

    return matrices;
end;

MatrixRepresentationDTPols := function( obj )
    local   G,  mpols,  d,  phi,  coeffs,  genpols,  y,  i,  matrices,  
            matrix,  c,  v,  r;

    if IsGroup( obj ) then
        G := obj;
        phi := UnitriangularMatrixRepresentation( G );
    else
        phi := obj;
        G := Source( phi );
    fi;

    Print( "#  Hall multiplication polynomials\n" );
    mpols := DTPolynomialsByPcpGroup( G );
    d     := Length( mpols.indeterminates ) - 1;

    Print( "#  coefficients\n" );
    coeffs := PolynomialCoefficientsMatrixRep( phi, 
                      mpols.indeterminates{[1..d]} );

    Print( "#  Replace y\n" );
    genpols := [];
    y := mpols.indeterminates[d+1];
    for i in [1..d] do
        genpols[ i ] := List( mpols.polynomials[i], 
                               p->Value( p, [y], [1] ) );
    od;

    return MatricesCoeffsGens( coeffs, genpols );
end;

ExtendMatRepresentation := function( G, mpols, genpols, matrices )
    local   d,  k,  coeffs,  t,  queue,  c,  v,  r,  i,  matrix;

    d := Length( GeneratorsOfGroup( G ) );
    k := d-Length(matrices)+1;
    coeffs := PolynomialCoefficientsMatrices( matrices,
                      mpols.indeterminates{[k..d]} );

    t := mpols.indeterminates[d+1];
    queue := ShallowCopy( coeffs.coefficients );
    while Length( queue ) > 0 do

        c := queue[ Length(queue) ]; Unbind( queue[ Length(queue) ] );
        v := Value( c, mpols.indeterminates{[k..d]}, 
                    genpols[k-1]{[2..d-k+2]} );
        if IsCyclotomic( v ) then v := v + 0 * t; fi;
        
        r := InsertPolynomial( coeffs.coefficients, v );
        if r * 0 <> r then
            Add( queue, r );
        fi;
        
        for i in [k..d] do 
            v := Value( c, mpols.indeterminates{[i..d]}, genpols[i] );
            if IsCyclotomic( v ) then v := v + 0 * t; fi;
    
            r := InsertPolynomial( coeffs.coefficients, v );
            if r * 0 <> r then
                Add( queue, r );
            fi;
        od;
    od;

    matrices := [];
    matrix := [];
    for c in coeffs.coefficients do
        v := Value( c, mpols.indeterminates{[k..d]}, 
                    genpols[k-1]{[2..d-k+2]} );
        if IsCyclotomic( v ) then v := v + 0 * t; fi;
        
        r := ReducePolynomial( coeffs.coefficients, v );
        if r[2] * 0 <> r[2] then Error(); fi;
        Add( matrix, r[1] );
    od;
    Add( matrices, TransposedMat(matrix) );
    for i in [k..d] do
        matrix := [];
        for c in coeffs.coefficients do
            v := Value( c, mpols.indeterminates{[i..d]}, genpols[i] );
            if IsCyclotomic( v ) then v := v + 0 * t; fi;
                
            r := ReducePolynomial( coeffs.coefficients, v );
            if r[2] * 0 <> r[2] then Error(); fi;
            Add( matrix, r[1] );
        od;
        Add( matrices, TransposedMat(matrix) );
    od;

    return matrices;
end;

ExtendDirectSum := function( matrices )
    local   d,  em,  extended,  m;
    
    d := Length( matrices[1] );

    em := IdentityMat( d+2 );
    em[d+1][d+2] := 1;

    extended := [ em ];
    for m in matrices do
        em := IdentityMat( d+2 );
        em{[1..d]}{[1..d]} := m;
        Add( extended, em );
    od;
    return extended;
end;

BuildMatRepresentation := function( G )
    local   t,  mpols,  d,  y,  genpols,  matrices,  i;

    t := Runtime();
    Print( "Computing Hall polynomials \c" );
    mpols := DTPolynomialsByPcpGroup( G );
    Print( "(", Runtime()-t, " msec)\n" );

    t := Runtime();
    Print( "Computing representation \c" );
    d := Length( mpols.indeterminates ) - 1;
    y := mpols.indeterminates[d+1];
    genpols := List( [1..d], i->List( mpols.polynomials[i], 
                       p->Value( p, [y], [1] ) ) );;

    matrices := [  [[ 1, 0 ], [ 1, 1 ]] ];

    for i in [d-1,d-2..1] do
        #Print( "Step ", i, ": \c" );
        if genpols[i]{[2..d-i+1]} = mpols.indeterminates{[i+1..d]} then
            matrices := ExtendDirectSum( matrices );
        else
            matrices := ExtendMatRepresentation( G,mpols,genpols,matrices );
        fi;
        #Print( " dim = ", Length(matrices[1]),
        #       " (", Runtime()-t, " msec)\n" ); t := Runtime();
    od;

    return matrices;
end;

ExtendDualRepresentationSlow := function( G, mpols, genpols, k, module )
    local   d,  t,  queue,  c,  v,  r,  i;

    d := Length( mpols.indeterminates ) - 1;

    t := mpols.indeterminates[d+1];
    queue := ShallowCopy( module.functions );
    while Length( queue ) > 0 do

        c := queue[ Length(queue) ]; Unbind( queue[ Length(queue) ] );
        v := Value( c, mpols.indeterminates{[k..d]}, 
                    genpols[k-1]{[2..d-k+2]} );
        if IsCyclotomic( v ) then v := v + 0 * t; fi;
        
        r := InsertPolynomial( module.functions, v );
        if r * 0 <> r then
            Add( queue, r );
        fi;
 
#        for i in [k..d] do 
#            v := Value( c, mpols.indeterminates{[i..d]}, genpols[i] );
#            if IsCyclotomic( v ) then v := v + 0 * t; fi;
#    
#            r := InsertPolynomial( module.functions, v );
#            if r * 0 <> r then
#                Add( queue, r );
#            fi;
#        od;
    od;
    return module;
end;

ExtendDualRepresentation := function( G, mpols, genpols, k, module )
    local   d,  t,  queue,  c,  r;

    d := Length( mpols.indeterminates ) - 1;

    t := mpols.indeterminates[d+1];
    queue := ShallowCopy( module.functions );

    for c in queue do

#        Print( c, "\n" );
        repeat
            c := Value( c, mpols.indeterminates{[k..d]}, 
                        genpols[k-1]{[2..d-k+2]} ) - c;
#            Print( c, "\n" );
            if IsCyclotomic( c ) then c := c + 0 * t; fi;
        
            r := InsertPolynomial( module.functions, c );
        until r * 0 = r;
#        Print( "\n" );
    od;
    return module;
end;

MatRepByModule := function( G, mpols, genpols, module )
    local   d,  t,  matrices,  matrix,  c,  v,  r,  i;

    d := Length( mpols.indeterminates ) - 1;
    t := mpols.indeterminates[ d+1 ];

    matrices := [];
    matrix := [];
    for c in module.functions do
        v := Value( c, mpols.indeterminates{[2..d]}, 
                    genpols[1]{[2..d]} );
        if IsCyclotomic( v ) then v := v + 0 * t; fi;
        
        r := ReducePolynomial( module.functions, v );
        if r[2] * 0 <> r[2] then Error(); fi;
        Add( matrix, r[1] );
    od;
    Add( matrices, TransposedMat(matrix) );
    for i in [2..d] do
        matrix := [];
        for c in module.functions do
            v := Value( c, mpols.indeterminates{[i..d]}, genpols[i] );
            if IsCyclotomic( v ) then v := v + 0 * t; fi;
                
            r := ReducePolynomial( module.functions, v );
            if r[2] * 0 <> r[2] then Error(); fi;
            Add( matrix, r[1] );
        od;
        Add( matrices, TransposedMat(matrix) );
    od;
    return matrices;
end;

BuildDualRepresentation := function( G )
    local   t,  mpols,  d,  y,  genpols,  module,  i,  mrep;

    t := Runtime();
    Print( "Computing Hall polynomials \c" );
    mpols := DTPolynomialsByPcpGroup( G );
    Print( "(", Runtime()-t, " msec)\n" );

    d := Length( mpols.indeterminates ) - 1;
    y := mpols.indeterminates[d+1];
    genpols := List( [1..d], i->List( mpols.polynomials[i], 
                       p->Value( p, [y], [1] ) ) );;

    module := rec( functions := [ 1+y*0 ] );
    InsertPolynomial( module.functions, mpols.indeterminates[d] );

    for i in [d-1,d-2..1] do
        #Print( "Step ", i, ": \c" );    t := Runtime();

        if genpols[i]{[2..d-i+1]} = mpols.indeterminates{[i+1..d]} then
            InsertPolynomial( module.functions, mpols.indeterminates[i] );
        else
            module := ExtendDualRepresentation( G,mpols,genpols,i+1,module );
        fi;
        #Print( " dim = ", Length(module.functions),
        #        " (", Runtime()-t, " msec)\n" ); t := Runtime();
    od;

    Print( "Computing representation " ); t := Runtime();
    mrep := MatRepByModule( G, mpols, genpols, module );
    Print( " (", Runtime()-t, " msec)\n" );

    return [ module, mrep ];
end;

