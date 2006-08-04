#############################################################################
##
#F  TrailingMonomial( pol )
##                               
TrailingMonomial := function( pol )

    pol := ExtRepPolynomialRatFun( pol );
    if Length( pol ) = 0  then
        return fail;
    fi;
    return pol[1];
end;

PositionOfMonomial := function( extrep, monom )
    local   a,  m,  z;
    
    a := 0;  z := Length( extrep ) / 2 + 1;
    while  a+1 < z do               #  extrep[l] < m  and   m  <= list[h]
        m := Int( (a + z) / 2 );    #  a < m < z
        if MonomialTotalDegreeLess( extrep[ 2 * m - 1 ], monom ) then
            a := m;
        else       
            z := m;
        fi;
    od;
    return 2 * z - 1;
end;

CoefficientOfMonomial := function( p, m )
    local   pos;

    if IsRat( p ) then
        if m = [] then return p; fi;
        return 0;
    fi;

    p := ExtRepPolynomialRatFun( p );
    if p = [] then return 0; fi;

    pos := PositionOfMonomial( p, m );
    if pos < Length(p) and p[ pos ] = m then
        return p[ pos + 1 ];
    else
        return 0;
    fi;
end;
  
SmallestTrailingMonomial := function( M, i )
    local   stm,  stc,  pos,  j,  m,  c,  tmp;
    
    #  find the smallest trailing monomial
    stm := fail;
    for j in [i+1..Length( M )] do
        if M[j] <> 0 * M[j] then
            m := TrailingMonomial( M[j] );
            c := AbsInt( CoefficientOfMonomial( M[j], m ) );
            if stm = fail or 
               (MonomialTotalDegreeLess( m, stm ) or c < stc) then
                stm := m;
                stc := c;
                pos := j;
            fi;
        fi;
    od;

    if stm <> fail then
        tmp := M[ i ]; M[ i ] := M[ pos ]; M[ pos ] := tmp;
    fi;
    return stm;
end;

PolynomialGaussIntegral := function( M )
    local   i,  monom,  c,  done,  j,  d;
    
    M := ShallowCopy( M );
    for i in [1..Length( M )] do
        repeat 
            monom := SmallestTrailingMonomial( M, i );
            if monom = fail then 
                break;
            fi;
            c     := CoefficientOfMonomial( M[i], monom );
            done := true;
            for j in [i+1..Length(M)] do
                d    := CoefficientOfMonomial( M[j], monom );
                if d mod c <> 0 then
                    done := false;
                fi;
                M[j] := M[j] - Int( d/c ) * M[i];
            od;
        until done;
    od;    

    i := 1; while i <= Length(M) and M[i] <> 0 * M[i] do
        i := i + 1;
    od;
    
    return M{[1..i-1]};
end;        

PolynomialGaussRational := function( M )
    local   i,  monom,  c,  j,  d;
    
    M := ShallowCopy( M );
    for i in [1..Length( M )] do
        monom := SmallestTrailingMonomial( M, i );
        if monom = fail then break; fi;
        c    := CoefficientOfMonomial( M[i], monom );
        M[i] := M[i] / c;

        for j in [1..i-1] do
            c    := CoefficientOfMonomial( M[j], monom );
            M[j] := M[j] - c * M[i];
        od;
        for j in [i+1..Length(M)] do
            c    := CoefficientOfMonomial( M[j], monom );
            M[j] := M[j] - c * M[i];
        od;
    od;    

    i := 1; while i <= Length(M) and M[i] <> 0 * M[i] do
        i := i + 1;
    od;
    
    return M{[1..i-1]};
end;        

ReducePolynomial := function( M, p )
    local   coeffs,  i,  monom,  c;

    coeffs := [];
    for i in [1..Length( M )] do
        monom := TrailingMonomial( M[i] );

        c := CoefficientOfMonomial( p, monom ) 
                / CoefficientOfMonomial( M[i], monom );
        Add( coeffs, c );
        p := p - c * M[i];
    od;
    return [coeffs, p];
end;

InsertPolynomial := function( M, p )
    local   i;

    p := ReducePolynomial( M, p ); p := p[2];
    if p <> 0 * p then
        i := 1; while i <= Length(M) and MonomialTotalDegreeLess( 
                TrailingMonomial( M[i] ), TrailingMonomial( p ) ) do
            i := i + 1;
        od;
        M{[i+1..Length(M)+1]} := M{[i..Length(M)]};
        M[i] := p;
    fi;

    return p;
end;

ReducePolynomialIntegral := function( M, p )
    local   coeffs,  i,  monom,  rep,  c;

    coeffs := [];
    for i in [1..Length( M )] do
        monom := TrailingMonomial( M[i] );

        c := CoefficientOfMonomial( p, monom ) 
                / CoefficientOfMonomial( M[i], monom );
        
        if not IsInt( c ) then
            M[i] := M[i] / DenominatorRat( c );
            c := NumeratorRat( c );
        fi;

        Add( coeffs, c );
        p := p - c * M[i];
    od;
    Print( coeffs, "\n" );
    return [coeffs, p];
end;

