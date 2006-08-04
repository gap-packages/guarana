#############################################################################
##
#F  Substitute( pol, indets, values )
##
Substitute := Value;
                               
MySubstitute := function( f, indets, values )
    local   n,  one,  null,
            g,
            i,  j,
            pos,
            monom,  vmonom;
    
    if Length( indets ) <> Length( values ) then
        Error( "number of indeterminates must be the number of values" );
    fi;
    
    n := Length( indets );
    
    # if there is nothing to do, do nothing.
    if n = 0 then
        return f;
    fi;
    
    if  f = f * 0 then return 0; fi;

    # prepare lists for easy access
    indets := List( indets, x -> ExtRepPolynomialRatFun( x )[1][1] );
    values := ShallowCopy( values );
    SortParallel( indets, values );
    
    one  := One( values[1] );
    null := 0 * one;
    
    # unwrap f
    f := ExtRepPolynomialRatFun(f);

    
    # run through the monomials
    g := 0 * one;
    for i in [1,3..Length(f)-1] do
        monom  := f[i];
        vmonom := one;
        for j in [1,3..Length(monom)-1 ] do
            pos := PositionSorted( indets, monom[j] );
            if pos > n then
                Error( "cannot substitute for indeterminate ", monom[j] );
            fi;
            if values[ pos ] = null then
                vmonom := null;
                break;
            fi;

            vmonom := vmonom * values[ pos ]^monom[ j+1 ];
        od;
        if vmonom <> null then
            g := g + f[i+1] * vmonom;
        fi;
    od;
    return g;
end;
