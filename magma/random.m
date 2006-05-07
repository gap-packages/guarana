/*###########################################################################
#W random.m                GUARANA package                     Bjoern Assmann
##
## Magma code for testing the runtime of collection from the left.  
##
#H  @(#)$Id$
##
#Y 2006
*/

RuntimeCftl := function( G, range )
    g := Random( G, range );
    h := Random( G, range );
    t := Cputime();
    k := g*h;
    return Cputime( t );
end function;

RuntimesCftl := function( G, range, no )
    res := [];
    for i in [1..no] do 
        t := RuntimeCftl( G, range );
        Append( ~res, t );
    end for;
    sum := &+ res;
    average := sum/no;
    r := <range, average, res>;
    print r;
    return r;
end function;

ranges := [8,16,32,64,128,256];

RuntimesCftlByRanges := function( G, ranges, no )
    results := [];
    for range in ranges do
        r := RuntimesCftl( G, range, no );
        Append( ~results, r );
    end for;
    return results;
end function;    

// - save important examples as text files
// - write test function in magma
