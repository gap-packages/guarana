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
        print i;
        t := RuntimeCftl( G, range );
        print "Maximum memory usage ", GetMaximumMemoryUsage();
        print "Current memory usage ", GetMemoryUsage();
        print "time ", t;
        Append( ~res, t );
    end for;
    sum := &+ res;
    average := sum/no;
    r := <range, average, res>;
    print r;
    print "\n";
    return r;
end function;

RuntimesCftlByRanges := function( G, ranges, no )
    results := [];
    for range in ranges do
        r := RuntimesCftl( G, range, no );
        Append( ~results, r );
    end for;
    return results;
end function;    

/*
// set max memory in the cshell to 1GB
// and virtual memory to 256MB or 512MB
limit memoryuse 1048576
limit vmemoryuse 262144
limit vmemoryuse 524288
limit vmemoryuse 1048576

// Set memory limit to 1GB. 
SetMemoryLimit( 1024^3 );
load "../random.m";
load "Tr_4_O1.m";
load "Tr_5_O1.m";
load "Tr_6_O1.m";
load "Tr_7_O1.m";

ranges := [ 8, 16, 32, 64 ];
ranges := [100];
no := 10;
G := Tr_5_O1;
RuntimesCftlByRanges( G, ranges, no );

*/
