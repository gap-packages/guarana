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
        sum := &+ res;
        average := sum/i;
        print "Average so far : ", average;
    end for;
    sum := &+ res;
    average := sum/no;
    r := <res,range, average>;
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

load "Tr_2_O1.m";
load "Tr_3_O1.m";
load "Tr_4_O1.m";
load "Tr_5_O1.m";
load "Tr_6_O1.m";
load "Tr_7_O1.m";
load "Tr_8_O1.m";

load "Tr_2_O2.m";
load "Tr_3_O2.m";
load "Tr_4_O2.m";
load "Tr_5_O2.m";
load "Tr_6_O2.m";
load "Tr_7_O2.m";

load "F_22_Aut1.m";
load "F_23_Aut1.m";
load "F_24_Aut1.m";
load "F_25_Aut1.m";
load "F_26_Aut1.m";
load "F_27_Aut1.m";
load "F_28_Aut1.m";

load "F_32_Aut2.m";
load "F_33_Aut2.m";
load "F_34_Aut2.m";
load "F_35_Aut2.m";
load "F_36_Aut2.m";


SetMemoryLimit( 1024^3 );
load "../random.m";
load "F_22_Aut1.m"

ranges := [1000];
no := 1000;
G := F_22_Aut1;
times := RuntimesCftlByRanges( G, ranges, no );


SetMemoryLimit( 1024^3 );
load "../random.m";
load "Tr_2_O2.m";
load "Tr_3_O2.m";
load "Tr_4_O2.m";
ranges := [10];
no := 1000;
G := Tr_3_O2;
times := RuntimesCftlByRanges( G, ranges, no );


SetMemoryLimit( 1024^3 );
load "../random.m";
load "Tr_6_O1.m";
ranges := [1..10];
G := Tr_6_O1;
no := 1000;
times := RuntimesCftlByRanges( G, ranges, no );


*/

PrintRuntimes := function( l )
    for a in l do 
        print a[2], a[3];
    end for;
    return 1;
end function;
        
