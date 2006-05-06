#############################################################################
#W print.gi                GUARANA package                     Bjoern Assmann
##
## Functions to print polycyclic presentations to files  
##
#H  @(#)$Id$
##
#Y 2006
##

GUARANA.PrintPcpMagmaStyle := function( pcp )
    local collector, gens, rods, n, conj, g, h;

    collector := Collector( G );
    gens := Pcp( G );
    rods := RelativeOrdersOfPcp( gens );
    n := Length( gens );

    Print( "<" );
    for g in [1..n] do Print( " ", gens[g] ); od;
    Print( " | \n\n" );
    for g in [1..n] do
        if rods[g] <> 0 then
            ##  print the power relation for g.
            Print( "    ", gens[g], "^", rods[g], " = ", 
                   gens[g]^rods[g], "\n" );
        fi;
    od;
    if rods <> 0 * rods then Print( "\n" ); fi;

    for h in [1..n] do
        for g in [1..h-1] do
            conj := gens[h]^gens[g];
            if conj <> gens[h] then
                ##  print the conjuagte relation for h^g.
                Print( "    ", gens[h], "^", gens[g], " = ", 
                       gens[h]^gens[g], "\n" );
		Print( "    ", gens[h], "^", gens[g]^-1, " = ",
		       gens[h]^gens[g]^-1, "\n" );
            fi;
        od;
    od;
    Print( ">\n" );
end;
