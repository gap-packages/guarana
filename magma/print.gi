#############################################################################
#W print.gi                GUARANA package                     Bjoern Assmann
##
## Functions to print polycyclic presentations to files  
##
#H  @(#)$Id$
##
#Y 2006
##

GUARANA.Print2FilePcpMagmaStyle := function( G, grpName, file )
    local gens, rods, n, conj, g, h;
    
    gens := Pcp( G );
    rods := RelativeOrdersOfPcp( gens );
    n := Length( gens );

    PrintTo( file, grpName, " := PolycyclicGroup<" );
    for g in [1..n] do 
	if g <> n then
	    AppendTo( file, " ", gens[g], "," );
	else 
	    AppendTo( file,  " ", gens[g] );
	fi;
    od;
    AppendTo( file,  " | \n\n" );
    for g in [1..n] do
        if rods[g] <> 0 then
            ##  print the power relation for g.
	    if gens[g]^rods[g] = gens[g]^0 then 
                AppendTo( file,  "    ", gens[g], "^", rods[g], ",\n" );
	    else
                AppendTo( file,  "    ", gens[g], "^", rods[g], " = ", 
                      gens[g]^rods[g], ",\n" );
	    fi;
        fi;
    od;
    if rods <> 0 * rods then AppendTo( file,  "\n" ); fi;

    for h in [1..n] do
        for g in [1..h-1] do
            conj := gens[h]^gens[g];
            ##  print the conjuagte relation for h^g.
	    if rods[g] <> 0 then 
		if h = n and g=n-1 then
                    AppendTo( file,  "    ", gens[h], "^", gens[g], " = ", 
                              gens[h]^gens[g], "\n" );
	        else 
                    AppendTo( file,  "    ", gens[h], "^", gens[g], " = ", 
                              gens[h]^gens[g], ",\n" );
	        fi;
	    else  
		if h=n and g=n-1 then 
                    AppendTo( file,  "    ", gens[h], "^", gens[g], " = ", 
                              gens[h]^gens[g], ",\n" );
                    AppendTo( file,  "    ", gens[h], "^", gens[g]^-1, " = ",
		              gens[h]^(gens[g]^-1), "\n" );
                else
                    AppendTo( file,  "    ", gens[h], "^", gens[g], " = ", 
                              gens[h]^gens[g], ",\n" );
		    AppendTo( file,  "    ", gens[h], "^", gens[g]^-1, " = ",
		              gens[h]^(gens[g]^-1), ",\n" );
	        fi;
            fi;
        od;
    od;
    AppendTo( file,  ">;\n" );
end;

# - save important examples as text files
# - write test function in magma
