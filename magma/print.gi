#############################################################################
#W print.gi                GUARANA package                     Bjoern Assmann
##
## Functions to print polycyclic presentations to files  
##
#H  @(#)$Id$
##
#Y 2006
##

GUARANA.PrintPcpMagmaStyle := function( G, grpName )
    local gens, rods, n, conj, g, h;

    gens := Pcp( G );
    rods := RelativeOrdersOfPcp( gens );
    n := Length( gens );

    Print( grpName, " := PolycyclicGroup<" );
    for g in [1..n] do 
	if g <> n then
	    Print( " ", gens[g], "," );
	else 
	    Print( " ", gens[g] );
	fi;
    od;
    Print( " | \n\n" );
    for g in [1..n] do
        if rods[g] <> 0 then
            ##  print the power relation for g.
	    if gens[g]^rods[g] = gens[g]^0 then 
                Print( "    ", gens[g], "^", rods[g], ",\n" );
	    else
                Print( "    ", gens[g], "^", rods[g], " = ", 
                      gens[g]^rods[g], ",\n" );
	    fi;
        fi;
    od;
    if rods <> 0 * rods then Print( "\n" ); fi;

    for h in [1..n] do
        for g in [1..h-1] do
            conj := gens[h]^gens[g];
            if h = n and g = h-1 then 
                ##  print the conjuagte relation for h^g.
                Print( "    ", gens[h], "^", gens[g], " = ", 
                       gens[h]^gens[g], ",\n" );
		Print( "    ", gens[h], "^", gens[g]^-1, " = ",
		       gens[h]^(gens[g]^-1), "\n" );
            else
                ##  print the conjuagte relation for h^g.
                Print( "    ", gens[h], "^", gens[g], " = ", 
                       gens[h]^gens[g], ",\n" );
		Print( "    ", gens[h], "^", gens[g]^-1, " = ",
		       gens[h]^(gens[g]^-1), ",\n" );
            fi;
        od;
    od;
    Print( ">;\n" );
end;

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
            if h = n and g = h-1 then 
                ##  print the conjuagte relation for h^g.
                AppendTo( file,  "    ", gens[h], "^", gens[g], " = ", 
                       gens[h]^gens[g], ",\n" );
		AppendTo( file,  "    ", gens[h], "^", gens[g]^-1, " = ",
		       gens[h]^(gens[g]^-1), "\n" );
            else
                ##  print the conjuagte relation for h^g.
                AppendTo( file,  "    ", gens[h], "^", gens[g], " = ", 
                       gens[h]^gens[g], ",\n" );
		AppendTo( file,  "    ", gens[h], "^", gens[g]^-1, " = ",
		       gens[h]^(gens[g]^-1), ",\n" );
            fi;
        od;
    od;
    AppendTo( file,  ">;\n" );
end;

# - write function for writing in files
# - save important examples as text files
# - write random by range function in magma
# - write test function in magam
