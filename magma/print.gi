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

GUARANA.PrintExamClassesMagmaStyle := function( exam_string, ranges_n, 
                                                ranges_c )
    local G, s1, s2, s3, grp_name, file_name, s4, n, c;
    if exam_string = "Tr_n_O1" then 
        for n in ranges_n do 
            # get group
            G := GUARANA.Tr_n_O1( n )[1];
            # get name of group and magma file
            s1 := "Tr_";
            s2 := String( n );
            s3 := "_O1";
            grp_name := Concatenation( s1, s2, s3 );
            file_name := Concatenation( grp_name, ".m" );
            GUARANA.Print2FilePcpMagmaStyle( G, grp_name, file_name );
        od;
    elif exam_string = "Tr_n_O2" then 
        for n in ranges_n do 
            # get group
            G := GUARANA.Tr_n_O2( n )[1];
            # get name of group and magma file
            s1 := "Tr_";
            s2 := String( n );
            s3 := "_O2";
            grp_name := Concatenation( s1, s2, s3 );
            file_name := Concatenation( grp_name, ".m" );
            GUARANA.Print2FilePcpMagmaStyle( G, grp_name, file_name );
        od;
    elif exam_string = "F_nc_Aut1" then 
        for n in ranges_n do 
            for c in ranges_c do 
                # get group
                G := GUARANA.F_nc_Aut1( n, c )[1];
                # get name of group and magma file
                s1 := "F_";
                s2 := String( n );
                s3 := String( c );
                s4 := "_Aut1";
                grp_name := Concatenation( s1, s2, s3, s4 );
                file_name := Concatenation( grp_name, ".m" );
                GUARANA.Print2FilePcpMagmaStyle( G, grp_name, file_name );
            od;
        od;
    elif exam_string = "F_nc_Aut2" then 
        for n in ranges_n do 
            for c in ranges_c do 
                # get group
                G := GUARANA.F_nc_Aut2( n, c )[1];
                # get name of group and magma file
                s1 := "F_";
                s2 := String( n );
                s3 := String( c );
                s4 := "_Aut2";
                grp_name := Concatenation( s1, s2, s3, s4 );
                file_name := Concatenation( grp_name, ".m" );
                GUARANA.Print2FilePcpMagmaStyle( G, grp_name, file_name );
            od;
        od;
    fi;
end;

if false then 
    exam_string := "Tr_n_O2";
    ranges_n := [2..3];
    ranges_c := [0];
    GUARANA.PrintExamClassesMagmaStyle( exam_string, ranges_n, ranges_c );

    exam_string := "F_nc_Aut1";
    ranges_n := [2];
    ranges_c := [2..8];
    GUARANA.PrintExamClassesMagmaStyle( exam_string, ranges_n, ranges_c );

    exam_string := "F_nc_Aut2";
    ranges_n := [3];
    ranges_c := [2..6];
    GUARANA.PrintExamClassesMagmaStyle( exam_string, ranges_n, ranges_c );
fi;


# - save important examples as text files
# - write test function in magma
