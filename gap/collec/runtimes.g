#############################################################################
##
#W runtimes.g             GUARANA package                     Bjoern Assmann
##
#H  @(#)$Id$
##
#Y 2006
##
##
GUARANA.ProdRuntimeVers := function( args )
    local g,h;

    g := args[1];
    h := args[2];
    return g*h;
end;

GUARANA.AverageRuntimeCollec := function( malCol, ranges, no )
    local results, times, g, h, time, sum, average, range, i;

    results := [];
    for range in ranges do 
        times := [];
        for i in [1..no] do
    
            g := Random( malCol, range );
            h := Random( malCol, range );
        
            # compute time for product
            time := GUARANA.CompleteRuntime2( GUARANA.ProdRuntimeVers, 
                                          [g,h] ).time;

            Add( times, time );
        od;

        # compute average
        sum := Sum( times );
        Print( "times : " , times, "\n" );
        average :=  Int( sum/no ) ;
        Print( "average : ", average, "\n" );
        Add( results, [range, average] );
    od;
    return results;
end;

# a function that gives the test results for a class of examples
# like Tr_1(n,Q) etc.

## IN
## range_n ............................. which n (degree of matrix groups
##                                       should be covered).
##
GUARANA.AverageRuntimesCollec_Tr_n_O1 := function( range_n, ranges, no )
    local results, res, malCol, info, r, n;

    results := [];
    Add( results, [ "Tr_n_O1" ] );
    for n in range_n do 
        res := [];
        malCol := GUARANA.MalcevColl_Tr_n_O1( n );
        info := [ "dim = ", n, "class = ", n-1 ] ;
        Add( res,info );
        r := GUARANA.AverageRuntimeCollec( malCol, ranges, no );
        Add( res, r );
        Print( r, " " );
        Add( results, res );
    od;
    return results;
end;

GUARANA.AverageRuntimesCollec_Tr_n_O2 := function( range_n, ranges, no )
    local results, res, malCol, info, r, n;

    results := [];
    Add( results, [ "Tr_n_O2" ] );
    for n in range_n do 
        res := [];
        malCol := GUARANA.MalcevColl_Tr_n_O2( n );
        info := [ "dim = ", n, "class = ", n-1 ] ;
        Add( res,info );
        r := GUARANA.AverageRuntimeCollec( malCol, ranges, no );
        Add( res, r );
        Print( r, " " );
        Add( results, res );
    od;
    return results;
end;

GUARANA.AverageRuntimesCollec_F_nc_Aut1 := function( range_n, range_c, 
                                                     ranges, no )
    local results, res, malCol, info, r, c, n;

    results := [];
    Add( results, [ "F_nc_Aut1" ] );
    for c in range_c do 
        for n in range_n do 
            res := [];
            malCol := GUARANA.MalcevColl_F_nc_Aut1( n, c );
            info := [ "dim = ", n, "class = ", c ] ;
            Add( res,info );
            r := GUARANA.AverageRuntimeCollec( malCol, ranges, no );
            Add( res, r );
            Print( r, " " );
            Add( results, res );
        od;
    od;
    return results;
end;

GUARANA.AverageRuntimesCollec_F_nc_Aut2 := function( range_n, range_c, 
                                                     ranges, no )
    local results, res, malCol, info, r, c, n;

    results := [];
    Add( results, [ "F_nc_Aut2" ] );
    for c in range_c do 
        for n in range_n do 
            res := [];
            malCol := GUARANA.MalcevColl_F_nc_Aut2( n, c );
            info := [ "dim = ", n, "class = ", c ] ;
            Add( res,info );
            r := GUARANA.AverageRuntimeCollec( malCol, ranges, no );
            Add( res, r );
            Print( r, " " );
            Add( results, res );
        od;
    od;
    return results;
end;

GUARANA.Latex_GroupClass2Latex := function( str, n, c )
    local ss, s1, s2, s3, s4,s2b;

    if str = "Tr_n_O1" then
        s1 := "$G(\\mathrm{Tr}_";
        s2 := String( n );
        s3 := "(\\mathcal{O}_1))$";
        ss := Concatenation( s1, s2, s3 );
    elif str = "Tr_n_O2" then 
        s1 := "$G(\\mathrm{Tr}_";
        s2 := String( n );
        s3 := "(\\mathcal{O}_2))$";
        ss := Concatenation( s1, s2, s3 );
    elif str = "F_nc_Aut1" then 
        s1 := "$G( \\langle \\bar{\\varphi_1} \\rangle \\rtimes F_{";
        s2 := String( n );
        s3 := String( c );
        s4 := "})$";
        ss := Concatenation( s1, s2, s3, s4 );
    elif str = "F_nc_Aut2" then 
        s1 := "$G( \\langle \\bar{\\varphi_2} \\rangle \\rtimes F_{";
        s2 := String( n );
        s3 := String( c );
        s4 := "})$";
        ss := Concatenation( s1, s2, s3, s4 );
    elif str = "Tr_1_n_O1" then 
        s1 := "$G(\\mathrm{Tr}_1(";
        s2 := String( n );
        s3 := ", \\mathcal{O}_1))$";
        ss := Concatenation( s1, s2, s3 );
    elif str = "Tr_1_n_O2" then 
        s1 := "$G(\\mathrm{Tr}_1(";
        s2 := String( n );
        s3 := ", \\mathcal{O}_2))$";
        ss := Concatenation( s1, s2, s3 );
    elif str = "F_nc" then
        s1 := "$F_{";
        s2 := String( n );
        s2b := ",";
        s3 := String( c );
        s4 := "}$";
        ss := Concatenation( s1, s2, s2b, s3, s4 );
    fi;
    return ss;
end;

GUARANA.Latex_RuntimesToLatex := function( list )
    local groupClass, res, k, ll, n, c, stringLine, le, i, j;

    groupClass := list[1][1];
    res := "";
    k := Length( list );
    for i in [2..k] do 
        ll := list[i];
        n := ll[1][2];
        c := ll[1][4];
        stringLine := GUARANA.Latex_GroupClass2Latex( groupClass, n, c );

        le := Length( ll[2] );
        for j in [1..le] do 
            Append( stringLine, " & & " );
            Append( stringLine, String( ll[2][j][2] ) );
        od;
        Append( stringLine, " \\\\ \\hline \n" );
        Append( res, stringLine );
    od;
    return res;

end;

# Example
if false then 
    ranges := [2, 4, 8 ];
    ranges_n := [3..5];
    no := 100;
    res := GUARANA.AverageRuntimesCollec_Tr_n_O1( ranges_n, ranges, no );
    latex_code := GUARANA.Latex_RuntimesToLatex( res );
    Print( latex_code );

    range_n := [2];
    range_c := [4,5,6,7];
    no := 100;
    ranges := [2,4,8];
    res := GUARANA.AverageRuntimesCollec_F_nc_Aut1( range_n, range_c, ranges,
                                                    no );
    latex_code := GUARANA.Latex_RuntimesToLatex( res );
    Print( latex_code );
fi;

GUARANA.Latex_GenerateRuntimesTable := function( class_string, ranges, range_n, range_c, no )
    local res, latex_code;

    if class_string = "Tr_n_O1" then  
        res := GUARANA.AverageRuntimesCollec_Tr_n_O1( range_n, ranges, no );
    elif class_string = "Tr_n_O2" then  
        res := GUARANA.AverageRuntimesCollec_Tr_n_O2( range_n, ranges, no );
    elif class_string = "F_nc_Aut1" then
        res:=GUARANA.AverageRuntimesCollec_F_nc_Aut1(range_n,range_c,ranges,no);
    elif class_string = "F_nc_Aut2" then
        res:=GUARANA.AverageRuntimesCollec_F_nc_Aut2(range_n,range_c,ranges,no);
    fi;
    latex_code := GUARANA.Latex_RuntimesToLatex( res );
    Print( latex_code );
    return latex_code;
end;

# Example
if false then 
    class_string := "Tr_n_O1";
    ranges := [1,10,100,1000];
    range_n := [4..8];
    range_c := [0];
    no := 100;
    GUARANA.Latex_GenerateRuntimesTable( class_string, ranges, range_n,
    range_c, no );

    class_string := "Tr_n_O2";
    ranges := [1,10,100,1000];
    range_n := [2..7];
    range_c := [0];
    no := 100;
    GUARANA.Latex_GenerateRuntimesTable( class_string, ranges, range_n,
    range_c, no );

    class_string := "F_nc_Aut1";
    ranges := [1,10,100,1000];
    range_n := [2];
    range_c := [2..8];
    no := 100;
    GUARANA.Latex_GenerateRuntimesTable( class_string, ranges, range_n,
    range_c, no );
    
    class_string := "F_nc_Aut2";
    ranges := [1,10,100,1000];
    range_n := [3];
    range_c := [6];
    no := 100;
    latex_times := GUARANA.Latex_GenerateRuntimesTable( class_string, ranges, 
                                            range_n, range_c, no );

    class_string := "F_nc_Aut2";
    ranges := [1000];
    range_n := [3];
    range_c := [6];
    no := 1000;
    latex_times := GUARANA.Latex_GenerateRuntimesTable( class_string, ranges, 
                                            range_n, range_c, no );

GUARANA.Latex_GenerateRuntimesTable( "Tr_n_O1", [2,4,8], [2..4], [1], 100 );
GUARANA.Latex_GenerateRuntimesTable( "Tr_n_O2", [2,4,8], [2..4], [1], 100 );
GUARANA.Latex_GenerateRuntimesTable( "F_nc_Aut1", [2,4,8], [2..3], [2..4], 100 );
GUARANA.Latex_GenerateRuntimesTable( "F_nc_Aut2", [2,4,8], [2..3], [2..4], 100 );


fi;

PrintList := function( list )
    local a;
    for a in list do 
        Print( a[1], " ", a[2], "\n" );
    od;
end;

#############################################################################
##
## Methods to get the data for the graph in the BCH paper.
##
## We consider the group Tr_6_O1
##
##
if false then 
ranges := [1..200];
no := 100;

# compute the times
times := GUARANA.AverageRuntimesCollec_Tr_n_O1( [6], ranges, no );

# print in format range, time, \n
# to a file
LogTo( "times_Tr_6_O1.txt" );
PrintList( times[2][2] );
LogTo();


fi;

if false then
ranges := [1..10];
no := 100;
times := GUARANA.AverageRuntimesCollec_Tr_n_O1( [6], ranges, no );
LogTo( "times_smallRange_Tr_6_O1.txt" );
PrintList( times[2][2] );
LogTo();
fi;


#############################################################################
##
## Functions for computing the time needed for the complete setup
## of the Malcev collector
##

## Example
## GUARANA.RuntimesSetupCollector_Tr_n_O1( [2..4] );
##
GUARANA.RuntimesSetupCollector_Tr_n_O1 := function( range_n )
    local results, res, info, args, time, n, hl;
    
    results := [];
    Add( results, [ "Tr_n_O1" ] );
    for n in range_n do 
        # get pcs and subgroups C,N etc
        args := GUARANA.Tr_n_O1( n );

        # collect some meta data
        res := [];
        hl := HirschLength( args[1] );
        info := [ "dim = ", n, "class = ", n-1, "hl = ", hl  ] ;
        Add( res,info );

        # compute time needed for the setup
        time := GUARANA.CompleteRuntime2( MalcevCollectorConstruction, 
                                          args ).time;
        Add( res, time );
        Add( results, res );
    od;
    return results;
end;

GUARANA.RuntimesSetupCollector_Tr_n_O2 := function( range_n )
    local results, res, info, args, time, n, hl;
    
    results := [];
    Add( results, [ "Tr_n_O2" ] );
    for n in range_n do 
        # get pcs and subgroups C,N etc
        args := GUARANA.Tr_n_O2( n );

        # collect some meta data
        res := [];
        hl := HirschLength( args[1] );
        info := [ "dim = ", n, "class = ", n-1, "hl = ", hl  ] ;
        Add( res,info );

        # compute time needed for the setup
        time := GUARANA.CompleteRuntime2( MalcevCollectorConstruction, 
                                          args ).time;
        Add( res, time );
        Add( results, res );
    od;
    return results;
end;

GUARANA.RuntimesSetupCollector_F_nc_Aut1 := function( range_n, range_c )
    local results, args, res, hl, info, time, c, n;
    
    results := [];
    Add( results, [ "F_nc_Aut1" ] );
    for c in range_c do 
        for n in range_n do 
            # get pcs and subgroups C,N etc
            args := GUARANA.F_nc_Aut1( n,c );

            # collect some meta data
            res := [];
            hl := HirschLength( args[1] );
            info := [ "dim = ", n, "class = ", c, "hl = ", hl  ] ;
            Add( res,info );

            # compute time needed for the setup
            time := GUARANA.CompleteRuntime2( MalcevCollectorConstruction, 
                                          args ).time;
            Add( res, time );
            Add( results, res );
        od;
    od;
    return results;
end;


GUARANA.RuntimesSetupCollector_F_nc_Aut2 := function( range_n, range_c )
    local results, args, res, hl, info, time, c, n;
    
    results := [];
    Add( results, [ "F_nc_Aut2" ] );
    for c in range_c do 
        for n in range_n do 
            # get pcs and subgroups C,N etc
            args := GUARANA.F_nc_Aut2( n,c );

            # collect some meta data
            res := [];
            hl := HirschLength( args[1] );
            info := [ "dim = ", n, "class = ", c, "hl = ", hl  ] ;
            Add( res,info );

            # compute time needed for the setup
            time := GUARANA.CompleteRuntime2( MalcevCollectorConstruction, 
                                          args ).time;
            Add( res, time );
            Add( results, res );
        od;
    od;
    return results;
end;

GUARANA.Latex_RuntimesSetupCollectorToLatex := function( list )
    local groupClass, res, k, ll, n, c, stringLine, hl, i;

    groupClass := list[1][1];
    res := "";
    k := Length( list );
    for i in [2..k] do 
        ll := list[i];

        # add group name
        n := ll[1][2];
        c := ll[1][4];
        stringLine := GUARANA.Latex_GroupClass2Latex( groupClass, n, c );

        # add hirsch length
        hl := ll[1][6];
        Append( stringLine, " & " );
        Append( stringLine, String( hl ) );

        # add runtime
        Append( stringLine, " & " );
        Append( stringLine, String( ll[2] ) );

        # add new line
        Append( stringLine, " \\\\ \\hline \n" );

        Append( res, stringLine );
    od;
    return res;

end;

GUARANA.Latex_GenerateCompleteSetupTable 
                := function( class_string, range_n, range_c )
    local res, latex_code;

    if class_string = "Tr_n_O1" then  
        res := GUARANA.RuntimesSetupCollector_Tr_n_O1( range_n );
    elif class_string = "Tr_n_O2" then  
        res := GUARANA.RuntimesSetupCollector_Tr_n_O2( range_n );
    elif class_string = "F_nc_Aut1" then
        res := GUARANA.RuntimesSetupCollector_F_nc_Aut1( range_n,range_c );
    elif class_string = "F_nc_Aut2" then
        res := GUARANA.RuntimesSetupCollector_F_nc_Aut2( range_n,range_c );
    fi;
    latex_code := GUARANA.Latex_RuntimesSetupCollectorToLatex( res );
    Print( latex_code );
    return latex_code;
end;

# Example
if false then 
    class_string := "Tr_n_O1";
    range_n := [2..8];
    range_c := [0];
    latex_code := GUARANA.Latex_GenerateCompleteSetupTable( class_string, range_n, range_c );

    class_string := "Tr_n_O2";
    range_n := [2..7];
    range_c := [0];
    latex_code := GUARANA.Latex_GenerateCompleteSetupTable( class_string, range_n, range_c );

    class_string := "F_nc_Aut1";
    range_n := [2];
    range_c := [2..8];
    latex_code := GUARANA.Latex_GenerateCompleteSetupTable( class_string, range_n, range_c );

    class_string := "F_nc_Aut1";
    range_n := [2];
    range_c := [8];
    latex_code := GUARANA.Latex_GenerateCompleteSetupTable( class_string, range_n, range_c );

    class_string := "F_nc_Aut2";
    range_n := [3];
    range_c := [2..6];
    latex_code := GUARANA.Latex_GenerateCompleteSetupTable( class_string, range_n, range_c );

    class_string := "F_nc_Aut2";
    range_n := [3];
    range_c := [6];
    latex_code := GUARANA.Latex_GenerateCompleteSetupTable( class_string, range_n, range_c );
fi;


#############################################################################
##
#E
