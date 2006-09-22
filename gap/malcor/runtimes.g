#############################################################################
##
#W runtimes.g             GUARANA package                     Bjoern Assmann
##
#H  @(#)$Id$
##
#Y 2006
##
##

## Aim: 
## - given a group compute the runtimes for 
##   bch setup, log and exp pols.
## 
## - given a group string like "Tr_n_O1" and a range for n
##   compute the above for all groups (together with HirschLength and 
##   nilpotency class.
## 
## - a function that converts this into latex code

GUARANA.RuntimeMalcevBCHSetUp := function( N )
    return GUARANA.CompleteRuntime2( MalcevObjectByTGroup, N ).time;
end;

GUARANA.RuntimeLogAndExpPols := function( malObj )
    return GUARANA.CompleteRuntime1( AddLogAndExpPolynomials, malObj );
end;
#Example
if false then 
    mo_F_28 := GUARANA.Get_FNG_MalcevObject( 2, 8 );
    GUARANA.RuntimeLogAndExpPols( mo_F_29 );
fi;

GUARANA.RuntimeBCHSetupAndLogAndExpPols := function( N )
    local res, time_setup, malObj, time_pols;

    res := GUARANA.CompleteRuntime2( MalcevObjectByTGroup, N );
    time_setup := res.time;
    malObj := res.result;
    time_pols := GUARANA.CompleteRuntime1( AddLogAndExpPolynomials, malObj );
    return [ time_setup, time_pols ];
end;


## IN
## range_n ............................. which n (degree of matrix groups
##                                       should be covered).
##
GUARANA.AverageRuntimesSetup_Tr_1_n_O1 := function( range_n )
    local results, res, N, info, r, n;

    results := [];
    Add( results, [ "Tr_1_n_O1" ] );
    for n in range_n do 
        res := [];
        N := GUARANA.Tr_n_O1( n )[4];
        info := [ "dim = ", n, "class = ", n-1, "hl = ", HirschLength(N) ] ;
        Add( res,info );
        r := GUARANA.RuntimeBCHSetupAndLogAndExpPols( N );
        Add( res, r );
        Add( results, res );
    od;
    return results;
end;

GUARANA.AverageRuntimesSetup_Tr_1_n_O2 := function( range_n )
    local results, res, N, info, r, n;

    results := [];
    Add( results, [ "Tr_1_n_O2" ] );
    for n in range_n do 
        res := [];
        N := GUARANA.Tr_n_O2( n )[4];
        info := [ "dim = ", n, "class = ", n-1, "hl = ", HirschLength(N) ] ;
        Add( res,info );
        r := GUARANA.RuntimeBCHSetupAndLogAndExpPols( N );
        Add( res, r );
        Add( results, res );
    od;
    return results;
end;

GUARANA.AverageRuntimesSetup_F_nc := function( n, range_c )
    local results, res, N, info, r,c;

    results := [];
    Add( results, [ "F_nc" ] );
    for c in range_c do 
        res := [];
        N := GUARANA.Examples_FreeNilpotentGrp( n, c );
        info := [ "n = ", n, "class = ", c, "hl = ", HirschLength(N) ] ;
        Add( res,info );
        r := GUARANA.RuntimeBCHSetupAndLogAndExpPols( N );
        Add( res, r );
        Add( results, res );
    od;
    return results;
end;


GUARANA.Latex_RuntimesSetupToLatex := function( list )
    local groupClass, res, k, ll, n, c, stringLine, le, i, j,hl;

    groupClass := list[1][1];
    res := "";
    k := Length( list );
    for i in [2..k] do 
        ll := list[i];
        n := ll[1][2];
        c := ll[1][4];
        hl := ll[1][6];
        stringLine := GUARANA.Latex_GroupClass2Latex( groupClass, n, c );
        Append( stringLine, " & " );
        Append( stringLine, String( hl ) );
        Append( stringLine, " & " );
        Append( stringLine, String( c ) );
        Append( stringLine, " & " );

        le := Length( ll[2] );
        for j in [1..le] do 
            Append( stringLine, " & " );
            Append( stringLine, String( ll[2][j] ) );
        od;
        Append( stringLine, " \\\\ \\hline \n" );
        Append( res, stringLine );
    od;
    return res;

end;

# Example
if false then 
    ranges_n := [2..5];
    res := GUARANA.AverageRuntimesSetup_Tr_1_n_O1( ranges_n );
    latex_code := GUARANA.Latex_RuntimesSetupToLatex( res );
    Print( latex_code );
fi;

GUARANA.Latex_GenerateRuntimesSetupTable := function( class_string, range_n, range_c )
    local res, latex_code;

    if class_string = "Tr_1_n_O1" then  
        res := GUARANA.AverageRuntimesSetup_Tr_1_n_O1( range_n );
    elif class_string = "Tr_1_n_O2" then  
        res := GUARANA.AverageRuntimesSetup_Tr_1_n_O2( range_n );
    elif class_string = "F_nc" then
        res := GUARANA.AverageRuntimesSetup_F_nc( range_n[1], range_c );
    fi;
    latex_code := GUARANA.Latex_RuntimesSetupToLatex( res );
    Print( latex_code );
    return latex_code;
end;

# Example
if false then 
    range_n := [2..8];
    range_n := [7];
    range_c := [0];
    class_string := "Tr_1_n_O1";
    latex_code1 := GUARANA.Latex_GenerateRuntimesSetupTable( class_string, range_n, range_c );
    Print( latex_code1 );

    range_n := [2..7];
    range_c := [0];
    class_string := "Tr_1_n_O2";
    latex_code2 := GUARANA.Latex_GenerateRuntimesSetupTable( class_string, range_n, range_c );
    Print( latex_code2 );
    class_string := "Tr_1_n_O2";
    latex_code := GUARANA.Latex_GenerateRuntimesSetupTable( class_string, range_n, range_c );
    Print( latex_code );

    range_n := [3];
    range_c := [2..6];
    class_string := "F_nc";
    latex_code := GUARANA.Latex_GenerateRuntimesSetupTable( class_string, range_n, range_c );
    Print( latex_code );
fi;

#############################################################################
##
## Runtimes for Log
##
GUARANA.AverageRuntimeLog := function(  malcevObject, noTests,range, method )
    local results, g, r, i,aver;

	SetLogAndExpMethod( malcevObject, method );
    results := [];
    for i in [1..noTests] do
        g := RandomGrpElm( malcevObject, range );
        r := GUARANA.CompleteRuntime1( Log, g ); 
        Add( results, r );
    od;
    aver := Int( Sum( results )/noTests );
    return [aver, results];
end;

if false then
    # get malcev setup for F_nc
    method := "pols";
    #method := "simple";
    range := 1024;
    noTests := 100;
    n := 3;
    c := 6;
    mo_Fnc := GUARANA.Get_FNG_MalcevObject( n, c );

    res := GUARANA.AverageRuntimeLog( mo_Fnc, noTests, range, method );
fi;


#############################################################################
##
#E
