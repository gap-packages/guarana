#############################################################################
#
# Code for testing the fast multiplication algorithms and for the
# determination of runtimes.

SC_RandomPcpElement := function( G, range )
    local pcp,n,ll,vec,g;
    pcp := Pcp( G );
    n := Length( pcp );
    ll := [-range..range];
    vec := List( [1..n], x-> RandomList( ll ) );
    g := MappedVector( vec, pcp );
    return g;
end;


SC_StandardDeviationAndMean := function( list )
    local n,sum,mean,vec,diff,x,deviation,variance;
    n := Length( list );
    sum := Sum( list );
    mean := sum/n;

    #compute deviation
    vec := List( [1..n], x -> mean );
    diff := vec - list;
    x := ScalarProduct( diff, diff );
    variance := x/n;
    deviation := RootInt( Int( variance ) );

    return rec( mean := mean, deviation := deviation, 
                variance := variance ); 
    
end;

SC_ComputeAverageRuntimeFastMultiplication := function( FCR, range, no )
    local results,i,g1,g2,time,sum,average,collectFunc,GG;

    # decide which function to use
    if FCR.info = "bch" then 
        collectFunc := BCH_FastMultiplicationAbelianSemiTGroup_RuntimeVersion;
        GG := FCR.recNewParent.GG;
    elif FCR.info = "mat" then 
        collectFunc := MAT_FastMultiplicationAbelianSemiTGroup_RuntimeVersion;
        GG := FCR.GG;
    else 
        return fail;
    fi;

    results := [];
    # do 'no' random multiplications 
    for i in [1..no] do
        #Print( "Multiplication no ", i, "\n" );
        g1 := SC_RandomPcpElement( GG, range );
        #Print( "Used g1 : ", g1, "\n" );
        g2 := SC_RandomPcpElement( GG, range );
        #Print( "Used g2 : ", g2, "\n" );
        time := POL_CompleteRuntime2( collectFunc, [FCR,g1,g2]  ).time;
        Add( results, time );
        #Print( "Result : ", res.result , "\n" );
    od;

    #compute average
    sum := Sum( results );
    average := StringTime( Int( sum/no ) );
       
    #HR = human readable
    return rec( results := results, average := Int( sum/no), 
                average_hr := average );

end;

SC_CollectionFromTheLeft_RuntimeVersion := function( args )
    local g1,g2,g;
    g1 := args[1];
    g2 := args[2];

    g := Exponents( g1*g2 );
  
    return g;
end;

##Mark#
#
#T := SC_Exams( 7 );
#FCR := MAT_FastConjugationRec( T.G, T.N );;
#SC_CompareAverageRuntimeCollectionFromTheLeftVersusMalcev( FCR, 2, 50 );

#############################################################################
#
# compare malcev with collection from the left
SC_CompareAverageRuntimeCFTLVersusMalcev := function( FCR, range, no )
    local results,results2,i,g1,g2,sum,sum2,deviation,deviation2,test,res,res2,
          av1,av2,mean,mean2,collectFunc,GG;

    # decide which function to use
    if FCR.info = "bch" then 
        collectFunc := BCH_FastMultiplicationAbelianSemiTGroup_RuntimeVersion;
        GG := FCR.recNewParent.GG;
    elif FCR.info = "mat" then 
        collectFunc := MAT_FastMultiplicationAbelianSemiTGroup_RuntimeVersion;
        GG := FCR.GG;
    else 
        return fail;
    fi;

    results := [];
    results2 := [];
    
    # do "no" random multiplications 
    for i in [1..no] do
        Print( "Multiplication no ", i, "\n" );
        g1 := SC_RandomPcpElement( GG, range );
        Print( "Used g1 : ", g1, "\n" );
        g2 := SC_RandomPcpElement( GG, range );
        Print( "Used g2 : ", g2, "\n" );
        
        # Malcev computation
        Print( "Malcev: \n" );
        res2 := POL_CompleteRuntime2( 
                      collectFunc,
                      [FCR,g1,g2]  );
        Add( results2, res2.time );
        Print( "Result : ", res2.result , "\n" );
       
        # collection from the left
        Print( "Collection from the left \n" );
        res := POL_CompleteRuntime2( 
                    SC_CollectionFromTheLeft_RuntimeVersion,
                    [g1,g2]  );
        Add( results, res.time );

        # test if exponent vector is the same
        if res.result<> res2.result then
            Error( "Mist\n" );
        fi;
        Print( "\n" );
    od;

    #compute average and mean
    av1 := SC_StandardDeviationAndMean( results );
    mean := StringTime( Int( av1.mean ) );
    deviation := StringTime( Int( av1.deviation ) );
    av2 := SC_StandardDeviationAndMean( results2 );
    mean2 := StringTime( Int( av2.mean ) );
    deviation2 := StringTime( Int( av2.deviation ) );

    Print( "\n" );
    return rec( results := results, mean := mean, deviation := deviation,
                variance := av1.variance,
                results_fast := results2, mean_fast := mean2,
                deviation_fast := deviation2,
                variance_fast := av2.variance  );
end;

#############################################################################
#
# computes the runtime as a function in the range.
# The range varies between base^0,base^1,base^2, ..., upperBound
# 
SC_MalcevRuntimeByRange := function( FCR, base, upperBound )
    local no, result, range,time,counter;
   
    no := 1;
    result := [];
    counter := 1;
    range := base;
    while range <= upperBound do 
       time := SC_ComputeAverageRuntimeFastMultiplication( 
                                         FCR,range , no ).average;
       Add( result, [range, time] );
       Print( "range: ", range, " time: ", time, "\n" );
       range := range*base;
    od;
    return result;
end; 


#############################################################################
#
#F prepare plot data
#
SC_List2PlotStream := function( ll, filename )
    local str, out, i,n;
    str := "";; 
    out := OutputTextString(str,true);;
 
    n := Length( ll );
    for i in [1..n] do
        AppendTo( out, i, " " );
        AppendTo( out, ll[i], "\n" );
    od;
    CloseStream( out );
    PrintTo( filename, str );

    return str;

end;
