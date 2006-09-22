# functions for testing 
# in particular for testing the performance

#############################################################################
# Test functions:
#

# test free nilpotent groups 

# start determines the the position in the list ll from which
# it contains TGroup records
#
# Example
# 
# exams_unitr_2 := BCH_Get_Unitriangular_TGroupRecords( 10, 2 );
# BCH_ComputeTimeBCHSetup( exams_unitr_2, 2, 10, recBCH9 );
#
BCH_ComputeTimeBCHSetup := function( ll, start, stop, recBCH )
    local results,i,time1,time2,recLieAlg;
    results := [];

    # compute lie algebra via bch method for all examples
    for i in [start..stop] do
        time1 := POL_CompleteRuntime2(
                    BCH_SetUpLieAlgebraRecordByMalcevbasis,
                    ll[i]  ).time;
        recLieAlg := BCH_SetUpLieAlgebraRecordByMalcevbasis( ll[i] );
        time2 := POL_CompleteRuntime2(
                    BCH_ComputeStructureConstants,
                    [recBCH, recLieAlg]  ).time;
        Add( results, time1+time2 );
        recLieAlg := BCH_SetUpLieAlgebraRecordByMalcevbasis( ll[i] );
        BCH_ComputeStructureConstants( [recBCH, recLieAlg] );
        Print( TestJacobi( recLieAlg.scTable ) );
    od;
    return results;
end;

# results: On Thinkpad machine
# Free nilpotent on 2 gens from class 1 to 9
# [ 2, 4, 10, 20, 48, 118, 146, 250, 994 ]

# Free nilpotent group on 3 gens from class 1 to 6
# [ 0, 3, 7, 25, 145, 472 ]

# unitriangular group of dim 2 to 10 over maximal order of  x^2-3
#[ 37, 17, 63, 68, 209, 573, 2259, 7483, 26084 ]



# start determines the the position in the list ll from which
# it contains TGroup records
#
# rep_method .... method which is used to compute the matrix 
#                 representation. It is either
#                 deGraaf_Nickel or
#                 Nickel
#
# Example of use:
# exams_Fnc := BCH_Get_FNG_TGroupRecords( 3, 6 );
# MAT_ComputeTimeMatSetup( exams_Fnc, 2,9, "Nickel" );
#
# exams_unitr_2 := BCH_Get_Unitriangular_TGroupRecords( 7, 3 );
# MAT_ComputeTimeMatSetup(exams_unitr_2, 2, 5, "Nickel" );
#
MAT_ComputeTimeMatSetup := function( ll, start, stop, rep_method )
    local results,n,i,time;
    results := [];

    # compute lie algebra via mat method for all examples
    for i in [start..stop] do
        time := POL_CompleteRuntime2(
                        MAT_LieAlgebra,
                        [ll[i].N, rep_method ]  ).time;
        Add( results, time );
    od;
    return results;
end;

if false
    # for the BCH paper with Steve
    exams_Fnc := BCH_Get_FNG_TGroupRecords( 3, 6 );
    # compute times for F_32 to F_35
    MAT_ComputeTimeMatSetup( exams_Fnc, 3,6, "Nickel" );

    # to be exact substract the conjugator time ?
    # probably not because the representations is not integral anyway.


    # compute setup time for Tr_1_n_O2 for n = 2..6
    exams_unitr_2 := BCH_Get_Unitriangular_TGroupRecords( 7, 3 );
    MAT_ComputeTimeMatSetup(exams_unitr_2, 2, 6, "Nickel" );

    
fi;
    



BCH_AbstractLog_Simple_ByExponent_RuntimeVersion := function( args )
    local recLieAlg, recBCH, exp;
    recLieAlg := args[1];
    recBCH := args[2];
    exp := args[3];

    return BCH_AbstractLog_Simple_ByExponent( recLieAlg, recBCH, exp );
end;

# function to compute the time which is needed to compute logarithms
# 
# ranges .... a list of integers which give the ranges which should be tested
# no .......  number of tests which should be made to compute the average
#             runtime
#
# Example of use
# ranges := [2,4,8,16,32,64];
# 
# for BCH
# recLieAlg_bch := BCH_LieAlgebraByTGroupRec( recBCH9,exams_Fnc[7] );;
# BCH_ComputeTime_Log( recLieAlg_bch, recBCH9, ranges, 10, 1 );
#
# for NickelMat
# recLieAlg_nickelMat := MAT_LieAlgebra( [exams_Fnc[7].N, "Nickel"] );;
# BCH_ComputeTime_Log( recLieAlg_nickelMat, recBCH9, ranges, 10, 2 );
#
BCH_ComputeTime_Log  := function( recLieAlg, recBCH, ranges, no, method )
    local results, domain, times,exp, time,sum, average,range,i,hl;

    results := [];
    if method = 1 then 
        hl := HirschLength( recLieAlg.recTGroup.NN );
    elif method = 2  then
       hl := HirschLength( recLieAlg.N );
    fi;

    for range in ranges do 
        domain := [-range..range];
        times := [];
        for i in [1..no] do
            exp := List( [1..hl], x -> Random( domain ) );
            if method = 1 then 
                time := POL_CompleteRuntime2(
                            BCH_AbstractLog_Simple_ByExponent_RuntimeVersion,
                           [recLieAlg, recBCH, exp ]  ).time;
            elif method = 2 then  
                time := POL_CompleteRuntime2(
                            MAT_PcpExp2LieAlg_TestVersion ,
                           [exp, recLieAlg ]  ).time;
            elif method = 3 then 
                time := POL_CompleteRuntime2(
                            BCH_Logarithm_Symbolic,
                            [recLieAlg, exp ]  ).time;
            fi;
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


BCH_Abstract_Exponential_ByVector_RuntimeVersion := function( args )
    local recBCH, recLieAlg, vec;
    recBCH := args[1];
    recLieAlg := args[2];
    vec := args[3];
    return BCH_Abstract_Exponential_ByVector( recBCH, recLieAlg, vec );
end;

# Example of use 
# exams_Fnc := BCH_Get_FNG_TGroupRecords( 2, 9 );;
#
# ranges := [2,4,8,16,32,64];
# 
# for BCH
# recLieAlg_bch := BCH_LieAlgebraByTGroupRec( recBCH9,exams_Fnc[7] );;
# BCH_ComputeTime_Exp( recLieAlg_bch, recBCH9, ranges, 10, 1 );
#
# for NickelMat
# recLieAlg_nickelMat := MAT_LieAlgebra( [exams_Fnc[7].N, "Nickel"] );;
# BCH_ComputeTime_Exp( recLieAlg_nickelMat, recBCH9, ranges, 10, 2 );
#
BCH_ComputeTime_Exp  := function( recLieAlg, recBCH, ranges, no, method )
    local results, domain1, domain2, times,vec, time,sum, average,range,i,hl;

    results := [];
    if method = 1 then 
        hl := HirschLength( recLieAlg.recTGroup.NN );
    elif method = 2  then
       hl := HirschLength( recLieAlg.N );
    fi;

    for range in ranges do 
        domain1 := [-range..range];
        domain2 := [1..range];
        times := [];
        for i in [1..no] do
            vec := List( [1..hl], x-> Random( domain1 )/Random( domain2) );
            if method = 1 then 
                time := POL_CompleteRuntime2(
                            BCH_Abstract_Exponential_ByVector_RuntimeVersion,
                            [ recBCH, recLieAlg, vec ] ).time;
            elif method = 2 then  
                # does not work in because exp(vec) does not have to lie 
                # in the group.
                time := POL_CompleteRuntime2(
                            MAT_LieAlg2Pcp_TestVersion,
                            [vec, recLieAlg ]  ).time;
            fi;
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



BCH_Abstract_Exponential_ByElm_RunTimeVersion := function( args )
    local recBCH, recLieAlg,x;
    recBCH := args[1];
    recLieAlg := args[2];
    x := args[3];
    return BCH_Abstract_Exponential_ByElm(recBCH,recLieAlg,x ); 
end;


# Example of use 
# exams_F2c := BCH_Get_FNG_TGroupRecords( 2, 9 );;
# exams_F3c := BCH_Get_FNG_TGroupRecords( 3, 6 );;
#
# ranges := [2,4,8,16,32,64];
# ranges := [1024, 2048, 4096];
# ranges := [1024];
# 
# for BCH
# recLieAlg_bch := BCH_LieAlgebraByTGroupRec( recBCH9,exams_Fnc[8] );;
# BCH_ComputeTime_ExpAndLog( recLieAlg_bch, recBCH9, ranges, 10, 1 );
#
# for NickelMat
# recLieAlg_nickelMat := MAT_LieAlgebra( [exams_Fnc[8].N, "Nickel"] );;
# BCH_ComputeTime_ExpAndLog( recLieAlg_nickelMat, recBCH9, ranges, 10, 2 );
#
BCH_ComputeTime_ExpAndLog  := function( recLieAlg, recBCH, ranges, no, method )
    local results, domain, times,vec, time_log,time_exp,sum, 
          average_log,average_exp,range,i,hl,log;

    results := [];
    if (method = 1 or method = 3) then 
        hl := HirschLength( recLieAlg.recTGroup.NN );
    elif method = 2  then
       hl := HirschLength( recLieAlg.N );
    fi;

    for range in ranges do 
        domain := [-range..range];
        times := [[],[]];
        for i in [1..no] do
            exp := List( [1..hl], x-> Random( domain ));
            if method = 1 then 
                time_log := POL_CompleteRuntime2(
                              BCH_AbstractLog_Simple_ByExponent_RuntimeVersion,
                              [recLieAlg, recBCH, exp ]  ).time;
                log := BCH_AbstractLog_Simple_ByExponent_RuntimeVersion(
                              [recLieAlg, recBCH, exp ]  );
                time_exp := POL_CompleteRuntime2(
                              BCH_Abstract_Exponential_ByElm_RunTimeVersion,
                              [recBCH, recLieAlg, log] ).time;
            elif method = 2 then  
                time_log := POL_CompleteRuntime2(
                             MAT_PcpExp2LieAlg_TestVersion ,
                             [exp, recLieAlg ]  ).time;
                log := MAT_PcpExp2LieAlg_TestVersion([exp, recLieAlg] );
                time_exp := POL_CompleteRuntime2(
                             MAT_LieAlg2Pcp_TestVersion,
                            [log, recLieAlg ]  ).time; 
            elif method = 3 then  
                time_log := POL_CompleteRuntime2(
                             BCH_Logarithm_Symbolic ,
                             [recLieAlg, exp ]  ).time;
                log :=  BCH_Logarithm_Symbolic([recLieAlg,exp] );
                time_exp := POL_CompleteRuntime2(
                             BCH_Exponential_Symbolic,
                            [ recLieAlg,log ]  ).time; 
            fi;
            
            Add( times[1], time_log );
            Add( times[2], time_exp );
        od;
        # compute average
        sum := Sum( times[1] );
        Print( "times log : " , times[1], "\n" );
        average_log :=  Int( sum/no ) ;
        Print( "average log: ", average_log, "\n" );
        sum := Sum( times[2] );
        Print( "times exp : " , times[2], "\n" );
        average_exp :=  Int( sum/no ) ;
        Print( "average exp: ", average_exp, "\n" );
        Add( results, [range, average_log,average_exp] );
    od;

    return results;
end;


# function to compute the time which is needed to compute logarithms
# and exponentials for several lie algebras
# 
# ranges .... a list of integers which give the ranges which should be tested
# no .......  number of tests which should be made to compute the average
#             runtime
#
# Example of use:
# 
# Get the groups F_21 up to F_29
# exams_F2c := BCH_Get_FNG_TGroupRecords( 2, 9 );;
# exams_F3c := BCH_Get_FNG_TGroupRecords( 3, 5 );;
# ranges := [2,4,8,16,32,64,128];
# ranges := [64, 128, 256, 512,1024,2048,4096];
# ranges := [1024,2048,4096]; 
# ranges := [1024,2048];
#
# For BCH
# recLieAlgs_bch_F2c := List( [2..Length( exams_F2c )], x-> BCH_LieAlgebraByTGroupRec( recBCH9,exams_F2c[x] ));;
# res_bch_F2c := BCH_ComputeTime_ExpAndLog_List(recLieAlgs_bch_F2c,2,5,recBCH9,ranges,200,1);  
#  
#
# recLieAlgs_bch_F3c := List( [2..Length( exams_F3c )], x-> BCH_LieAlgebraByTGroupRec( recBCH9,exams_F3c[x] ));;
# res_bch_F3c := BCH_ComputeTime_ExpAndLog_List(recLieAlgs_bch_F3c,2,5,recBCH9,ranges,200,1); 
#
# For NickelMat
# recLieAlgs_nickelMat_Fnc := List( [2..9], x-> MAT_LieAlgebra( [exams_Fnc[x].N, "Nickel"] ) );;
# res_mat := BCH_ComputeTime_ExpAndLog_List(recLieAlgs_nickelMat_Fnc,2,7,recBCH9,ranges,10,2);  
# 
#[ [ "example 2", [ [ 1024, 0, 0 ], [ 2048, 0, 0 ], [ 4096, 0, 0 ] ] ],
#  [ "example 3", [ [ 1024, 1, 1 ], [ 2048, 1, 1 ], [ 4096, 1, 1 ] ] ],
#  [ "example 4", [ [ 1024, 4, 4 ], [ 2048, 4, 5 ], [ 4096, 4, 5 ] ] ],
#  [ "example 5", [ [ 1024, 45, 51 ], [ 2048, 46, 54 ], [ 4096, 48, 61 ] ] ],
#  [ "example 6", [ [ 1024, 210, 244 ], [ 2048, 233, 259 ], [ 4096, 237, 272 ]
#         ] ],
#  [ "example 7", [ [ 1024, 1653, 1842 ], [ 2048, 1797, 1989 ], [ 4096, 1841,
#              2027 ] ] ],
#  [ "example 8", [ [ 1024, 6217, 6772 ], [ 2048, 6618, 7201 ],
#          [ 4096, 6695, 7291 ] ] ] ]
#
# recLieAlgs_nickelMat_F3c := List( [2..Length( exams_F3c )], x-> MAT_LieAlgebra( [exams_F3c[x].N, "Nickel"] ) );;
# res_mat_F3c := BCH_ComputeTime_ExpAndLog_List(recLieAlgs_nickelMat_Æ3c,2,5,recBCH9,ranges,200,2);  
# 
#
#
# For deGraaf_Nickel
# recLieAlgs_deGraafMat_Fnc := List( [2..7], x-> MAT_LieAlgebra( [exams_Fnc[x].N, "deGraaf_Nickel"] ) );;
# res_mat_deGraaf := BCH_ComputeTime_ExpAndLog_List(recLieAlgs_deGraafMat_Fnc,2,5,recBCH9,ranges,10,2); 
#
#
# Triangular examples:
# exams_unitr_2 := BCH_Get_Unitriangular_TGroupRecords( 10, 2 );;
# 
# recLieAlgs_nickelMat_unitr_2 := List( [2..Length( exams_unitr_2 )-3], x-> MAT_LieAlgebra( [exams_unitr_2[x].N, "Nickel"] ) );;
# res_mat_unitr_2 := BCH_ComputeTime_ExpAndLog_List(recLieAlgs_nickelMat_unitr_2,2,5,recBCH9,ranges,200,2); 
# 
#
## recLieAlgs_bch_unitr_2 := List( [2..Length( exams_unitr_2)-3], x-> BCH_LieAlgebraByTGroupRec( recBCH9,exams_unitr_2[x] ));;
# res_bch_unitr_2 := BCH_ComputeTime_ExpAndLog_List(recLieAlgs_bch_unitr_2,2,5,recBCH9,ranges,200,1);  
#  
# Symbolically
# i := 0;
# for i in [1..6] do BCH_AddStarLogAndExpPols( [recLieAlgs_bch_unitr_2[i], recBCH9] ); od; 
# res_bch_symbolic_unitr_2 := BCH_ComputeTime_ExpAndLog_List(recLieAlgs_bch_unitr_2,6,6,recBCH9,ranges,200,3);  


BCH_ComputeTime_ExpAndLog_List 
           := function( recLieAlgs, start, stop, recBCH, ranges, no, method )
    local results, recLieAlg, i,result, string;
    results := [];
    for i in [start..stop] do 
        recLieAlg := recLieAlgs[i];
        result := BCH_ComputeTime_ExpAndLog( recLieAlg, 
                                              recBCH, ranges, no, method );
        string := Concatenation( "example ", String( i ) );
        Add( results, [string, result ] );
    od;
    return results;

end;


BCH_Test_ExpOfLogByTGroupList := function( recBCH, ll, start, noTests, range )
    local n,i,recLieAlg;
    n := Length( ll );
    for i in [start..n] do
        recLieAlg := BCH_SetUpLieAlgebraRecordByMalcevbasis( ll[i] );
        BCH_ComputeStructureConstants( [recBCH, recLieAlg] );
        BCH_Test_ExpOfLog ( recBCH, recLieAlg, noTests,range );
        Print( i );
    od;
    return 0;
end;
