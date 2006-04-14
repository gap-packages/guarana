Read( "read.g" );

# get some Tgroup records of free nilpotent groups
ll := BCH_Get_FNG_TGroupRecords( 2, 5 );; 

# set up Lie algebra record and compute structure constants for F_2,4
i := 5;
recLieAlg := BCH_SetUpLieAlgebraRecordByMalcevbasis( ll[i] );
BCH_ComputeStructureConstants( [recBCH9, recLieAlg] );


# Use Lie algebras for collection:
Ex := List( [1..8], x-> SC_Exams(x) );
FCR_bch := [];
FCR_mat := [];


# work with an individual example
i := 5;
FCR_bch[i] := BCH_FastConjugationRec( Ex[i].G, Ex[i].N, recBCH9 );
FCR_mat[i] := MAT_FastConjugationRec( Ex[i].G, Ex[i].N);

# get all examples
FCR_bch := List( [1..8], x-> BCH_FastConjugationRec(Ex[x].G, Ex[x].N, recBCH9));
FCR_mat := List( [1..8], x-> MAT_FastConjugationRec(Ex[x].G, Ex[x].N));


# compare performance with malcev
SC_ComputeAverageRuntimeFastMultiplication( FCR_bch[i], 100, 10 );
SC_ComputeAverageRuntimeFastMultiplication( FCR_mat[i], 100, 10 );

# test if collection works fine
SC_CompareAverageRuntimeCFTLVersusMalcev( FCR_bch[i], 2, 100 );


g1 := Random( FCR5.recNewParent.GG );
g2 := Random( FCR5.recNewParent.GG );
BCH_FastMultiplicationAbelianSemiTGroup ( FCR5, recBCH9, g1, g2 );

times := [];
for i in [1..8] do 
   l1 :=  SC_MalcevRuntimeByRange( FCR_bch[i], 500, 2500 );
   l2 :=  SC_MalcevRuntimeByRange( FCR_mat[i], 500, 2500 );
   Add( times, l1 );
   Add( times, l2 );
od;


# work with bigger groups
x := Indeterminate( Rationals );
pol :=  x^3 - x^2 + 4;
R := PresentTriang( 8, pol );
Ex_big := SC_Exams_Help1( R );
FCR_big :=  BCH_FastConjugationRec( Ex_big.G, Ex_big.N, recBCH9 );
g1 := Random( FCR_big.recNewParent.GG );
g2 := Random( FCR_big.recNewParent.GG );



R9 := PresentTriang( 9, pol );
Ex_big9 := SC_Exams_Help1( R9 );
FCR_big9 :=  BCH_FastConjugationRec( Ex_big9.G, Ex_big9.N, recBCH9 );
SC_MalcevRuntimeByRange( FCR_big9, 2, 1024 );

R10 := PresentTriang( 10, pol );;
Ex_big10 := SC_Exams_Help1( R10 );;
FCR_big10 :=  BCH_FastConjugationRec( Ex_big10.G, Ex_big10.N, recBCH9 );;
time;
SC_MalcevRuntimeByRange( FCR_big10, 2, 1024 );


# symbolic collection ( assume FCR is given )
n := HirschLength( FCR.recNewParent.NN );
vars_x := BCH_RationalVariableList( n, "x" );
vars_y := BCH_RationalVariableList( n, "y" );





