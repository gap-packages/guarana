#
gap> LoadPackage("fga", false);
true
gap> exam := GUARANA.SomePolyMalcevExams( 1 );;
gap> malCol := MalcevCollectorConstruction( exam );
<<Malcev collector>>
  F : [ 2, 2, 2 ]
  C : <<Malcev object of dimension 3>>
  N : <<Malcev object of dimension 6>>
gap> GUARANA.Tests_G_Collection( malCol, 2 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 
gap> GUARANA.Tests_G_Inversion( malCol, 20 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 

#
gap> exam := GUARANA.SomePolyMalcevExams( 2 );;
gap> malCol := MalcevCollectorConstruction( exam );
<<Malcev collector>>
  F : [ 2, 2, 2, 2 ]
  C : <<Malcev object of dimension 4>>
  N : <<Malcev object of dimension 12>>
gap> GUARANA.Tests_G_Collection( malCol, 2 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 
gap> GUARANA.Tests_G_Inversion( malCol, 20 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 

#
gap> exam := GUARANA.SomePolyMalcevExams( 3 );;
gap> malCol := MalcevCollectorConstruction( exam );
<<Malcev collector>>
  F : [ 2, 2, 2, 2, 2 ]
  C : <<Malcev object of dimension 5>>
  N : <<Malcev object of dimension 20>>
gap> GUARANA.Tests_G_Collection( malCol, 2 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 
gap> GUARANA.Tests_G_Inversion( malCol, 20 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 

#
gap> exam := GUARANA.SomePolyMalcevExams( 4 );;
gap> malCol := MalcevCollectorConstruction( exam );
<<Malcev collector>>
  F : [ 2, 2, 2 ]
  C : <<Malcev object of dimension 3>>
  N : <<Malcev object of dimension 9>>
gap> GUARANA.Tests_G_Collection( malCol, 2 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 
gap> GUARANA.Tests_G_Inversion( malCol, 20 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 

#
gap> exam := GUARANA.SomePolyMalcevExams( 5 );;
gap> malCol := MalcevCollectorConstruction( exam );
<<Malcev collector>>
  F : [ 2, 2, 2, 2 ]
  C : <<Malcev object of dimension 4>>
  N : <<Malcev object of dimension 18>>
gap> GUARANA.Tests_G_Collection( malCol, 2 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 
gap> GUARANA.Tests_G_Inversion( malCol, 20 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 

#
gap> exam := GUARANA.SomePolyMalcevExams( 6 );;
gap> malCol := MalcevCollectorConstruction( exam );
<<Malcev collector>>
  F : [  ]
  C : <<Malcev object of dimension 1>>
  N : <<Malcev object of dimension 8>>
gap> GUARANA.Tests_G_Collection( malCol, 2 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 
gap> GUARANA.Tests_G_Inversion( malCol, 20 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 

#
gap> exam := GUARANA.SomePolyMalcevExams( 7 );;
gap> malCol := MalcevCollectorConstruction( exam );
<<Malcev collector>>
  F : [  ]
  C : <<Malcev object of dimension 1>>
  N : <<Malcev object of dimension 14>>
gap> GUARANA.Tests_G_Collection( malCol, 2 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 
gap> GUARANA.Tests_G_Inversion( malCol, 20 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 

#
gap> exam := GUARANA.SomePolyMalcevExams( 8 );;
gap> malCol := MalcevCollectorConstruction( exam );
<<Malcev collector>>
  F : [  ]
  C : <<Malcev object of dimension 1>>
  N : <<Malcev object of dimension 32>>
gap> GUARANA.Tests_G_Collection( malCol, 2 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 
gap> GUARANA.Tests_G_Inversion( malCol, 20 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 

#
gap> exam := GUARANA.SomePolyMalcevExams( 9 );;
gap> malCol := MalcevCollectorConstruction( exam );
<<Malcev collector>>
  F : [ 2, 2, 2 ]
  C : <<Malcev object of dimension 4>>
  N : <<Malcev object of dimension 12>>
gap> GUARANA.Tests_G_Collection( malCol, 2 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 
gap> GUARANA.Tests_G_Inversion( malCol, 20 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 

#
gap> exam := GUARANA.Tr_n_O1( 1 );;
gap> malCol := MalcevCollectorConstruction( exam );
<<Malcev collector>>
  F : [ 2 ]
  C : <<Malcev object of dimension 1>>
  N : <<Malcev object of dimension 0>>
gap> GUARANA.Tests_G_Collection( malCol, 2 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 
gap> GUARANA.Tests_G_Inversion( malCol, 20 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 

#
gap> exam := GUARANA.Tr_n_O2( 1 );;
gap> malCol := MalcevCollectorConstruction( exam );
<<Malcev collector>>
  F : [ 2 ]
  C : <<Malcev object of dimension 1>>
  N : <<Malcev object of dimension 0>>
gap> GUARANA.Tests_G_Collection( malCol, 2 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 
gap> GUARANA.Tests_G_Inversion( malCol, 20 );
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 
56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 
82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 
