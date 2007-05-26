
ll := GUARANA.Tr_n_O1( 3 );
malCol := MalcevCollectorConstruction( ll );


exps_g := [ 1, 1, 1, -3, -2, 1, -2, -1, 0, 3, -1,3 ];
exps_h := [ 1, 0, 1, -1, 0, 2, 0, 4, -1, 5, 9,-5 ];

g := MalcevGElementByExponents( malCol, exps_g );
h := MalcevGElementByExponents( malCol, exps_h );
k := g*h;

Random( malCol, 10 );


GUARANA.AverageRuntimeCollec( malCol, [2,4,8], 1000 );   

    
