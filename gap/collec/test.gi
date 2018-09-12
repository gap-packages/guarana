#############################################################################
##
#W test.gi              GUARANA package                     Bjoern Assmann
##
##
#H  @(#)$Id$
##
#Y 2006
##
##

GUARANA.Test_CN_Collection := function( malCol, range )
    local g, h, exp_g, exp_h, exp_gh, indeces, pcpCN, g2, h2, exp_gh2;
   
    # get random elements in CN
    g := Random( malCol, "CN", range );
    h := Random( malCol, "CN", range );
    exp_g := Exponents( g );
    exp_h := Exponents( h );

    # compute gh with Malcev
    exp_gh := Exponents( g*h );
    #Print( exp_gh );

    # compute gh with Cfts
    indeces := Concatenation( malCol!.indeces[2], malCol!.indeces[3] );
    pcpCN := Pcp( malCol!.G ){indeces};
    g2 := GUARANA.GrpElmByExpsAndPcs( pcpCN, exp_g ); 
    h2 := GUARANA.GrpElmByExpsAndPcs( pcpCN, exp_h ); 
    exp_gh2 := Exponents( g2*h2 ){indeces};
    #Print( exp_gh2 );

    return exp_gh2 = exp_gh;
end;

GUARANA.Tests_CN_Collection := function( malCol, range )
    local no, test, i;
    
    no := 100;
    for i in [1..no] do
	test := GUARANA.Test_CN_Collection( malCol, range );
	if test = true then
	    #Print( "yeah" );
	else 
	    Error( " " );
	fi;
    od;
end;

GUARANA.Test_MO_CconjF_ElmConjByFiniteElm := function( malCol, range )
    local g, ind_f, ind_c, exp_f, exp_c, c_f, exp_c_f, pcpG, indeces, pcpCN, f2, c2, c_f_2, exp_c_f_2;

    # get random exp_c, exp_f 
    g := Random( malCol, range );
    ind_f := malCol!.indeces[1];
    ind_c := malCol!.indeces[2];
    exp_f := g!.exps_f;
    exp_c := Exponents( g ){ind_c};

    # compute c^f using Malcev 
    c_f := GUARANA.MO_CconjF_ElmConjByFiniteElm( malCol, exp_c, exp_f );
    exp_c_f := Exponents( c_f );
    Print( exp_c_f );

    # compute c^f using collection from the left
    pcpG := Pcp( malCol!.G );  
    indeces := Concatenation( malCol!.indeces[2], malCol!.indeces[3] );
    pcpCN := pcpG{indeces};
    f2 := GUARANA.GrpElmByExpsAndPcs( pcpG, exp_f );
    c2 := GUARANA.GrpElmByExpsAndPcs( pcpCN, exp_c );
    c_f_2 := c2^f2;
    exp_c_f_2 := Exponents( c_f_2 ){indeces};
    Print( exp_c_f_2 );

    # compare
    return exp_c_f = exp_c_f_2;
end;

GUARANA.Test_CN_ConjugationByFiniteElm := function( malCol, range )
    local g, exp_g, n, h, exp_h, exp_f, g_f, exp_g_f, indeces, pcpCN, g2, f2, g2_f2, exp_g2_f2;

    # get random element in CN
    g := Random( malCol, "CN", range ); 
    exp_g := Exponents( g );

    # get random exp vector of finite part 
    n := malCol!.lengths[1];
    h := Random( malCol,  1 ); 
    exp_h := Exponents( h );
    exp_f := exp_h{[1..n]}; 

    # compute g^f with Malcev 
    g_f := GUARANA.MO_CN_ConjugationByFiniteElm( g, exp_f );
    exp_g_f := Exponents( g_f ); 
    Print( exp_g_f, "\n" );

    # comput g^f with Cftl
    indeces := Concatenation( malCol!.indeces[2], malCol!.indeces[3] );
    pcpCN := Pcp( malCol!.G ){indeces};
    g2 := GUARANA.GrpElmByExpsAndPcs( pcpCN, exp_g ); 
    f2 := GUARANA.GrpElmByExpsAndPcs( Pcp(malCol!.G), exp_f );
    g2_f2 := g2^f2;
    exp_g2_f2 := Exponents( g2_f2 ){indeces};
    Print( exp_g2_f2 );

    return exp_g_f = exp_g2_f2; 
end;

GUARANA.Test_G_Collection := function( malCol, range )
    local g, h, exp_g, exp_h, exp_gh, pcpG, g2, h2, exp_gh2;
   
    # get random elements in G
    g := Random( malCol, range );
    h := Random( malCol, range );
    exp_g := Exponents( g );
    exp_h := Exponents( h );

    # comput gh with Malcev
    exp_gh := Exponents( g*h );
    #Print( exp_gh, "\n" );

    # comput gh with Cfts
    pcpG := Pcp( malCol!.G );
    g2 := GUARANA.GrpElmByExpsAndPcs( pcpG, exp_g );
    h2 := GUARANA.GrpElmByExpsAndPcs( pcpG, exp_h );
    exp_gh2 := Exponents( g2*h2 );
    #Print( exp_gh2, "\n" );

    return exp_gh2 = exp_gh;
end;

GUARANA.Tests_G_Collection := function( malCol, range )
    local no, test, i;
    
    no := 100;
    for i in [1..no] do
	test := GUARANA.Test_G_Collection( malCol, range );
	if test = true then
	    Print( i, " " );
	else 
	    Error( " " );
	fi;
    od;
    Print("\n");
end;

GUARANA.Test_G_Inversion := function( malCol, range )
    local g, g_inv, g_inv_inv;
 
    # get random element in G
    g := Random( malCol, range );

    # comput g^-1
    g_inv := g^-1;

    # compute g^-1^-1
    g_inv_inv :=  g_inv^-1; 

    return g = g_inv_inv;
end;

GUARANA.Tests_G_Inversion := function( malCol, range )
    local no, test, i;
    
    no := 100;
    for i in [1..no] do
	test := GUARANA.Test_G_Inversion( malCol, range );
	if test = true then
	    Print( i, " " );
	else 
	    Error( " " );
	fi;
    od;
    Print("\n");
end;

GUARANA.Test_SmallExams := function()
    local ind, exams, malCol, rangeCftl, rangeInv, i;

    # get indeces malcev records 
    ind :=  [1..9];
    exams := [];
    for i in ind do 
	    exams[i] := GUARANA.SomePolyMalcevExams( i );
    od;
    exams[10] := GUARANA.Tr_n_O1( 1 );
    exams[11] := GUARANA.Tr_n_O2( 1 );
    Add( ind, 10 );
    Add( ind, 11 );

    for i in ind do  
	    Print( "Testing example ", i, "\n" );
	    # get malcev record
	    malCol := MalcevCollectorConstruction( exams[i] );

        # compare Malcev collection with collection from the left
	    Print( "Comparing with Cftl \n" );
	    rangeCftl := 2;
	    GUARANA.Tests_G_Collection( malCol, rangeCftl );
	    Print( "\n" );

        # check inversion 
	    Print( "Checking Inversion \n" );
	    rangeInv := 20;
	    GUARANA.Tests_G_Inversion( malCol, rangeInv );
	    Print( "\n" );
    od;
    return 0;
end;

GUARANA.Test_BigExams := function()
    local ind_Tr_n_O1, ind_Tr_n_O2, exam, malCol, rangeInv, i;
    
    # get indeces of examples 
    ind_Tr_n_O1 := [1..9];
    ind_Tr_n_O2 := [1..9];
    
    for i in ind_Tr_n_O1 do  
	    Print( "Testing Tr_n_O1 with n = ", i, "\n" );
	    # get malcev record
	    exam := GUARANA.Tr_n_O1( i );
	    malCol := MalcevCollectorConstruction( exam );

        # check inversion 
	    Print( "Checking Inversion \n" );
	    rangeInv := 20;
	    GUARANA.Tests_G_Inversion( malCol, rangeInv );
	    Print( "\n" );
    od;

    for i in ind_Tr_n_O2 do  
	    Print( "Testing Tr_n_O2 with n = ", i, "\n" );
	    # get malcev record
	    exam := GUARANA.Tr_n_O2( i );
	    malCol := MalcevCollectorConstruction( exam );

        # check inversion 
	    Print( "Checking Inversion \n" );
	    rangeInv := 20;
	    GUARANA.Tests_G_Inversion( malCol, rangeInv );
	    Print( "\n" );
    od;

    return 0;
end;

#############################################################################
##
#E
