#############################################################################
##
#W elms.gi               GUARANA package                     Bjoern Assmann
##
##
#H  @(#)$Id$
##
#Y 2006
##
##

#############################################################################
##
## Methods for constructing Malcev CN Elements
##
InstallGlobalFunction( MalcevCNElementBy2Exponents,
function( malCol, exps_c, exps_n )
    local c, n, elm;

    # get malcev element of C and N
    c := MalcevGenElementByExponents( malCol!.mo_CC, exps_c );
    n := MalcevGenElementByExponents( malCol!.mo_NN, exps_n );

    elm := rec( malCol := malCol,
                c := c,
                n := n );
    return Objectify( malCol!.cn_elms_type , elm );
end);

InstallGlobalFunction( MalcevCNElementBy2Coefficients,
function( malCol, coeffs_c, coeffs_n )
    local c, n, elm;

    # get malcev element of C and N
    c := MalcevGenElementByCoefficients( malCol!.mo_CC, coeffs_c );
    n := MalcevGenElementByCoefficients( malCol!.mo_NN, coeffs_n );

    elm := rec( malCol := malCol,
                c := c,
                n := n );
    return Objectify( malCol!.cn_elms_type , elm );
end);

InstallGlobalFunction( MalcevCNElementByExponents,
function( malCol, exps_cn )
    local hlCN_N, hlC, hlN, hlCN, exps_n, exps_c, c, n, elm;

    # setup
    hlCN_N := Length( malCol!.indeces[2] );
    hlC := HirschLength( malCol!.CC );
    hlN := HirschLength( malCol!.NN ); 
    hlCN := hlCN_N + hlN;

    # check input
    if Length( exps_cn ) <> hlCN then
        Error( "Input vector has not correct lenghts" );
    fi;

    # get exps_c, exps_n 
    exps_n := exps_cn{ [1+hlCN_N..hlCN] };
    exps_c := List( [1..hlC], x-> 0 );
    exps_c{[1..hlCN_N]} := exps_cn{[1..hlCN_N]};

    # get malcev element of C and N
    c := MalcevGenElementByExponents( malCol!.mo_CC, exps_c );
    n := MalcevGenElementByExponents( malCol!.mo_NN, exps_n );

    elm := rec( malCol := malCol,
                c := c,
                n := n );
    return Objectify( malCol!.cn_elms_type , elm );
end);

GUARANA.RandomCNElement := function( malCol, range )
    local hl, exps, a;

    # get HirschLength of CN
    hl := malCol!.lengths[2] + malCol!.lengths[3];

    range := 10;
    exps := List( [1..hl], x-> Random( [-range..range] ) );

    a := MalcevCNElementByExponents( malCol, exps );
    return a;
end;

InstallOtherMethod( Random, 
               "for Malcev CN Elements (Guarana)", 
               true, 
               [IsMalcevCollectorRep, IsString ], 
               0,
function( malCol, info )
    local range;

    range := 10;
    if info = "CN" then 
        return GUARANA.RandomCNElement( malCol, range );
    fi;
end);

#############################################################################
##
#M Print Malcev CN elements
##
InstallMethod( PrintObj, 
               "for Malcev gen elements (Guarana)", 
               true, 
               [GUAR_IsCNElement ], 
               0,
function( elm )
    Print( "c \n", elm!.c, "\n" );
    Print( "n \n", elm!.n );
end );

# Exponents ( with respect of pcs of CN )
#       GetTail ( in the intersection of C and N )
#       map it to N
#       get CN/N part
# 

#############################################################################
##
#E
