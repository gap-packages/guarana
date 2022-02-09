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
                n := n,
                exps := "unknown yet" );
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
                n := n,
                exps := "unknown yet" );
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
                n := n,
                exps := exps_cn );
    return Objectify( malCol!.cn_elms_type , elm );
end);

InstallGlobalFunction( MalcevCNElementByExponentsNC,
function( malCol, exps_cn )
    local hlCN_N, hlC, hlN, hlCN, exps_n, exps_c, c, n, elm;

    # setup
    hlCN_N := Length( malCol!.indeces[2] );
    hlC := HirschLength( malCol!.CC );
    hlN := HirschLength( malCol!.NN ); 
    hlCN := hlCN_N + hlN;

    # get exps_c, exps_n 
    exps_n := exps_cn{ [1+hlCN_N..hlCN] };
    exps_c := List( [1..hlC], x-> 0 );
    exps_c{[1..hlCN_N]} := exps_cn{[1..hlCN_N]};

    # get malcev element of C and N
    c := MalcevGenElementByExponents( malCol!.mo_CC, exps_c );
    n := MalcevGenElementByExponents( malCol!.mo_NN, exps_n );

    elm := rec( malCol := malCol,
                c := c,
                n := n,
                exps := exps_cn );
    return Objectify( malCol!.cn_elms_type , elm );
end);

InstallGlobalFunction( MalcevCNElementBy2GenElements,
function( malCol, c, n )
    local elm;
    elm := rec( malCol := malCol,
                c := c,
                n := n,
                exps := "unknown yet" );
    return Objectify( malCol!.cn_elms_type , elm );
end);

InstallGlobalFunction( MalcevSymbolicCNElementByExponents,
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
    c := MalcevSymbolicGenElementByExponents( malCol!.mo_CC, exps_c );
    n := MalcevSymbolicGenElementByExponents( malCol!.mo_NN, exps_n );

    elm := rec( malCol := malCol,
                c := c,
                n := n,
                exps := exps_cn );
    return Objectify( malCol!.cn_elms_type , elm );
end);

GUARANA.CN_Identity := function( malCol )
    local hlCN, exps;
    hlCN := malCol!.lengths[2] + malCol!.lengths[3];
    exps := List( [1..hlCN], x-> 0 );
    return MalcevCNElementByExponents( malCol, exps );
end;


#############################################################################
##
## Methods for constructing Malcev G Elements
##
InstallGlobalFunction( MalcevGElementByExponents,
function( malCol, exps )
    local n, indeces, rels, exps_f, exps_cn, cn_elm, elm, i,r;

    # check input
    n := Sum( malCol!.lengths );
    if n <> Length( exps ) then
        Error( "Wrong length of exponent vector\n" );
    fi;

    indeces := malCol!.indeces;
    rels := RelativeOrdersOfPcp( Pcp( malCol!.G ) );
    exps_f := [];
    for i in indeces[1] do
        r := [0..rels[i]-1];
        if not exps[i] in r then 
            Error( "Input vector out of range" );
        fi;
        exps_f[i] := exps[i];
    od;

    exps_cn := exps{ Concatenation(indeces[2],indeces[3]) };
    cn_elm := MalcevCNElementByExponents( malCol, exps_cn );

    elm := rec( malCol := malCol,
                cn_elm := cn_elm,
                exps_f := exps_f,
                exps := exps );
    return Objectify( malCol!.g_elms_type , elm );
end);

InstallGlobalFunction( MalcevGElementByExponentsNC,
function( malCol, exps )
    local exps_cn, cn_elm, indeces, exps_f, elm;

    indeces := malCol!.indeces;
    exps_f := exps{indeces[1]};
    exps_cn := exps{ Concatenation(indeces[2],indeces[3]) };
    cn_elm := MalcevCNElementByExponentsNC( malCol, exps_cn );

    elm := rec( malCol := malCol,
                cn_elm := cn_elm,
                exps_f := exps_f,
                exps := exps );
    return Objectify( malCol!.g_elms_type , elm );
end);

InstallGlobalFunction( MalcevGElementByCNElmAndExps,
function( malCol, exps_f, cn_elm )
    local elm;
    elm := rec( malCol := malCol,
                cn_elm := cn_elm,
                exps_f := exps_f,
                exps := "unknown yet" );
    return Objectify( malCol!.g_elms_type , elm );
end);

GUARANA.G_Identity := function( malCol )
    local hlG_CN, exps_f, cn_elm;
    hlG_CN := malCol!.lengths[1];
    exps_f := List( [1..hlG_CN], x-> 0 );
    cn_elm := GUARANA.CN_Identity( malCol );
    return MalcevGElementByCNElmAndExps( malCol, exps_f, cn_elm );
end;

#############################################################################
##
## Methods for creating random elements

GUARANA.RandomCNElement := function( malCol, range )
    local hl, exps, a;

    # get HirschLength of CN
    hl := malCol!.lengths[2] + malCol!.lengths[3];

    exps := List( [1..hl], x-> Random( [-range..range] ) );

    a := MalcevCNElementByExponents( malCol, exps );
    return a;
end;

GUARANA.RandomCNElement2 := function( malCol, range )
    local hlC, hlN, exps_c, exps_n;

    hlC := HirschLength( malCol!.CC );
    hlN := HirschLength( malCol!.NN );

    exps_c := List( [1..hlC], x-> Random( [-range..range] ) );
    exps_n := List( [1..hlN], x-> Random( [-range..range] ) );

    return MalcevCNElementBy2Exponents( malCol, exps_c, exps_n );
end;

GUARANA.RandomGElement := function( malCol, range )
    local rels, domain, exps, i, indeces;

    rels := RelativeOrdersOfPcp( Pcp( malCol!.G ) );
    domain := [-range..range];
    exps := [];
    indeces := malCol!.indeces;

    for i in indeces[1] do
        exps[i] := Random( domain ) mod rels[i];
    od;
    for i in Concatenation( indeces[2], indeces[3] ) do
        exps[i] := Random( domain );
    od;

    return MalcevGElementByExponents( malCol, exps );
end;

InstallOtherMethod( Random, 
               "for Malcev Collectors (Guarana)", 
               true, 
               [IsMalcevCollectorRep ], 
               0,
function( malCol )
    local range;
    range := 10;
    return Random( malCol, range ) ;
end);

InstallOtherMethod( Random, 
               "for Malcev Collectors and Integers (Guarana)", 
               true, 
               [IsMalcevCollectorRep, IsInt ], 
               0,
function( malCol, range )
    return Random( malCol, "G", range );
end);

InstallOtherMethod( Random, 
               "for Malcev Collectors and strings (Guarana)", 
               true, 
               [IsMalcevCollectorRep, IsString ], 
               0,
function( malCol, info )
    local range;
    range := 10;
    return Random( malCol, info, range );
end);

InstallOtherMethod( Random, 
               "for Malcev Collectors, strings and integers (Guarana)", 
               true, 
               [IsMalcevCollectorRep, IsString, IsInt ], 
               0,
function( malCol, info, range )
    if info = "CN" then 
        return GUARANA.RandomCNElement( malCol, range );
    elif info = "CN2" then 
        return GUARANA.RandomCNElement2( malCol, range );
    elif info = "G" then
        return GUARANA.RandomGElement( malCol, range );
    fi;
end);

#############################################################################
##
#M Print Malcev CN elements
##
InstallMethod( PrintObj, 
               "for Malcev CN elements (Guarana)", 
               true, 
               [IsMalcevCNElement ], 
               0,
function( elm )
    Print( "c \n", elm!.c, "\n" );
    Print( "n \n", elm!.n, "\n" );
    Print( "Exponents: ", elm!.exps );
end );

#############################################################################
##
#M Print Malcev G elements
##
InstallMethod( PrintObj, 
               "for Malcev G elements (Guarana)", 
               true, 
               [IsMalcevGElement ], 
               0,
function( elm )
    Print( Exponents( elm ) );    
#    if IsString( elm!.exps ) then 
#        Print( "Exponent vector of finite part: ", elm!.exps_f, "\n" );
#        Print( "CN Element: \n", elm!.cn_elm );
#    else
#        Print(  elm!.exps );
#    fi;
end );

#############################################################################
##
## Computing normal forms for CN elements

## g = c n = c_1^x_1 ... c_k^x_k * tail * n
## c(g) = c_1^x_1 ... c_k^x_k 
## Normal form is with respect to the pcs of CN.
## g = c(g)n(g) 
##
##
## COMMENT
## You could also divide off in the Lie algebra 
## (Pols single versus generic would be suitable for that ).
##
GUARANA.CN_NormalForm := function( g )
    local malCol, exps_c, hlC, hlCN_N, exps_c_g, exps_tail, tail, 
          tail2, n_new, c_new, g_new;

    # catch case when element is known to be normalised
    malCol := g!.malCol;
    if not IsString( g!.exps ) then 
        return g; 
    fi;

    # get exponents of c(g) with respect to the pcs of C
    exps_c := Exponents( g!.c );
    hlC := HirschLength( malCol!.CC );
    hlCN_N := Length( malCol!.indeces[2] );
    exps_c_g := List( [1..hlC], x-> 0 );
    exps_c_g{[1..hlCN_N]} := exps_c{[1..hlCN_N]};

    # get tail in C
    exps_tail := List( [1..hlC], x-> 0 );
    exps_tail{[hlCN_N+1..hlC]} := exps_c{[hlCN_N +1..hlC]};
    if exps_tail = 0* exps_tail then 
        # elm is already in normal form
        g!.exps := Concatenation( exps_c_g{[1..hlCN_N]},
                                  Exponents( g!.n ) );
        return g;
    fi;
    tail := MalcevGenElementByExponents( malCol!.mo_CC, exps_tail );

    # map it to N
    tail2 := GUARANA.MapFrom_MOC_2_MON( malCol, tail );

    # multiply with n 
    n_new := tail2 * g!.n;

    c_new := MalcevGenElementByExponents( malCol!.mo_CC, exps_c_g );

    g_new := MalcevCNElementBy2GenElements( malCol, c_new, n_new ); 
    g_new!.exps := Concatenation( exps_c_g{[1..hlCN_N]},
                                  Exponents( g_new!.n ) );
    return g_new;
end;

InstallMethod( NormalForm, 
               "for Malcev CN elements (Guarana)", 
               true, 
               [IsMalcevCNElement ], 
               0,
function( elm )
    return GUARANA.CN_NormalForm( elm );
end );

InstallMethod( Normalise, 
               "for Malcev CN elements (Guarana)", 
               true, 
               [IsMalcevCNElement ], 
               0,
function( elm )
    local malCol, nf;

    # catch case when element is known to be normalised
    malCol := elm!.malCol;
    if not IsString( elm!.exps ) then 
        return elm; 
    fi;

    # compute normal form
    nf := GUARANA.CN_NormalForm( elm );

    # update elm
    elm!.c := nf!.c;
    elm!.n := nf!.n;
    elm!.exps := nf!.exps;
end );


## Note that this changes the input
InstallOtherMethod( Exponents, 
               "for Malcev CN elements (Guarana)", 
               true, 
               [IsMalcevCNElement ], 
               0,
function( g )
    if IsString( g!.exps ) then 
        Normalise( g );
        return g!.exps;
    else 
        return g!.exps;
    fi;
end);

InstallOtherMethod( Exponents, 
               "for Malcev G elements (Guarana)", 
               true, 
               [IsMalcevGElement ], 
               0,
function( g )
    local exps_f, exps_cn;
    if g!.exps = [] then
        return g!.exps;
    elif IsString( g!.exps ) then
        exps_f := g!.exps_f;
        exps_cn := Exponents( g!.cn_elm);
        g!.exps := Concatenation( exps_f, exps_cn );
        return g!.exps;
    else
        return g!.exps;
    fi;
end);

#############################################################################
##
#M x = y 
##
InstallOtherMethod( \=, 
               "for Malcev CN elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevCNElement, IsMalcevCNElement ],
		0, 
function( x, y )
    return Exponents( x ) = Exponents( y );
end);

InstallOtherMethod( \=, 
               "for Malcev G elments (Guarana)",
	       IsIdenticalObj,
	        [IsMalcevGElement, IsMalcevGElement ],
		0, 
function( x, y )
    return Exponents( x ) = Exponents( y );
end);

#############################################################################
##
## Methods for setting the malcev collector as used collector of pcp groups
##
InstallMethod( AttacheMalcevCollector, 
               "for pcp groups and Mal'cev collectors (Guarana)", 
               true, 
               [IsPcpGroup, IsMalcevCollectorRep], 
               0,
function( G, malCol )
    local coll;
    coll := Collector( G );
    coll![GUARANA.MALCEV_COLL_STORAGEPLACE] := malCol;
end);

InstallMethod( AttachedMalcevCollector, 
               "for pcp elements (Guarana)", 
               true, 
               [IsPcpElement ], 
               0,
function( g )
    local coll;
    coll := Collector( g );
    if IsBound( coll![GUARANA.MALCEV_COLL_STORAGEPLACE] ) then 
        return coll![GUARANA.MALCEV_COLL_STORAGEPLACE]; 
    else
        return fail;;
    fi;
end);

InstallMethod( IsMalcevPcpElement, 
               "for pcp elements (Guarana)", 
               true, 
               [IsPcpElement ], 
               0,
function( g )
    local coll;
    coll := Collector( g );
    if IsBound( coll![GUARANA.MALCEV_COLL_STORAGEPLACE] ) then 
        return true;
    else
        return false;
    fi;
end);


# turn true only for testing purpose
GUARANA.UseMalcevColl := false;
if GUARANA.UseMalcevColl then 

InstallMethod( \*, "for pcp elements",
               IsIdenticalObj, 
               [IsPcpElement ,#and IsMalcevPcpElement, 
                IsPcpElement], 30,
function( g, h )
    local malCol, gg, hh, exp_res, res;
    if IsMalcevPcpElement( g ) then 
        Print( "." );
        malCol := AttachedMalcevCollector( g ); 
        gg := MalcevGElementByExponents( malCol, Exponents( g ) );
        hh := MalcevGElementByExponents( malCol, Exponents( h ) );
        exp_res := Exponents( gg*hh);
        res := PcpElementByExponents( g!.collector, exp_res );
        return res;
    else
        TryNextMethod();
    fi;
end);

InstallMethod( COMM, "for pcp elements",
               IsIdenticalObj, 
               [IsPcpElement ,#and IsMalcevPcpElement, 
                IsPcpElement], 30,
function( g, h )
    local malCol, gg, hh, exp_res, res;
    if IsMalcevPcpElement( g ) then 
        Print( "," );
        malCol := AttachedMalcevCollector( g ); 
        gg := MalcevGElementByExponents( malCol, Exponents( g ) );
        hh := MalcevGElementByExponents( malCol, Exponents( h ) );
        exp_res := Exponents( g^-1*h^-1*g*h);
        res := PcpElementByExponents( g!.collector, exp_res );
        return res;
    else
        TryNextMethod();
    fi;
end);

InstallMethod( \^, "for pcp elements",
               [IsPcpElement, IsInt], 30,
function( g, n )
    local malCol, gg, hh, exp_res, res;
    if IsMalcevPcpElement( g ) then
        Print( ":" );
        malCol := AttachedMalcevCollector( g );
        gg := MalcevGElementByExponents( malCol, Exponents( g ) );
        exp_res := Exponents( gg^n);
        res := PcpElementByExponents( g!.collector, exp_res );
        return res;
    else
        TryNextMethod();
    fi;
end);

InstallMethod( \^, "for two pcp elements",
               IsIdenticalObj,
               [IsPcpElement, IsPcpElement], 30,
function( g, h )
    local malCol, gg, hh, exp_res, res;
    if IsMalcevPcpElement( g ) then
        Print( "-" );
        malCol := AttachedMalcevCollector( g );
        gg := MalcevGElementByExponents( malCol, Exponents( g ) );
        hh := MalcevGElementByExponents( malCol, Exponents( g ) );
        exp_res := Exponents( h^-1*g*h);
        res := PcpElementByExponents( g!.collector, exp_res );
        return res;
    else
        TryNextMethod();
    fi;
end);

InstallMethod( Inverse, "for pcp elements",
               [IsPcpElement], 30,
function( g )
    local malCol, gg, hh, exp_res, res;
    if IsMalcevPcpElement( g ) then
        Print( ";" );
        malCol := AttachedMalcevCollector( g );
        gg := MalcevGElementByExponents( malCol, Exponents( g ) );
        exp_res := Exponents( Inverse(gg ));
        res := PcpElementByExponents( g!.collector, exp_res );
        return res;
    else
        TryNextMethod();
    fi;
end);

InstallMethod( InverseMutable, "for pcp elements",
               [IsPcpElement], 30,
function( g )
    local malCol, gg, hh, exp_res, res;
    if IsMalcevPcpElement( g ) then
        Print( ";" );
        malCol := AttachedMalcevCollector( g );
        gg := MalcevGElementByExponents( malCol, Exponents( g ) );
        exp_res := Exponents( Inverse(gg ));
        res := PcpElementByExponents( g!.collector, exp_res );
        return res;
    else
        TryNextMethod();
    fi;
end);

fi;


