#############################################################################
##
#W collec.gi              GUARANA package                     Bjoern Assmann
##
## Methods for collection in CN.
## Methods for collection in G using the collection in CN
## as a black box algorithm.
##
#H  @(#)$Id$
##
#Y 2006
##
##

#############################################################################
##
#F GUARANA.CutExpVector( malcevRec, exp )
##
## IN
## malcevRec ............................ malcev record
## exp .................................. exponent vector of an elment
##                                        g in G.
##
## OUT
## [exp_f,exp_c,exp_n] where exp is the concantenation of these and
## exp_f specifies the exponents of the finite part G/CN
## exp_fa specifies the exponents of the free abelian part CN/N
## exp_n specifies the exponents of the nilpotent part  N.
##
GUARANA.CutExpVector := function( malcevRec, exp )
    local indeces, exp_f, exp_fa, exp_n;
    indeces := malcevRec.indeces;
    exp_f := exp{indeces[1]};
    exp_fa := exp{indeces[2]};
    exp_n := exp{indeces[3]};
    return [exp_f, exp_fa, exp_n];
end;

#############################################################################
##
#F GUARANA.RandomGrpElm( args )
##
## IN
## args[1]=malcevRec ............... malcev record
## args[2]=range ................... range of random element
## args[4]=form .................... string that is either 
##                                   "exp" or "elm" depending on 
##                                   how the random element should be
##                                   given. 
## args[3]=grp ..................... optional string specifying whether
##                                   and elment in G,CN or N should be given.
##                                   Default is "G".
##
## OUT
## Random group elment of N,CN or G given by a exponent vector or 
## as actual group element
##
GUARANA.RandomGrpElm := function( args )
    local malcevRec, range, form, grp, hl, domain, indeces, rels, exp, i;

    malcevRec := args[1];
    range := args[2];
    if IsBound( args[3]  ) then 
	grp := args[3];
    else
	grp := "G";
    fi;
    if IsBound( args[4] ) then 
        form := args[4];
    else
	form := "exp";
    fi;

    hl := HirschLength( malcevRec.G  );
    domain := [-range..range];
    indeces := malcevRec.indeces;
    rels := RelativeOrdersOfPcp( Pcp( malcevRec.G ));

    exp := List( [1..hl], x-> 0 );
    if grp = "N" then 
	for i in indeces[3] do
	    exp[i] := Random( domain );
	od;
    elif grp = "CN" then 
	for i in Concatenation( indeces[2], indeces[3] ) do
	    exp[i] := Random( domain );
	od;
    elif grp = "G" then 
	for i in Concatenation( indeces[2], indeces[3] ) do
	    exp[i] := Random( domain );
	od;
	for i in indeces[1] do
	    exp[i] := Random( domain ) mod rels[i];
	od;
    else
	Error( " " );
    fi;
	    
    if form = "exp" then 
	return exp;
    else 
	return GUARANA.GrpElmByExpsAndPcs( exp );
    fi;
end;

## Example
## ll := GUARANA.SomePolyMalcevExams( 3 );
## R := GUARANA.InitialSetupCollecRecord( ll );
## GUARANA.AddCompleteMalcevInfo( R );


## IN
## malcevRec ........................ malcev record
## exp_g,exp_h  ..................... exponent vectors of group elements
##                                    g,h in CN 
## 
## OUT 
## The exponent vector of g*h
##
GUARANA.Collection_CN := function( malcevRec, exp_g, exp_h )

    # test whether g,h in CN

    

end;

# 
# collection functions
# -computations with powers of autom of N. 
#  Two possibilities, n given by group exponent vector or by 
#  Lie algebra coeff vector
# -collection in CN
# -computations with consecutive powers of automorphisms. 
#  again n given in two ways. 
# -Powering in C. 
# -Inversion in CN
# -Powering in CN
# -Collection in G
#
# construct more examples
# 

#############################################################################
##
#E
