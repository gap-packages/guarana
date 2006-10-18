#############################################################################
##
#W collec.gi               GUARANA package                     Bjoern Assmann
##
## Methods for  symbolic collection ( in CN ).
##
#H  @(#)$Id$
##
#Y 2006
##
##

#TODO
# -compare with normal Mal'cev collector to see how it works
# -is there a problem with indeterminantes (in Rationals versus 
#  splitting field)?
# -we need a setter functions that swiches on, the right symbolic 
#  collection in C,N.
#  This also affects the setup of the malcev collector.
#
# 1. setup of the symbolic elements.
#    maybe in a file coalled elms.gi
#    There we should have symbolic CN elements
# 2. Then complete this file with the aim to install a "*" method
#    for those elements.
if false then 
    mc := ExamplesOfSomeMalcevCollectors(1 );   
    GUARANA.SC_FullSetup( mc ); 

    GUARANA.SC_GetGenericElms( mc );
    g := mc!.symCol.left;
    res := GUARANA.SC_ComputeCollectionFunc( mc, g );
    # gives an error
    Exponents( res ); 
    

fi;

GUARANA.SC_GetGenericElms := function( malCol )
    local left_vars, right_vars;

    if not IsBound( malCol!.symCol ) then 
        Error( "Symbolic collector is not set up." );
    fi;

    # get variables
    left_vars := malCol!.symCol.left_vars;
    right_vars := malCol!.symCol.right_vars;

    # store elms
    malCol!.symCol.left:=MalcevSymbolicCNElementByExponents(malCol,left_vars);
    malCol!.symCol.right:=MalcevSymbolicCNElementByExponents(malCol,right_vars);
    
end;


## IN
## malCol ..................... Malcev collector
##    n ....................... Malcev gen elment of N
##
## n^c given as Malcev gen elment of N where 
## c is the C-part of malCol!.symCol.right, 
## which is an internal generic element.
## Note that we use fact that c has only non-zero exponents in 
## the CN/N part.
##
GUARANA.SC_N_ConjugationByRightC_Elm := function( malCol, n )
    local coeffs;

    coeffs := Coefficients( n );

    # conjugate by semisimple action
    coeffs := coeffs * malCol!.symCol.fullSemisimpleAction;

    # conjugate by unipotent action
    coeffs := coeffs * malCol!.symCol.fullUnipotentAction;

    return MalcevSymbolicGenElementByCoefficients( malCol!.mo_NN, coeffs );
end;

##
## OUT
## Compute the collection function of g*h 
## where h is malCol!.symCol.right;
##
## Comment:
## g h = c(g)n(g) c(h)n(h)
##     = c(g)c(h0 n(g)^c(h) n(h)
##
GUARANA.SC_ComputeCollectionFunc := function( malCol, g )
    local h, c_new, n_c, n_new;

    # setup
    h := malCol!.symCol.right;

    c_new := g!.c * h!.c;
    n_c := GUARANA.SC_N_ConjugationByRightC_Elm( malCol, g!.n );
    n_new := n_c * h!.n;
    return MalcevCNElementBy2GenElements( malCol, c_new, n_new );
end;

# TODO:
# -CN elements should become symbolic (should pass the test "IsSymbolic"
# -MalcevCNElementBy2GenElements should become symbolic
#
#############################################################################
##
#E
