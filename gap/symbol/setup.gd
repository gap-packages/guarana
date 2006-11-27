#############################################################################
##
#W setup.gi               GUARANA package                     Bjoern Assmann
##
## Methods for setting up the data structures that are needed for 
## symbolic collection ( in CN ).
##
#H  @(#)$Id$
##
#Y 2006
##
##

#############################################################################
##
## global flag that is used for computations with symbolic 
## group elements whose indeterminants lie in an extension field.
## The crucial point is that these indeterminants don't interact 
## smoothly in GAP with polynomials that describe log, exp and star,
## since these have indeterminants over the rationals.
##
GUARANA.COMP_OVER_EXT_FIELD := false;

DeclareGlobalFunction( "AddSymbolicCollector" );

#############################################################################
##
#E
