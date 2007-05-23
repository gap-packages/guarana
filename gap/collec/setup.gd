#############################################################################
##
#W setup2.gd               GUARANA package                     Bjoern Assmann
##
## Methods for setting up the data structures that are needed for 
## Mal'cev collection.
##
#H  @(#)$Id$
##
#Y 2006
##
##

DeclareCategory( "IsMalcevCollector", IsObject );
DeclareRepresentation( "IsMalcevCollectorRep", IsComponentObjectRep, 
                       [ "G", 
                         "indeces",
                         "lengths",
                         "N",
                         "NN",
                         "C",
                         "CC",
                         "mo_NN",
                         "mo_CC"  ] );

BindGlobal( "MalcevCollectorFamily", 
     NewFamily( "MalcevCollector", IsMalcevCollectorRep ) );

DeclareGlobalFunction( "MalcevCollectorConstruction"); 

#############################################################################
##
#E
