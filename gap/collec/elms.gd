#############################################################################
##
#W elms.gd               GUARANA package                     Bjoern Assmann
##
## Methods for setting up the data structures that are needed for 
## Mal'cev collection.
##
#H  @(#)$Id$
##
#Y 2006
##
##

DeclareCategory( "IsMalcevCNElement", IsObject );
DeclareCategoryFamily( "IsMalcevCNElement" );
DeclareCategoryCollections( "IsMalcevCNElement" );

DeclareRepresentation( "IsMalcevCNElementRep", 
                        IsComponentObjectRep, 
            [ "malcevCollector", 
              "c",
              "n" ] );

DeclareGlobalFunction( "MalcevCNElementByExponents" );
DeclareGlobalFunction( "MalcevCNElementBy2Coefficients" );
DeclareGlobalFunction( "MalcevCNElementBy2Exponents" );

#############################################################################
##
#E
