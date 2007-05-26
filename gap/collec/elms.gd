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
              "n",
              "exps"] );

DeclareGlobalFunction( "MalcevCNElementByExponents" );
DeclareGlobalFunction( "MalcevCNElementByExponentsNC" );
DeclareGlobalFunction( "MalcevCNElementBy2Coefficients" );
DeclareGlobalFunction( "MalcevCNElementBy2Exponents" );
DeclareGlobalFunction( "MalcevCNElementBy2GenElements" );
DeclareGlobalFunction( "MalcevSymbolicCNElementByExponents" );

DeclareCategory( "IsMalcevGElement", IsObject );
DeclareCategoryFamily( "IsMalcevGElement" );
DeclareCategoryCollections( "IsMalcevGElement" );

DeclareRepresentation( "IsMalcevGElementRep", 
                        IsComponentObjectRep, 
            [ "malcevCollector", 
              "cn_elm",
              "exps_f",
              "exps"] );

DeclareGlobalFunction( "MalcevGElementByExponents" );
DeclareGlobalFunction( "MalcevGElementByExponentsNC" );
DeclareGlobalFunction( "MalcevGElementByCNElmAndExps" );

DeclareOperation( "NormalForm" , [IsObject] );
DeclareOperation( "Normalise" , [IsObject] );

DeclareOperation( "AttacheMalcevCollector", [IsPcpGroup,IsMalcevCollectorRep] );
DeclareAttribute( "AttachedMalcevCollector", IsPcpElement );
DeclareProperty( "IsMalcevPcpElement", IsPcpElement );
GUARANA.MALCEV_COLL_STORAGEPLACE := 1110;

#############################################################################
##
#E
