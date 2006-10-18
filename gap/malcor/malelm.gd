#############################################################################
##
#W malelm.gd              GUA_package                     Bjoern Assmann
##
#H  @(#)$Id$
##
#Y 2006
##
##

DeclareCategory( "IsMalcevElement", IsObject );

DeclareCategory( "IsMalcevGenElement", IsMalcevElement );
DeclareCategoryFamily( "IsMalcevGenElement" );
DeclareCategoryCollections( "IsMalcevGenElement" );

DeclareRepresentation( "IsMalcevGenElementRep", 
                        IsComponentObjectRep, 
            [ "malcevObject", 
              "grp_elm",
              "lie_elm" ] );

DeclareGlobalFunction( "MalcevGenElementConstruction" );
DeclareGlobalFunction( "MalcevGenElementByExponents" );
DeclareGlobalFunction( "MalcevGenElementByCoefficients" );
DeclareGlobalFunction( "MalcevGenElementByLieElement" );
DeclareGlobalFunction( "MalcevGenElementByGrpElement" );
DeclareGlobalFunction( "MalcevSymbolicGenElementByExponents" );
DeclareGlobalFunction( "MalcevSymbolicGenElementByCoefficients" );

#############################################################################
##
## Lie elements
##
DeclareCategory( "IsMalcevLieElement", IsMalcevElement );
DeclareCategoryFamily( "IsMalcevLieElement" );
DeclareCategoryCollections( "IsMalcevLieElement" );

DeclareRepresentation( "IsMalcevLieElementRep",
                        IsComponentObjectRep,
			[ "malcevObject",
			  "coefficients",
			  "word",
			  "weight", 
			  "name" ] );

DeclareGlobalFunction( "MalcevLieElementConstruction" );
DeclareGlobalFunction( "MalcevLieElementByCoefficients" );
DeclareGlobalFunction( "MalcevLieElementByWord" );
DeclareGlobalFunction( "MalcevSymbolicLieElementByCoefficients" );

#############################################################################
##
## Group elements
##
DeclareCategory( "IsMalcevGrpElement", IsMalcevElement );
DeclareCategoryFamily( "IsMalcevGrpElement" );
DeclareCategoryCollections( "IsMalcevGrpElement" );

DeclareRepresentation( "IsMalcevGrpElementRep",
                        IsComponentObjectRep,
			[ "malcevObject",
			  "coefficients",
			  "weight", 
			  "name" ] );

DeclareGlobalFunction( "MalcevGrpElementConstruction" );
DeclareGlobalFunction( "MalcevGrpElementByExponents" );
DeclareGlobalFunction( "MalcevSymbolicGrpElementByExponents" );

DeclareProperty( "IsSymbolicElement", IsMalcevElement );

DeclareAttribute( "Weight", IsMalcevElement );

DeclareOperation( "LieElement", [IsMalcevGenElement] );
DeclareOperation( "GrpElement", [IsMalcevGenElement] );

#############################################################################
##
#E
