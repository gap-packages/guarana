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

DeclareProperty( "IsSymbolicElement", IsMalcevElement );

DeclareAttribute( "Weight", IsMalcevElement );


#############################################################################
##
#E
