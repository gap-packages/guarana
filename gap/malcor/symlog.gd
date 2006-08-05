#############################################################################
##
#W symlog.gd              GUARANA package                     Bjoern Assmann
##
## Methods for symbolic log and exp. 
##
#H  @(#)$Id$
##
#Y 2006
##
##
DeclareOperation( "SetLogAndExpMethod", [IsObject, IsString] );
DeclareOperation( "SetStarMethod", [IsObject, IsString] );
DeclareOperation( "SetMultiplicationMethod", [IsObject, IsString] );


DeclareOperation( "LogMethod" , [IsObject] );
DeclareOperation( "ExpMethod" , [IsObject] );
DeclareOperation( "StarMethod" , [IsObject] );
DeclareOperation( "MultiplicationMethod" , [IsObject] );

DeclareGlobalFunction( "AddLogAndExpPolynomials" );
DeclareGlobalFunction( "AddStarPolynomials" );
DeclareGlobalFunction( "AddDTPolynomials" );
#############################################################################
##
#E
