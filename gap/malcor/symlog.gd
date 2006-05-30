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
DeclareOperation( "SetLogMethod", [IsObject, IsString] );
DeclareOperation( "SetExpMethod", [IsObject, IsString] );

DeclareOperation( "LogMethod" , [IsObject] );
DeclareOperation( "ExpMethod" , [IsObject] );

DeclareGlobalFunction( "AddLogAndExpPolynomials" );

#############################################################################
##
#E
