#############################################################################
##
#W    init.g                 Guarana package                  Bjoern Assmann
##
##    @(#)$Id$
##

#############################################################################
##
#V GUARANA is a holder for private functions and variables
##
BindGlobal("GUARANA", rec());

#############################################################################
##
#R  Read the declaration files.
##
GuaranaPkgName := "guarana";

ReadPackage( GuaranaPkgName, "gap/malcor/setup.gd" );
ReadPackage( GuaranaPkgName, "gap/malcor/malelm.gd" );
ReadPackage( GuaranaPkgName, "gap/malcor/tstar.gd" );
ReadPackage( GuaranaPkgName, "gap/malcor/symlog.gd" );
ReadPackage( GuaranaPkgName, "gap/malcor/test.gd" );
ReadPackage( GuaranaPkgName, "gap/collec/setup.gd" );
ReadPackage( GuaranaPkgName, "gap/collec/elms.gd" );
ReadPackage( GuaranaPkgName, "gap/symbol/setup.gd" );

#############################################################################
##
#E

