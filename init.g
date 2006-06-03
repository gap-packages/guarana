#############################################################################
##
#W    init.g                 Guarana package                  Bjoern Assmann
##
##    @(#)$Id$
##

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
ReadPackage( GuaranaPkgName, "gap/collec/setup2.gd" );
ReadPackage( GuaranaPkgName, "gap/collec/elms.gd" );

############################################################################
#R  read other packages
##
RequirePackage( "polycyclic" );

#############################################################################
##
#E

