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

############################################################################
#R  read other packages
##
RequirePackage( "polycyclic" );

#############################################################################
##
#E

