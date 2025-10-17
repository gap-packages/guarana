LoadPackage( "guarana" );
LoadPackage( "nq" );
LoadPackage( "alnuth" );
LoadPackage( "fga" ); # needed for GUARANA.GetSomeAutomorphsimOfF_nc to compute automorphisms of free groups

dirs := DirectoriesPackageLibrary( "guarana", "tst" );
TestDirectory(dirs, rec(exitGAP := true));
FORCE_QUIT_GAP(1);
