LoadPackage( "guarana" );
LoadPackage( "nq" );
LoadPackage( "alnuth" );

dirs := DirectoriesPackageLibrary( "guarana", "tst" );
TestDirectory(dirs, rec(exitGAP := true));
FORCE_QUIT_GAP(1);
