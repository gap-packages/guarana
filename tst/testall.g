LoadPackage( "guarana" );
dirs := DirectoriesPackageLibrary( "guarana", "tst" );
ReadTest( Filename( dirs, "guarana1.tst" ) );
ReadTest( Filename( dirs, "guarana2.tst" ) );

