LoadPackage( "guarana" );
dirs := DirectoriesPackageLibrary( "guarana", "tst" );
if IsPackageMarkedForLoading("nq", "2.0") <> true then
 Print( "pkg nq not available - guarana1.tst not executed.\n" );
else
 ReadTest( Filename( dirs, "guarana1.tst" ) );
fi;
if IsPackageMarkedForLoading("alnuth", "3.0.0") <> true then
 Print( "pkg alnuth not available - guarana2.tst not executed.\n" );
else
 ReadTest( Filename( dirs, "guarana2.tst" ) );
fi;
