#############################################################################
##
##  PackageInfo.g        GAP4 Package `Guarana'               Bjoern Assmann
##  

SetPackageInfo( rec(

PackageName := "Guarana",
Subtitle := "Applications of Lie methods for computations with infinite polycyclic  groups",
Version := "0.94",
Date := "27/04/2012",

ArchiveURL := Concatenation([ 
"http://www-circa.mcs.st-andrews.ac.uk/~jjm/software/Guarana/Guarana-", 
~.Version]),
ArchiveFormats := ".tar.gz",


Persons := [

  rec(
      LastName      := "Assmann",
      FirstNames    := "Bjoern",
      IsAuthor      := true,
      IsMaintainer  := false,
      Email         := "bjoern@mcs.st-and.ac.uk",
      PostalAddress := Concatenation( [
            "Mathematical Institute\n",
            "University of St. Andrews\n",
            "North Haugh, St. Andrews\n Fife, KY 16 9SS, Scotland" ] ),
      Place         := "St. Andrews",
      Institution   := "University of St. Andrews"),

  rec(
      LastName      := "McDermott",
      FirstNames    := "John",
      IsAuthor      := false,
      IsMaintainer  := true,
      Email         := "jjm@mcs.st-and.ac.uk",
      PostalAddress := Concatenation( [
            "Mathematical Institute\n",
            "University of St. Andrews\n",
            "North Haugh, St. Andrews\n Fife, KY 16 9SS, Scotland" ] ),
      Place         := "St. Andrews",
      Institution   := "University of St. Andrews"),

],

Status := "deposited",
#CommunicatedBy := "Charles Wright (Eugene)",
#AcceptDate := "08/2005",

README_URL := "http://www-circa.mcs.st-andrews.ac.uk/~jjm/software/Guarana/README",
PackageInfoURL := "http://www-circa.mcs.st-andrews.ac.uk/~jjm/software/Guarana/PackageInfo.g",

AbstractHTML := "The Guarana package provides computational applications of the Mal'cev correspondence, in particular for collection in infinite polycyclic groups.", 

PackageWWWHome :="http://www-circa.mcs.st-andrews.ac.uk/~jjm/software/Guarana",

PackageDoc := rec(
  BookName  := "Guarana",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Applications of Lie methods for computations with infinite polycyclic groups",
  Autoload  := true ),

Dependencies := rec(
  GAP := ">= 4.4",
  NeededOtherPackages := [ [ "gapdoc",">=1.3"],
                           [ "polycyclic", ">=1.1" ], 
                           [ "polenta", ">=1.2.3" ],
                         # [ "radiroot", ">=2.0" ]
                         ],
  SuggestedOtherPackages := [ [ "nq", ">=2.0" ],
                              [ "alnuth", ">=3.0.0" ]
                            ], 
  ExternalConditions :=[], 
), 

AvailabilityTest := ReturnTrue,             
BannerString := Concatenation([ 
"Loading Guarana ",
~.Version,
" ... \n" ]),     
Autoload := false,
TestFile := "tst/testall.g",
Keywords := ["Mal'cev correspondence", "Collection", "Lie algebra", "Baker Campbell Haussdorff Formula", "polycyclic groups" ],    

));

#############################################################################
##
#E

