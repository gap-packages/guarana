#############################################################################
##
##  PackageInfo.g        GAP4 Package `Polenta'                Bjoern Assmann
##  

SetPackageInfo( rec(

PackageName := "Guarana",
Subtitle := "Applications of Lie methods for computations with infinite polycyclic  groups",
Version := "0.9",
Date := "08/04/2006",

ArchiveURL := Concatenation([ 
"http://www.cs.st-andrews.ac.uk/~bjoern/software/Guarana/Guarana-", 
~.Version]),
ArchiveFormats := ".tar.gz",


Persons := [

  rec(
      LastName      := "Assmann",
      FirstNames    := "Bjoern",
      IsAuthor      := true,
      IsMaintainer  := true,
      Email         := "bjoern@mcs.st-and.ac.uk",
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

README_URL := "http://www.cs.st-andrews.ac.uk/~bjoern/software/Guarana/README",
PackageInfoURL := "http://www.cs.st-andrews.ac.uk/~bjoern/software/Guarana/PackageInfo.g",

AbstractHTML := "The Guarana packages provives computational applications of the Mal'cev correspondence, in particular for collection in infinite polycyclic groups.", 

PackageWWWHome :="http://www.cs.st-andrews.ac.uk/~bjoern/software/Guarana",

PackageDoc := rec(          
  BookName  := "Guarana",
  ArchiveURLSubset := ["doc", "htm"],
  HTMLStart := "htm/chapters.htm",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Applications of Lie methods for computations with infinite polycyclic groups",
  Autoload  := true ),

Dependencies := rec(
  GAP := ">= 4.3fix4",
  NeededOtherPackages := [[ "polycyclic", ">=1.1" ], 
                          [ "polenta", ">=1.2.3" ],
                          [ "radiroot", ">=2.0" ]],
  SuggestedOtherPackages := [ ["nq", ">=2.0"]], 
  ExternalConditions :=[], 
), 

AvailabilityTest := ReturnTrue,             
BannerString := Concatenation([ 
"Loading Guarana ",
~.Version,
" ... \n" ]),     
Autoload := true,
TestFile := "tst/testall.g",
Keywords := ["Mal'cev correspondence", "Collection", "Lie algebra", "Baker Campbell Haussdorff Formula" ],    

));

#############################################################################
##
#E

