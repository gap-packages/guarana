#############################################################################
##
##  PackageInfo.g        GAP4 Package `Polenta'                Bjoern Assmann
##  

SetPackageInfo( rec(

PackageName := "Guarana",
Subtitle := "Applications of Lie methods for computations with infinite groups",
Version := "0.9",
Date := "08/04/2006",

ArchiveURL := Concatenation([ 
"http://www.icm.tu-bs.de/ag_algebra/software/assmann/Polenta/Polenta-", 
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

Status := "other",
#CommunicatedBy := "Charles Wright (Eugene)",
#AcceptDate := "08/2005",

README_URL := "http://www.icm.tu-bs.de/ag_algebra/software/assmann/Polenta/README",
PackageInfoURL := "http://www.icm.tu-bs.de/ag_algebra/software/assmann/Polenta/PackageInfo.g",

AbstractHTML := "The Guarana packages provives methods for fast collection in infinite polycyclic groups.", 

PackageWWWHome :="http://www.icm.tu-bs.de/ag_algebra/software/assmann/Polenta",

#PackageDoc := rec(          
#  BookName  := "Polenta",
#  ArchiveURLSubset := ["doc", "htm"],
#  HTMLStart := "htm/chapters.htm",
#  PDFFile   := "doc/manual.pdf",
#  SixFile   := "doc/manual.six",
#  LongTitle := "Polycyclic presentations for matrix groups",
#  Autoload  := true ),

Dependencies := rec(
  GAP := ">= 4.3fix4",
  NeededOtherPackages := [[ "polycyclic", ">=1.1" ]],
  SuggestedOtherPackages := [ ], 
  ExternalConditions :=[], 
), 

AvailabilityTest := ReturnTrue,             
BannerString := Concatenation([ 
"Loading Polenta ",
~.Version,
" ... \n" ]),     
Autoload := true,
TestFile := "tst/testall.g",
Keywords := ["Mal'cev correspondence", "Collection", "Lie algebra", "Baker Campbell Haussdorff Formula" ],    

));

#############################################################################
##
#E













