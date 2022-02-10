#############################################################################
##
##  PackageInfo.g        GAP4 Package `Guarana'               Bjoern Assmann
##  

SetPackageInfo( rec(

PackageName := "Guarana",
Subtitle := "Applications of Lie methods for computations with infinite polycyclic groups",
Version := "0.96.3",
Date := "11/02/2022", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [

  rec(
      LastName      := "Assmann",
      FirstNames    := "Björn",
      IsAuthor      := true,
      IsMaintainer  := false,
      Email         := "bjoern@mcs.st-and.ac.uk",
  ),

  rec(
      LastName      := "McDermott",
      FirstNames    := "John",
      IsAuthor      := false,
      IsMaintainer  := false,
      # Former maintainer; converted the manual to GAPDoc
      Email         := "jjm@mcs.st-and.ac.uk",
      PostalAddress := Concatenation( [
            "Mathematical Institute\n",
            "University of St. Andrews\n",
            "North Haugh, St. Andrews\n Fife, KY 16 9SS, Scotland" ] ),
      Place         := "St. Andrews",
      Institution   := "University of St. Andrews"),

  rec(
      LastName      := "GAP Team",
      FirstNames    := "The",
      IsAuthor      := false,
      IsMaintainer  := true,
      Email         := "support@gap-system.org",
  ),

],

Status := "deposited",
#CommunicatedBy := "Charles Wright (Eugene)",
#AcceptDate := "08/2005",

PackageWWWHome  := "https://gap-packages.github.io/guarana/",
README_URL      := Concatenation( ~.PackageWWWHome, "README.md" ),
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
SourceRepository := rec(
    Type := "git",
    URL := "https://github.com/gap-packages/guarana",
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/guarana-", ~.Version ),
ArchiveFormats := ".tar.gz",

AbstractHTML := "The Guarana package provides computational applications of the Mal'cev correspondence, in particular for collection in infinite polycyclic groups.", 

PackageDoc := rec(
  BookName  := "Guarana",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Applications of Lie methods for computations with infinite polycyclic groups",
  Autoload  := true ),

Dependencies := rec(
  GAP := ">= 4.7",
  NeededOtherPackages := [ [ "gapdoc",">=1.3"],
                           [ "polycyclic", ">=2.11" ],
                           [ "polenta", ">=1.2.3" ],
                         # [ "radiroot", ">=2.0" ]
                         ],
  SuggestedOtherPackages := [ [ "nq", ">=2.0" ],
                              [ "alnuth", ">=3.0.0" ]
                            ], 
  ExternalConditions :=[], 
), 

AvailabilityTest := ReturnTrue,             

TestFile := "tst/testall.g",
Keywords := ["Mal'cev correspondence", "Collection", "Lie algebra", "Baker Campbell Haussdorff Formula", "polycyclic groups" ],    

  AutoDoc := rec(
    TitlePage := rec(
      Copyright := "&copyright; 2007 Björn Assmann.",
    ),
  ),

));

#############################################################################
##
#E

