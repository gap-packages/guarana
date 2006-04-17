#############################################################################
##
#W    read.g                 Guarana Package                  Bjoern Assmannn
##
##    @(#)$Id$
##
##

#############################################################################
##
#V GUARANA is a holder for private function
##
BindGlobal("GUARANA", rec());

#############################################################################
##
#R  Read the install files.
##
ReadPackage( "guarana", "gap/bch.gi" );
ReadPackage( "guarana", "exams/triang.gi" );
ReadPackage( "guarana", "exams/tgrps.gi" );
ReadPackage( "guarana", "gap/pcs.gi" );
ReadPackage( "guarana", "exams/recs.gi" );

#############################################################################
##
#E
