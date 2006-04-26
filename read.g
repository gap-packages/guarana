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
ReadPackage( "guarana", "gap/help.gi" );
ReadPackage( "guarana", "data/bchdat.g" );
ReadPackage( "guarana", "gap/bch.gi" );
ReadPackage( "guarana", "exams/triang.gi" );
ReadPackage( "guarana", "exams/tgrps.gi" );
ReadPackage( "guarana", "gap/pcs.gi" );
ReadPackage( "guarana", "exams/sta.gi" );
ReadPackage( "guarana", "gap/setup.gi" );
ReadPackage( "guarana", "exams/recs.gi" );
ReadPackage( "guarana", "gap/symlog.gi" );
ReadPackage( "guarana", "gap/supple/almcom.gi" );
ReadPackage( "guarana", "gap/supple/almsup.gi" );
ReadPackage( "guarana", "exams/supple.gi" );

#############################################################################
##
#E
