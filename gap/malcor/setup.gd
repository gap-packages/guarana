#############################################################################
##
#W setup.gd               GUARANA package                     Bjoern Assmann
##
## Methods for symbolic log and exp. 
##
#H  @(#)$Id$
##
#Y 2006
##
##
DeclareCategory( "IsMalcevObject", IsObject );
DeclareRepresentation( "IsMalcevObjectRep", IsComponentObjectRep, 
                       [ "L", 
		         "dim",
			 "recTGroup",
			 "scTable",
			 "lieAlgebraType",
			 "weights",
			 "max_weight",
			 "malcevBasisInfo",
			 "log_method",
			 "exp_method",
			 "star_method",
             "mult_method",
			 "lie_fam",
			 "lie_elms_type" ] );

BindGlobal( "MalcevObjectFamily", 
     NewFamily( "MalcevObject", IsMalcevObjectRep ) );

DeclareGlobalFunction( "MalcevObjectConstruction"); 
DeclareGlobalFunction( "MalcevObjectByTGroup"); 

DeclareAttribute( "Exp", IsObject );

#############################################################################
##
#E
