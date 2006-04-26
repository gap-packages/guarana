#############################################################################
##
#W supple.gi              GUARANA package                     Bjoern Assmann
##
## Examples of polycyclic groups given together with a normal
## T-group N and a nilpotent almost supplement C. 
##
#H  @(#)$Id$
##
#Y 2006
##
#############################################################################
##
#F GUARANA.PolycyclicExmas( n )
## 
## IN 
## n ..................... integer between 1 and 13
##
## OUT
## Record containg 
## G ............ pc group
## 
GUARANA.PolycyclicExmas := function( n )
  local G, N, nat, G_img, G_img_fit, H_img, H, C, CN;

    G := ExamplesOfSomePcpGroups( n );
    N := FittingSubgroup( G );
    nat := GUARANA.Supple_NaturalHomomorphismNC( G, N );
    G_img := Image( nat );
    G_img_fit := FittingSubgroup( G_img );
    H_img := Centre( G_img_fit );
    H := PreImage( nat, H_img );
    C := GUARANA.NilpAlmSup( G, N, H );
    CN := ClosureGroup( N, C );

    return rec( G := G, H := H, CN := CN, C := C, N := N );
end;

#############################################################################
##
#E
