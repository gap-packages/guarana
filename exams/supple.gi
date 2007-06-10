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
#F GUARANA.PolycyclicExams( n )
## 
## IN 
## n ..................... integer between 1 and 13
##
## OUT
## Record containg 
## G ............ pc group
## K ............ subgroup of finite index of G that is polyZ
## N ............ Fitting subgroup of K
## H ............ subgroup of finite index of K such that H/N is nilpotent
## C ............ T-group that is a nilpotent almost supplement for N in K
## C_N .......... intersection of C and N
##
## n = 13 gives an interesting example to work with. 
## n = 11,12,14,15,16 give you a nilpotent group G
##                
## TODO
## This does not work for n=14, which is a bug in polycycyclic.
## 
##
GUARANA.PolycyclicExams := function( n )
  local G, K, N, nat, K_img, K_img_fit, H_img, H, C, CN, C_N;

    G := ExamplesOfSomePcpGroups( n );
    K := PolyZNormalSubgroup( G );
    N := FittingSubgroup( K );
    nat := GUARANA.Supple_NaturalHomomorphismNC( K, N );
    K_img := Image( nat );
    K_img_fit := FittingSubgroup( K_img );
    H_img := Centre( K_img_fit );
    H := PreImage( nat, H_img );
    C := GUARANA.Supple_NilpAlmSup( K, N, H );
    CN := ClosureGroup( N, C );
    C_N := Intersection( C, N );

    return rec( G := G, K := K, H := H, CN := CN, C := C, N := N, C_N := C_N );
end;

#############################################################################
##
#E
