#############################################################################
##
#W setup.gi               GUARANA package                     Bjoern Assmann
##
## Methods for setting up the data structures that are needed for 
## Mal'cev collection.
##
#H  @(#)$Id$
##
#Y 2006
##
##

GUARANA.PcpGroupByPcs := function (G, pcs )
  local sers, l, H, i;
    # built corresponding series
    sers := [];
    l := Length( pcs );
    for i in [1..l+1] do
	Add( sers, Subgroup( G, pcs{[i..l]} ));
    od;

    # get pcp 
    H := PcpGroupBySeries( sers );

    return H;
end;
    
#############################################################################
##
#F GUARANA.SetupCollecRecord( args )
##
## IN
## G ................... polycyclic group whose polycyclic presentation
##                       is given with respect to a nice polycyclic 
##                       sequence, that is a sequence that fullfills
##                       all requirements that are needed for the 
##                       Malcev collection.
##         pcs = (f_1,...f_r,c_{r+1},...,c_{r+s},n_{r+s+1},...,n_{r+s+t})
##         (n_{r+s+1},...,n_{r+s+t}) is a Mal'cev basis of a normal
##                       T group N. 
##                       C =< c_{r+1},...,c_{r+s} > is a T-group and an 
##                       almost nilpotent supplement to N in G. 
##                       CN/N is free abelian.
##                       G/CN is finite.
## indezes ............. list containing three lists l1,l2,l3
##                       l1 = [1...r]
##                       l2 = [r+1,...,r+s]
##                       l3 = [r+s+1,...,r+s+t]
## N ................... normal T-group of N. 
## NN .................. pcp group isomorphic to N. 
##                       Note that in the gap sense NN is not a subgroup 
##                       of G. N and NN must have the same underlying 
##                       Mal'cev basis. 
## C ..................  T-group as described above. 
## CC .................  pcp group isomorphic to C. 
##                       The polycyclic presentation of CC must be 
##                       given with respect to Mal'cev basis 
##                       starting with (c_{r+s},...,c_{r+s}).
##                       CC!.bijection must be a bijection between 
##                       CC and C. 
## 
GUARANA.SetupCollecRecord := function( args )
  local G, indezes, N, NN, C, CC;
    G := args[1];
    indezes := args[2];
    N := args[3];
    NN := args[4];
    C := args[5];
    CC := args[6];
    return rec( G = G, indezes = indezes, N = N, NN = NN, 
                C = C, CC = CC );
end;
    
#############################################################################
##
#E
