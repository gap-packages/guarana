#############################################################################
#W pcs.gi                  GUARANA package                     Bjoern Assmann
##
## Methods for the computation of a nice polycyclic sequence. 
## "Nice" means that we can do Mal'cev collection with respect to that
## polycyclic sequence. 
##
#H  @(#)$Id$
##
#Y 2006
##

############################################################################
##
#F GUARANA.TGroupRec( N )
## 
## IN
## N ................ a T-group given by a pcp
##
## OUT
## a record containing some information about N. 
## For example it contains a Mal'cev basis of N and 
## pcp with respect to that basis. 
## 
## Example
## gap> N := GUARANA.Examples_Unitriangular( 4, 2 );
## Pcp-group with orders [ 0, 0, 0, 0, 0, 0 ]
## gap> GUARANA.TGroupRec( N );
##  rec( N := Pcp-group with orders [ 0, 0, 0, 0, 0, 0 ], 
##    NN := Pcp-group with orders [ 0, 0, 0, 0, 0, 0 ], 
##   sers := [ Pcp-group with orders [ 0, 0, 0, 0, 0, 0 ], 
##    Pcp-group with orders [ 0, 0 ], Pcp-group with orders [  ] ], 
##    indices := [ [ 1 .. 4 ], [ 5, 6 ] ], class := 2, 
##   weights := [ 1, 1, 1, 1, 2, 2 ], largestAbelian := 2 )
##
## TODO
## - weights should be computed in a different way.
## - it should also be possible to start with a given Mal'cev basis
## - the current weights are not correct. 
## 
GUARANA.TGroupRec := function( N )
    local sers,i, NN,l,leFactors,s,indices,class,
    wei, weights,a, largestAbelian;
    sers := UpperCentralSeries( N );
    NN := PcpGroupBySeries( sers ); 
    #basisNN := GeneratorsOfGroup( NN );
    #iso := NN!.bijection;
    #basisN := GeneratorsOfGroup( Image( iso ) );

    # get indices of generators of the factors of the upper central series
    l := Length( sers );
    leFactors := [];
    for i in [1..l-1] do
        Add( leFactors, HirschLength( sers[i]) - 
	HirschLength( sers[i+1]));
    od;
    s := 1;
    indices := [];
    for i in [1..l-1] do
	Add( indices, [s..s+leFactors[i]-1] );
	s := s + leFactors[i];
    od;

    # nilpotency class of group
    class := Length( leFactors );

    # weights of generators
    wei := 1;
    weights := [];
    for a in leFactors do
	for i in [1..a] do
	    Add( weights, wei );
	od;
	wei := wei + 1;
    od;

    # get larges abelian group in the upper central series
    for i in [1..Length(sers)] do
	if IsAbelian( sers[i] ) then
	    largestAbelian := i;
	    break;
	fi;
    od;

    return rec( N := N, #basisN := basisN, 
    NN := NN,
	  sers := sers, indices := indices,
	   class := class, weights := weights,
	     largestAbelian := largestAbelian );
end;

