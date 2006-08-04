
DTPolynomialsByPcpGroup := function( G )
    local   coll;

    coll := Collector( G );

    if not IsBound( coll![PC_DEEP_THOUGHT_POLS] ) or 
       coll![PC_DEEP_THOUGHT_POLS] =[] then

        # Compute the deep thought polynomials
        coll![PC_DEEP_THOUGHT_POLS] := calcreps2(coll![PC_CONJUGATES], 8, 1);

        # Compute the orders of the genrators of dtrws
        CompleteOrdersOfRws(coll);

        # reduce the coefficients of the deep thought polynomials
        ReduceCoefficientsOfRws(coll);
    fi;

    return PolynomialsDTCollector( coll );

end;
