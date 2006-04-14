

#############################################################################
# Modification code of Werner
# for the computation of the constructive pc-sequence.
#
# Corresponding to the list of matrices <gens> we have now
# also a list of shadows <gens_s>.
# All computations which are done with <gens> are done with <gens_s> as
# well. 

##
##    Initialise a new level for the recursive sifting structure.
##
##    At each level we store matrices of the apprpriate weight and for each
##    matrix its inverse and the diagonal of the correct weight.
##
MAT_MakeNewLevel := function( w )
                                                                               
    InfoMatrixNq( "#I  MAT_MakeNewLevel( ", w, " ) called\n" );
    return rec( weight :=   w,
                matrices := [],
                matrices_s := [],
                inverses := [],
                inverses_s := [],
                diags :=    [] );
end;

MAT_SiftUpperUnitriMat := function() return 0; end;

##
##    When a matrix was added to a level of the sifting structure, then
##    commutators with this matrix and the generators of the group have to be
##    computed and sifted in order to keep the sifting structure from this
##    level on closed under taking commutators.
##
##    This function computes the necessary commutators and sifts each of
##    them.
##
MAT_FormCommutators := function( gens, gens_s, level, j )
    local   C,  Mj,  i, Mj_s, C_s;
                                                                               
    InfoMatrixNq( "#I  Forming commutators on level ", level.weight, "\n" );
                                                                               
    Mj := level.matrices[j];
    Mj_s := level.matrices_s[j];
    for i in [1..Length(gens)] do
        C := Comm( Mj,gens[i] );
        C_s := Comm( Mj_s, gens_s[i] );
        if not IsIdentityMat(C) then
            if not IsBound( level.nextlevel ) then
                level.nextlevel := MAT_MakeNewLevel( level.weight+1 );
            fi;
            MAT_SiftUpperUnitriMat( gens, gens_s,  level.nextlevel, C, C_s );
        fi;
    od;
    return;
end;

##
##    Sift the unitriangular matrix <M> through the recursive sifting
##    structure <level>.  It is assumed that the weight of <M> is equal to or
##    larger than the weight of <level>.
##
##    This function checks, if there is a matrix N at this level such that M
##    N^k for suitable k has higher weight than M.  If N does not exist, M is
##    added to <level> and all commutators of M with the generators of the
##    group are formed and sifted.  If there is such a matrix N, then M N^k
##    is sifted through the next level.
##
MAT_SiftUpperUnitriMat := function( gens, gens_s, level, M, M_s )
    local   w,  d,  h,  r,  R,R_s,  Ri,Ri_s,  c,  rr,  RR, RR_s;
                                                                               
    w := WeightUpperUnitriMat( M );
    if w > level.weight then
        if not IsBound( level.nextlevel ) then
            level.nextlevel := MAT_MakeNewLevel( level.weight+1 );
        fi;
        MAT_SiftUpperUnitriMat( gens, gens_s, level.nextlevel, M, M_s );
        return;
    fi;
    InfoMatrixNq( "#I  Sifting at level ", level.weight, " with " );
                                                                               
    d := UpperDiagonalOfMat( M, w+1 );
    h := 1; while h <= Length(d) and d[h] = 0 do h := h+1; od;
                                                                               
    while h <= Length(d) do
        if IsBound(level.diags[h]) then
            r  := level.diags[ h ];
            R  := level.matrices[ h ];
            R_s := level.matrices_s[h];
            Ri := level.inverses[ h ];
            Ri_s := level.inverses_s[h];
            c := Int( d[h] / r[h] );
            InfoMatrixNq( " ", c );
            if c <> 0 then
                d := d - c * r;
                if c > 0 then  M := Ri^c * M;  
                               M_s := Ri_s^c * M_s;
                else           M := R^(-c) * M;
                               M_s := R_s^(-c) * M_s;
                fi;
            fi;
            rr := r; r := d; d := rr;
            RR := R; R := M; M := RR;
            RR_s := R_s; R_s := M_s; M_s := RR_s;
            while r[h] <> 0 do
                c := Int( d[h] / r[h] );
            InfoMatrixNq( " ", c );
                if c <> 0 then
                    d := d - c  * r;
                    M := R^(-c) * M;
                    M_s := R_s^(-c) * M_s;
                fi;
                rr := r; r := d; d := rr;
                RR := R; R := M; M := RR;
                RR_s := R_s; R_s := M_s; M_s := RR_s;
            od;
            if d <> level.diags[ h ] then
                level.diags[ h ] := d;
                level.matrices[ h ] := M;
                level.matrices_s[ h ] := M_s;
                level.inverses[ h ] := M^-1;
                level.inverses_s[ h ] := M_s^-1;
                InfoMatrixNq( "\n" );
                MAT_FormCommutators( gens, gens_s, level, h );
 InfoMatrixNq( "#I  continuing reduction on level ", level.weight, " with " );
            fi;
            d := r;
            M := R;
            M_s := R_s;
        else
            level.matrices[ h ] := M;
            level.matrices_s[ h ] := M_s;
            level.inverses[ h ] := M^-1;
            level.inverses_s[ h ] := M_s^-1;
            level.diags[ h ]    := d;
            InfoMatrixNq( "\n" );
            MAT_FormCommutators( gens, gens_s, level, h );
           return;
        fi;
        while h <= Length(d) and d[h] = 0 do h := h+1; od;
    od;
    InfoMatrixNq( "\n" );
                                                                               
    if WeightUpperUnitriMat(M) < Length(M)-1 then
        if not IsBound( level.nextlevel ) then
            level.nextlevel := MAT_MakeNewLevel( level.weight+1 );
        fi;
        MAT_SiftUpperUnitriMat( gens, gens_s, level.nextlevel, M, M_s );
    fi;
end;
                                                                               
##
##    The subgroup U of GL(n,Z) of upper unitriangular matrices is a
##    nilpotent group.  The n-th term of the lower central series of U
##    consists of all unitriangular matrices of weight at least n-1.  This
##    defines a filtration on each subgroup of U.
##
##    This function computes this filtration for the unitriangular matrix
##    group  G.
##
MAT_SiftUpperUnitriMatGroup :=  function( gens, gens_s )
    local   firstlevel,  g, g_s, i;
                                                                               
    firstlevel := MAT_MakeNewLevel( 0 );
    for i in [1..Length( gens )] do
        g := gens[i];
        g_s := gens_s[i];
        MAT_SiftUpperUnitriMat( gens, gens_s, firstlevel, g, g_s );
    od;
    return firstlevel;
end;
                                                                               
##
##
##
MAT_PolycyclicGenerators :=  function( L )
    local   matrices, gens_s,  gens,  i,  l;
                                                                               
    matrices := Compacted( L.matrices );
    gens_s := Compacted( L.matrices_s );
    gens := [];
    for i in [1..Length(L.diags)] do
        if IsBound( L.diags[i] ) then Add( gens, [L.weight,i] ); fi;
    od;
                                                                               
    l := L;
    while IsBound( l.nextlevel ) do
        l := l.nextlevel;
                                                                               
        Append( matrices, Compacted( l.matrices ) );
        Append( gens_s, Compacted( l.matrices_s ) );
                                                                               
        for i in [1..Length(l.diags)] do
            if IsBound( l.diags[i] ) then Add( gens, [l.weight,i] ); fi;
        od;
                                                                               
    od;
                                                                               
    return rec( gens := gens, matrices := matrices, gens_s := gens_s );
end;

#
# end of modificated code of werner
############################################################################

