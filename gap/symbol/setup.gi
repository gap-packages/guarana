#############################################################################
##
#W setup.gi               GUARANA package                     Bjoern Assmann
##
## Methods for setting up the data structures that are needed for 
## symbolic collection ( in CN ).
##
#H  @(#)$Id$
##
#Y 2006
##
##

if false then 
#LoadPackage( "radiroot" );
#RereadPackage( "guarana", "gap/symbol/setup.gi" );

mc := ExamplesOfSomeMalcevCollectors( 1 );
GUARANA.SC_FullSetup( mc );

fi;

GUARANA.SC_ComputeJordanDecomp := function( malCol )
    local hl_CNN, semisimple_auts, unipotent_auts, mat, split, symCol, i, 
          c_s, c_u;

    # get HL of CN/N
    hl_CNN := malCol!.lengths[2];

    # compute jordan decomposition for all relevant matrices.
    semisimple_auts := [];
    unipotent_auts := [];
    for i in [1..hl_CNN] do 
        mat := malCol!.C_lieAuts[i];
        split := POL_MultiplicativeJordanDecomposition( mat );
        Add( semisimple_auts, split[1] );
        Add( unipotent_auts, split[2] );
    od;

    # check wheter the semisimple and unipotent automorphims commute
    for c_s in semisimple_auts do 
        for c_u in unipotent_auts do 
            if Comm( c_s, c_u ) <> c_s^0 then 
                Error( "semisimple and unipotent part does not commute." );
            fi;
        od;
    od;

    symCol := rec( semisimple_auts := semisimple_auts, 
                   unipotent_auts := unipotent_auts );

    malCol!.symCol := symCol;
end;

# compute the transformation matrix over some extension field F,
# that diagonalises the semisimple automorphisms.
GUARANA.SC_DiagonaliseSemisimple := function( malCol )
    local mats, res;

    # get semisimple automorphisms
    mats := malCol!.symCol.semisimple_auts;

    # diagonalize them 
    res := POL_DiagonalizeMatsSimultaneously( mats );

    # store the all needed information
    malCol!.symCol.semisimple_auts_diag := res.diagonalmats;
    malCol!.symCol.splitField := res.diag_rec.splitField;
    malCol!.symCol.T :=res.T; 
    malCol!.symCol.Tinv := res.Tinv;
    malCol!.symCol.res := res;

end;

# set up the variables that will represent the entries of the 
# from w_ij^y_i.
# Note:
# - each matrix should have its own set of variables,
#   because we want to distinguish between different y_i. 
# - the same EV should get the same variable.
# - EV 1 should not get anything.
## IN mat ................... diagonal matrix
##    i ..................... index in the pcs of CN.
##
## OUT
## A list of variables "w_ij" that corresponds to an entry w_ij^y_i
## Further a list with the corresponding eigenvalues.
## A diagonal matrix with the w_ij^k at the approbriate places.
##
GUARANA.SC_VariablesByDiagMat := function( mat, i, malCol )
    local deg, variables, eigenvalues, var_mat, F, ev, pos, string, w, j;

    deg := Length( mat );
    variables := [];
    eigenvalues := [];
    F := malCol!.symCol.splitField;
    var_mat := IdentityMat( deg, F );

    for j in [1..deg] do 
        ev := mat[j][j]; 
        if ev <> One( F ) then 
            pos := Position( eigenvalues, ev );
            if pos = fail then 
                Add( eigenvalues, ev );
                #create variable
                string := Concatenation( "w", String(i), String(j) );
                w := Indeterminate( F, string:new );
                Add( variables, w );
                #store it also in var_mat
                var_mat[j][j] := w;
            else 
                var_mat[j][j] := variables[pos];
            fi;
        fi;
    od;
    return rec( variables := variables, 
                eigenvalues := eigenvalues,
                mat := var_mat );
end;

GUARANA.SC_AddOmegaVariables := function( malCol )
    local hl_CN_N, var_recs, mat, var_rec, i;

    hl_CN_N := Length( malCol!.symCol.semisimple_auts_diag );
    var_recs := [];
    for i in [1..hl_CN_N] do 
        mat := malCol!.symCol.semisimple_auts_diag[i];
        var_rec := GUARANA.SC_VariablesByDiagMat( mat,i,malCol );
        Add( var_recs, var_rec );
    od;
    malCol!.symCol.omega_var_recs := var_recs;
end;

# compute and store the matrix that corresponds to conjugation
# with the semisimple part of phi(c_1)^y_1 * ... * phi(c_k)^y_k
# Don't forget to include the base change matrix T.
GUARANA.SC_SetFullSemisimpleAction := function( malCol )
    local hl_CN_N, dim_L, F, fullSSAction, omega_var_recs, T, Tinv, i;

    # setup
    hl_CN_N := Length( malCol!.symCol.omega_var_recs );
    dim_L := malCol!.lengths[3];
    F := malCol!.symCol.splitField;
    fullSSAction := IdentityMat( dim_L, F );
    omega_var_recs := malCol!.symCol.omega_var_recs;
    T := malCol!.symCol.T;
    Tinv := malCol!.symCol.Tinv;


    for i in [1..hl_CN_N] do
        fullSSAction := fullSSAction * omega_var_recs[i].mat;
    od;

    fullSSAction := Tinv * fullSSAction * T;

    malCol!.symCol.fullSemisimpleAction := fullSSAction;
end;

## IN malCol ................... Mal'cev collector 
##    left ..................... string that is the name for the variables
##                               of the left element, for example "x_"
##
##    right .................... string that is the name for the variables
##                               of the left element, for example "x_"
##
GUARANA.SC_SetCollectionVariableNames := function( malCol, left, right )
    local hl_CN, left_vars, right_vars, stringl, v, stringr, i,F;

    # setup
    hl_CN := Sum(  malCol!.lengths{[2,3]} );
    F := malCol!.symCol.splitField;
   
    left_vars := [];
    right_vars := [];
    for i in [1..hl_CN] do 
        stringl := Concatenation( left, String(i) );
        v := Indeterminate( F, stringl:new );
        Add( left_vars, v );
        stringr := Concatenation( right, String(i) );
        v := Indeterminate( F, stringr:new );
        Add( right_vars, v );
    od;
    malCol!.symCol.left_vars := left_vars;
    malCol!.symCol.right_vars := right_vars;
end;

## IN mat .................. unipotent matrix
##    var .................. variable
##
## OUT mat^var computed via exp( var* log(mat ) )
##
GUARANA.SC_SymbolicPowerUnipotentMat := function( mat, var )
    local log, res;

    log := POL_Logarithm( mat );
    res := POL_Exponential( var*log );
    return res;
end;

# compute and store the matrix that corresponds to conjugation
#with the unipotent part of phi(c_1)^y_1 * ... * phi(c_k)^y_k
# The variables y_1,...,y_k should have been stored previously 
# somewhere in the symbolic collector
#
GUARANA.SC_SetFullUnipotentAction := function( malCol )
    local dim_L, F, fullUAction, unipotent_auts, hl_CN_N, u_i, var_i, symP, i;

    # setup 
    dim_L := malCol!.lengths[3];
    F := malCol!.symCol.splitField;
    fullUAction := IdentityMat( dim_L, F );
    unipotent_auts := malCol!.symCol.unipotent_auts;
    hl_CN_N := Length( unipotent_auts );

    for i in [1..hl_CN_N] do
        u_i := unipotent_auts[i];
        u_i := u_i * One( F );
        var_i := malCol!.symCol.right_vars[i];
        symP := GUARANA.SC_SymbolicPowerUnipotentMat( u_i, var_i ); 
        fullUAction := fullUAction * symP;
    od;

    malCol!.symCol.fullUnipotentAction := fullUAction;
end;

GUARANA.SC_SetGenericElms := function( malCol )
    local left_vars, right_vars;

    if not IsBound( malCol!.symCol ) then 
        Error( "Symbolic collector is not set up." );
    fi;

    # get variables
    left_vars := malCol!.symCol.left_vars;
    right_vars := malCol!.symCol.right_vars;

    # store elms
    malCol!.symCol.left:=MalcevSymbolicCNElementByExponents(malCol,left_vars);
    malCol!.symCol.right:=MalcevSymbolicCNElementByExponents(malCol,right_vars);
    
end;

## IN
## malCol ..................... Malcev collector
##    n ....................... Malcev gen elment of N
##
## n^c given as Malcev gen elment of N where 
## c is the C-part of malCol!.symCol.right, 
## which is an internal generic element.
## Note that we use fact that c has only non-zero exponents in 
## the CN/N part.
##
GUARANA.SC_N_ConjugationByRightC_Elm := function( malCol, n )
    local coeffs;

    coeffs := Coefficients( n );

    # conjugate by semisimple action
    coeffs := coeffs * malCol!.symCol.fullSemisimpleAction;

    # conjugate by unipotent action
    coeffs := coeffs * malCol!.symCol.fullUnipotentAction;

    return MalcevSymbolicGenElementByCoefficients( malCol!.mo_NN, coeffs );
end;

## OUT
## Compute the collection function of g*h 
## where h is malCol!.symCol.right;
##
## Comment:
## g h = c(g)n(g) c(h)n(h)
##     = c(g)c(h0 n(g)^c(h) n(h)
##
GUARANA.SC_ComputeCollectionFuncByLeftElm := function( malCol, g )
    local h, c_new, n_c, n_new, res;

    # setup
    h := malCol!.symCol.right;
    GUARANA.COMP_OVER_EXT_FIELD := true;
    GUARANA.EXT_FIELD := malCol!.symCol.splitField;

    c_new := g!.c * h!.c;
    n_c := GUARANA.SC_N_ConjugationByRightC_Elm( malCol, g!.n );
    n_new := n_c * h!.n;
    res :=  MalcevCNElementBy2GenElements( malCol, c_new, n_new );
    Normalise( res );

    GUARANA.COMP_OVER_EXT_FIELD := false;
    return res;
end;

## computes the symbolic functions of the g*h where 
## g and h are the two generic elms stored in maCol!.symCol
##
GUARANA.SC_ComputeCollectionFunc := function( malCol )
    local left, res;

    left := malCol!.symCol.left;
    res := GUARANA.SC_ComputeCollectionFuncByLeftElm( malCol, left );
    malCol!.symCol.collFuncs := Exponents( res );
    return res;
end;

GUARANA.SC_FullSetup := function( malCol )
    GUARANA.SC_ComputeJordanDecomp( malCol );
    GUARANA.SC_DiagonaliseSemisimple( malCol );
    GUARANA.SC_AddOmegaVariables( malCol );
    GUARANA.SC_SetFullSemisimpleAction( malCol );
    GUARANA.SC_SetCollectionVariableNames( malCol, "x", "y" );
    GUARANA.SC_SetFullUnipotentAction( malCol );

    # use star polynomials for mulitplication in malcev collectors
    SetStarMethod( malCol!.mo_CC, "pols" );
    SetStarMethod( malCol!.mo_NN, "pols" );
    SetMultiplicationMethod( malCol!.mo_CC, GUARANA.MultMethodIsStar );
    SetMultiplicationMethod( malCol!.mo_NN, GUARANA.MultMethodIsStar );

    GUARANA.SC_SetGenericElms( malCol );
    GUARANA.SC_ComputeCollectionFunc( malCol );
end;

InstallGlobalFunction( AddSymbolicCollector,
function( malCol )
    GUARANA.SC_FullSetup( malCol );
end);

#############################################################################
##
#M  Value
##
## vals are indeterminants in an extension field of Q.
##                               
InstallOtherMethod(Value,"rat.fun., with one",
  true,[IsPolynomialFunction,IsList,IsList,IsRingElement, IsRingElement],0,
function(rf,inds,vals,one, oneExtField )
local i,fam,ivals,valextrep,d;

  if Length(inds)<>Length(vals) then
    Error("wrong number of values");
  fi;

  # convert indeterminates to numbers
  inds:= ShallowCopy( inds );
  for i in [1..Length(inds)] do
    if not IsPosInt(inds[i]) then
      inds[i]:=IndeterminateNumberOfUnivariateRationalFunction(inds[i]);
    fi;
  od;

  ivals:=[]; # values according to index

  fam:=CoefficientsFamily(FamilyObj(rf));


  valextrep:=function(f)
  local i,j,v,c,m,p;
    i:=1;
    v:=Zero(fam)*one;
    while i<=Length(f) do
      c:=f[i];
      m:=one;
      j:=1;
      while j<=Length(c) do
	if not IsBound(ivals[c[j]]) then
	  p:=Position(inds,c[j]);
	  if p<>fail then
	    ivals[c[j]]:=vals[p];
	  else
	    ivals[c[j]]:=UnivariatePolynomialByCoefficients(
		      fam,[Zero(fam),One(fam)],c[j]);
	  fi;
	fi;
	m:=m*(ivals[c[j]]*one)^c[j+1];
	j:=j+2;
      od;
      v:=v+(f[i+1]*oneExtField)*m;
      i:=i+2;
    od;
    return v;
  end;

  if HasIsPolynomial(rf) and IsPolynomial(rf) then
    return valextrep(ExtRepPolynomialRatFun(rf));
  else
    d:=valextrep(ExtRepDenominatorRatFun(rf));
    if IsZero(d) then
      Error("Denominator evaluates as zero");
    fi;
    return valextrep(ExtRepNumeratorRatFun(rf))/d;
  fi;

end );

# for runtimes for symbolic collection paper
GUARANA.SC_RuntimesFullSetup := function( exams, start, stop ) 
    local times, exam, rec1,rec2, time1,time2, mc,i,time;
    times := [];
    for i in [start..stop] do 
        exam := exams[i];
        # get malcev collector
        rec1 := GUARANA.CompleteRuntime2( MalcevCollectorConstruction, exam );
        mc := rec1.result;
        time1 := rec1.time;
        mc := MalcevCollectorConstruction( exam ); 
        # add symbolic information
        time2 := GUARANA.CompleteRuntime1( AddSymbolicCollector, mc );
        
        time := time1 + time2; 
        Print( "Index: ", i, " Time: ", time, "\n" );
        Add( times, time1 + time2 );
    od;
    return times;
end;


if false then 
 exams_F_2c := List( [1..6], x-> GUARANA.F_nc_Aut1( 2, x ) );
 GUARANA.SC_RuntimesFullSetup( exams_F_2c, 2, 6 );
  #[ 115, 264, 717, 5593, 45601 ]
  # n = 2,...,6
  
 exams_F_3c := List( [1..5], x-> GUARANA.F_nc_Aut2( 3, x ) );
 GUARANA.SC_RuntimesFullSetup( exams_F_3c, 2, 5 );
 #[ 293, 2388, 74402, 
 # n = 2,3,4

 exams_Tr_n_O1 := List( [1..8], x-> GUARANA.Tr_n_O1( x ) );
 GUARANA.SC_RuntimesFullSetup( exams_Tr_n_O1, 2, 8 );
 # [ 193, 360, 1759, 11228, 99285, 
 # n = 2,...,6

 exams_Tr_n_O2 := List( [1..6], x-> GUARANA.Tr_n_O2( x ) );
 GUARANA.SC_RuntimesFullSetup( exams_Tr_n_O2, 2, 6 );
 # [ 343, 2318, 112.763, 2.279.446,
 # n = 2,...,5


fi;
#############################################################################
##
#E
