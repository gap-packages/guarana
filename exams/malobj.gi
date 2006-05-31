#############################################################################
##
#W mo_exs.gi                GUARANA package                     Bjoern Assmann
##
## Some Examples of Mal'cev records.  
##
#H  @(#)$Id$
##
#Y 2006
##
##
ExamplesOfSomeMalcevObjects := function( x )
    if x = 0 then 
        return GUARANA.Get_Unitriangular_MalcevObject( 1, 2);
    fi;
    if x in [1..5] then
        return GUARANA.Get_FNG_MalcevObject( 2, x );
    elif x in [6..10] then
        return GUARANA.Get_FNG_MalcevObject( 3, x-5 );
    elif x in [11..15] then 
        return GUARANA.Get_Unitriangular_MalcevObject(x-9 ,2 ); 
    else 
        Error( " " );
    fi;
end;

#############################################################################
##
#E 
