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
    local x;
    if x in [1..5] then
        return GUARANA.Get_FNG_MalcevObjects( 2, x );
    elif x in [6..10] then
        return GUARANA.Get_FNG_MalcevObjects( 3, x-5 );
    elif x in [11..15] then 
        return GUARANA.Get_Unitriangular_MalcevObjects(x-9 ,2 ); 
    else 
        Error( " " );
    fi;
fi;

#############################################################################
##
#E 
