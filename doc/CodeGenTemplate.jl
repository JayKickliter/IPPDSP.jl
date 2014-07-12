for ( julia_fun, ippf_prefix )  in  [ (   :"insert julia function name",  "ippsFunctionBaseName"  ) ]
    for ( "TypeSignatures" ) in "AnArrayOfTuples"
                        
        
        julia_fun! = symbol(string(julia_fun, '!')) # in-place version if it makes sense
        
        @eval begin
            function $(julia_fun!)( buffer, "julia function arguments" )
                sigLen = length( signal )
                outLen = length( buffer )                
                @ippscall( $ippfsr,  (  "the c function argument types" ),
                                        "julia varables to pass to the c function")
                buffer                                                                                                               
            end # function         

            $(julia_fun)( "julia function arguments" ) = $(julia_fun!)( "create a buffer", "julia function arguments" )

        end # eval
    end # type loop
end # function loop