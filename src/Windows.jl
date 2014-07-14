# Window functions

# without α
for ( julia_fun, ippf_prefix ) in   [   (   :bartlett,      "ippsWinBartlett"   ),  
                                        (   :blackman,      "ippsWinBlackmanStd"),
                                        (   :hamming,       "ippsWinHamming"    ),
                                        (   :hann,          "ippsWinHann"       ) ]

    julia_fun! = symbol(string(julia_fun, '!'))

    for ( T ) in [  :IPP16f,
                    :IPP32f, 
                    :IPP64f,
                    :IPP16fc,
                    :IPP32fc,
                    :IPP64fc  ]

        ippf  = string( ippf_prefix, IPPSuffix( T ) )

        @eval begin
            function $(julia_fun!)( y::Array{$T}, x::Array{$T}  )
                n = length( y )
                length(x) == n || throw( DimensionMismatch( "Inconsistent array lengths." ))
                if n > 0
                    @ippscall( $ippf, ( Ptr{$T},    Ptr{$T},    IPPInt  ), 
                                        x,          y,          n       )
                end
                
                return y
            end
            
            $(julia_fun!)( x::Array{$T} ) = $(julia_fun!)( x, x )       
            $(julia_fun)(  x::Array{$T} ) = $(julia_fun!)( similar(x), x )
            
            function $(julia_fun!)( y::Array{$T}, n::Integer )
                ny = length(y)
                ny == n || throw( DimensionMismatch( "Trying to fill a $ny element buffer with $n elements." ))
                x  = ones( $T, n )
                $(julia_fun!)( y, x )
            end
            
            $(julia_fun)( T::Type, n::Integer ) = $(julia_fun!)( Array( T, n ), n )
            $(julia_fun)( n::Integer ) = $(julia_fun!)( Array( IPP64f, n ), n )
        end
    end
end

# with α
for ( julia_fun, ippf_prefix ) in   [   (   :blackman,      "ippsWinBlackman"   ),                          
                                        (   :kaiser,        "ippsWinKaiser"     ) ]
                            

    julia_fun! = symbol(string(julia_fun, '!'))

    for ( T ) in [  :IPP16f, 
                    :IPP32f,  
                    :IPP64f, 
                    :IPP16fc, 
                    :IPP32fc, 
                    :IPP64fc  ]

        ippf  = string( ippf_prefix, IPPSuffix( T ) )

        @eval begin
            function $(julia_fun!)( y::Array{$T}, x::Array{$T}, α::FloatingPoint  )
                n = length( y )
                length(x) == n || throw( DimensionMismatch( "Inconsistent array lengths." ))
                if n > 0
                    @ippscall( $ippf, ( Ptr{$T},    Ptr{$T},    IPPInt, IPP64f  ), 
                                        x,          y,          n,      α       )
                end
                
                return y
            end
            
            $(julia_fun!)( x::Array{$T}, α::FloatingPoint ) = $(julia_fun!)( x, x, α )       
            $(julia_fun)(  x::Array{$T}, α::FloatingPoint ) = $(julia_fun!)( similar(x), x, α )
            
            function $(julia_fun!)( y::Array{$T}, n::Integer, α::FloatingPoint )
                ny = length(y)
                ny == n || throw( DimensionMismatch( "Trying to fill a $ny element buffer with $n elements." ))
                x  = ones( $T, n )
                $(julia_fun!)( y, x, α )
            end
            
            $(julia_fun)( T::Type, n::Integer, α::FloatingPoint ) = $(julia_fun!)( Array( T, n ), n, α )
            $(julia_fun)( n::Integer, α::FloatingPoint )          = $(julia_fun!)( Array( IPP64f, n ), n, α )
        end
    end
end
