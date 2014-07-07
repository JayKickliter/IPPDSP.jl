module WindowType
    const bartlett       = 0 
    const ippWinBlackman = 1
    const ippWinHamming  = 2
    const ippWinHann     = 3
    const ippWinRect     = 4
end

export  conv,
        xcorr,
        autocorr,
        autocorrb,
        autocorru,
        filt,
        FIRStateSize

FIRFilterTypes =    [   ( :Int16,               :Int16,                 "_16s"       ),
                        ( :Float32,             :Float32,               "_32f"       ),
                        ( :(Complex{Float32}),  :(Complex{Float32}),    "_32fc"      ),
                        ( :Float64,             :Float64,               "_64f"       ),
                        ( :(Complex{Float64}),  :(Complex{Float64}),    "_64fc"      ),
                        ( :Int32,               :Int16,                 "32s_16s"    ),
                        ( :(Complex{Int32}),    :(Complex{Int16}),      "32sc_16sc"  ),
                        ( :Float32,             :Int16,                 "32f_16s"    ),
                        ( :(Complex{Float32}),  :(Complex{Int16}),      "32fc_16sc"  ),
                        ( :Float64,             :Int16,                 "64f_16s"    ),
                        ( :Float64,             :Float32,               "64f_32f"    ),
                        ( :Float64,             :Int32,                 "64f_32s"    ),
                        ( :(Complex{Float64}),  :(Complex{Int16}),      "64fc_16sc"  ),
                        ( :(Complex{Float64}),  :(Complex{Int32}),      "64fc_32sc"  ),
                        ( :(Complex{Float64}),  :(Complex{Float32}),    "64fc_32fc"  )   ]




for ( julia_fun, ippf_prefix, types ) in    [   (   :conv,      "ippsConv", [   ( :Float32, "32f"     ),
                                                                                ( :Float64, "64f"     ),
                                                                                ( :Int16,   "16s_Sfs" ) ] ) ]
    julia_fun! = symbol(string(julia_fun, '!'))

    for ( T, ipp_suffix ) in types

        ippf  = string( ippf_prefix, '_', ipp_suffix )
        
        if eval(T)<:Integer # integer versions require a scaling factor
            @eval begin
                function $(julia_fun!)( y::Array{$T}, x1::Array{$T}, x2::Array{$T}, scale::Integer  )
                    ny  = length( y )
                    nx1 = length( x1 )
                    nx2 = length( x2 )
                    ny == nx1+nx2-1 || throw( DimensionMismatch( "length(y) must be length(x1)+length(x2)-1" ))
                    if ny > 3
                        @ippscall( $ippf, ( Ptr{$T},    IppInt,     Ptr{$T},    IppInt,     Ptr{$T},    IppInt  ),
                                            x1,         nx1,        x2,         nx2,        y,          scale   )
                    end
                    return y
                end
            
                $(julia_fun)(  x1::Array{$T}, x2::Array{$T}, scale::Integer ) = $(julia_fun!)( similar(x1, length(x1) + length(x2) - 1), x1, x2, scale )
            
            end
        else
            @eval begin
                function $(julia_fun!)( y::Array{$T}, x1::Array{$T}, x2::Array{$T}  )
                    ny  = length( y )
                    nx1 = length( x1 )
                    nx2 = length( x2 )
                    ny == nx1+nx2-1 || throw( DimensionMismatch( "length(y) must be length(x1)+length(x2)-1" ))
                    if ny > 3
                        @ippscall( $ippf, ( Ptr{$T},    IppInt,     Ptr{$T},    IppInt,     Ptr{$T} ),
                                            x1,         nx1,        x2,         nx2,        y       )
                    end
                    return y
                end
             
                $(julia_fun)(  x1::Array{$T}, x2::Array{$T} ) = $(julia_fun!)( similar(x1, length(x1) + length(x2) - 1), x1, x2 )
            end
        end
    end
end




for ( julia_fun, ippf_prefix, types ) in    [   (   :xcorr,      "ippsCrossCorr", [     ( :Float32,             "32f"       ),
                                                                                        ( :Float64,             "64f"       ),
                                                                                        ( :(Complex{Float32}),  "32fc"      ),
                                                                                        ( :(Complex{Float64}),  "64fc"      ),
                                                                                        ( :Int16,               "16s_Sfs"   ) ] ) ]

    julia_fun! = symbol(string(julia_fun, '!'))

    for ( T, ipp_suffix ) in types

        ippf  = string( ippf_prefix, '_', ipp_suffix )
        
        if eval(T)==Int16 # scaled version
            @eval begin
                function $(julia_fun!)( y::Array{$T}, x1::Array{$T}, x2::Array{$T}; scale::Integer = 0, lowlag::Integer = 0)
                    ny  = length( y )
                    nx1 = length( x1 )
                    nx2 = length( x2 )
                    ny > 1 && nx1 > 1 && nx2 > 1 || throw() # TODO: Check constraints                    
                    if ny > 3
                        @ippscall( $ippf, ( Ptr{$T},    IppInt,     Ptr{$T},    IppInt,     Ptr{$T},    IppInt, IppInt,  IppInt ),
                                            x1,         nx1,        x2,         nx2,        y,          ny,     lowlag,  scale  )
                    end
                    return y
                end
            
                $(julia_fun)(  x1::Array{$T}, x2::Array{$T}; args... ) = $(julia_fun!)( similar(x1, length(x1) + length(x2) - 1), x1, x2; args... )
            
            end
        else
            @eval begin
                function $(julia_fun!)( y::Array{$T}, x1::Array{$T}, x2::Array{$T}; lowlag::Integer = 0)
                    ny  = length( y )
                    nx1 = length( x1 )
                    nx2 = length( x2 )
                    ny > 1 && nx1 > 1 && nx2 > 1 || throw() # TODO: Check constraints                    
                    if ny > 3
                        @ippscall( $ippf, ( Ptr{$T},    IppInt,     Ptr{$T},    IppInt,     Ptr{$T},    IppInt, IppInt  ),
                                            x1,         nx1,        x2,         nx2,        y,          ny,     lowlag  )
                    end
                    return y
                end
            
                $(julia_fun)(  x1::Array{$T}, x2::Array{$T}; args... ) = $(julia_fun!)( similar(x1, length(x1) + length(x2) - 1), x1, x2; args... )
            end
        end
    end
end




for ( julia_fun, ippf_prefix )  in  [   (   :autocorr,  "ippsAutoCorr"          ),      
                                        (   :autocorrb, "ippsAutoCorr_NormA"    ),    
                                        (   :autocorru, "ippsAutoCorr_NormB"    )   ]
                                            
    julia_fun! = symbol(string(julia_fun, '!'))

    for ( T, ippf_suffix ) in [ (   :Int16,                 "16s_Sfs"   ),
                                (   :Float32,               "32f"       ), 
                                (   :Float64,               "64f"       ),
                                (   :(Complex{Float16}),    "16fc"      ),
                                (   :(Complex{Float32}),    "32fc"      ),
                                (   :(Complex{Float64}),    "64fc"      ) ]

        ippf  = string( ippf_prefix, '_', ippf_suffix )
        
        if eval(T) == Int16
            @eval begin
                function $(julia_fun!)( y::Array{$T}, x::Array{$T}; scale::Integer = 0 )
                    ny = length( y )
                    nx = length( x )
                    ny == nx || throw( DimensionMismatch("length(y) != length(x) ") )
                    if ny > 0
                        @ippscall( $ippf, ( Ptr{$T},    IppInt, Ptr{$T},    IppInt, IppInt  ), 
                                            x,          nx,     y,          ny,     scale   )
                    end
            
                    return y
                end
            
                $(julia_fun)(  x::Array{$T}; args... ) = $(julia_fun!)( similar(x, length(x)*2-1), x; args... )

            end
        else
            @eval begin
                function $(julia_fun!)( y::Array{$T}, x::Array{$T}  )
                    ny = length( y )
                    nx = length( x )
                    ny == nx|| throw( DimensionMismatch("length(y) != length(x) ") ) #todo: check this
                    if ny > 0
                        @ippscall( $ippf, ( Ptr{$T},    IppInt, Ptr{$T},    IppInt  ), 
                                            x,          nx,     y,          ny      )
                    end
            
                    return y
                end
                
                $(julia_fun)(  x::Array{$T} ) = $(julia_fun!)( similar(x), x )

            end            
        end
    end
end




for ( julia_fun, ippf_prefix1, ippf_prefix2 )  in  [ (   :FIRGetStateSize,  "ippsFIR",  "GetStateSize"  ) ]
    for ( T1, T2, ippf_suffix ) in FIRFilterTypes
    
        ippfsr  = string( ippf_prefix1,         ippf_prefix2, ippf_suffix )
        ippfmr  = string( ippf_prefix1, "MR",   ippf_prefix2, ippf_suffix )
        
        @eval begin
            function $(julia_fun)( ::Type{$T1}, ::Type{$T2}, tapslen::Integer  )
                buffersize = IppInt[0]  # the function will fill this value
                tapslen > 0 || throw( ) # todo: check this
                @ippscall( $ippfsr, ( IppInt,     Ptr{IppInt} ),
                                      tapslen,    buffersize    )

                return buffersize[1]
            end
        end
        
        if ippf_suffix != "_16s"
            @eval begin
                function $(julia_fun)( ::Type{$T1}, ::Type{$T2}, tapslen::Integer, upfactor::Integer, downfactor::Integer )
                    buffersize  = IppInt[0]  # the function will fill this value
                    tapslen > 0 || throw( ) # todo: check this
                    @ippscall( $ippfmr, (   IppInt,     IppInt,     IppInt,         Ptr{IppInt} ),
                                            tapslen,    upfactor,   downfactor,     buffersize  )
                    return buffersize[1]
                end                
            end
        end
    end
end

type FIRFilterSR{Tt, Ts}
    taps::Array{Tt, 1}
    delay_line::Array{Ts, 1}
    state::Array{Ptr{Uint8}, 1}
    buffer::Array{Uint8}
end

function FIRFilterSR{Tt, Ts}( ::Type{Ts}, taps::Array{Tt, 1} )
    n_taps      = length( taps )
    delay_line  = zeros( Ts, n_taps )
    state       = Array(Ptr{Uint8}, 1)
    buffer_size = FIRGetStateSize( Tt, Ts, n_taps )
    buffer      = zeros( Uint8, buffer_size )
    FIRFilterSR( taps, delay_line, state, buffer )
end



type FIRFilterMR{Tt, Ts}
    taps::Array{Tt, 1}
    delay_line::Array{Ts, 1}
    state::Array{Ptr{Uint8}, 1}
    buffer::Array{Uint8}
    interpolation::Integer
    decimation::Integer
end



#  Documentation is horrible, see this link: https://software.intel.com/en-us/forums/topic/307712
#
#  Name:         ippsFIRGetStateSize, ippsFIRMRGetStateSize,
#                ippsFIRInit, ippsFIRMRInit
#  Purpose:      ippsFIRGetStateSize - calculates the size of the FIR State
#                                                                   structure;
#                ippsFIRInit - initialize FIR state - set taps and delay line
#                using external memory buffer;
#  Parameters:
#      pTaps       - pointer to the filter coefficients;
#      tapsLen     - number of coefficients;
#      pDlyLine    - pointer to the delay line values, can be NULL;
#      ppState     - pointer to the FIR state created or NULL;
#      upFactor    - multi-rate up factor;
#      upPhase     - multi-rate up phase;
#      downFactor  - multi-rate down factor;
#      downPhase   - multi-rate down phase;
#      pStateSize  - pointer where to store the calculated FIR State structure
#                                                             size (in bytes);
#   Return:
#      status      - status value returned, its value are
#         ippStsNullPtrErr       - pointer(s) to the data is NULL
#         ippStsFIRLenErr        - tapsLen <= 0
#         ippStsFIRMRFactorErr   - factor <= 0
#         ippStsFIRMRPhaseErr    - phase < 0 || factor <= phase
#         ippStsNoErr            - otherwise

for ( julia_fun, ippf_prefix1, ippf_prefix2 )  in  [ (   :FIRInit,  "ippsFIR",  "Init"  ) ]
    for ( Tt, Ts, ippf_suffix ) in FIRFilterTypes
    
        ippfsr  = string( ippf_prefix1,         ippf_prefix2, ippf_suffix )
        ippfmr  = string( ippf_prefix1, "MR",   ippf_prefix2, ippf_suffix )
        
        @eval begin
            function $(julia_fun)( taps::Array{$Tt, 1}, delay_line::Array{$Ts, 1}, state::Array{Ptr{Uint8},1}, buffer::Array{Uint8, 1} )
                ntaps      = length( taps )
                ndelayline = length( delay_line ) 
                @ippscall( $ippfsr, (   Ptr{Uint8},     Ptr{$Tt},   IppInt, Ptr{$Ts},   Ptr{Uint8}  ),
                                        pointer(state), taps,       ntaps,  delay_line, buffer      )         
                nothing                                                                                                               
            end
        end
    end
end


FIRInit( fir::FIRFilterSR ) = FIRInit( fir.taps, fir.delay_line, fir.state, fir.buffer )

#  Names:         ippsFIR
#  Purpose:       FIR filter. Vector filtering
#  Parameters:
#      pSrcDst     - pointer to the input/output vector in in-place operation
#      pSrc        - pointer to the input vector
#      pDst        - pointer to the output vector
#      numIters    - number iterations (for single-rate equal length data vector)
#      pState      - pointer to the filter state
#      scaleFactor - scale factor value
#  Return:
#      ippStsContextMatchErr  - wrong state identifier
#      ippStsNullPtrErr       - pointer(s) to the data is NULL
#      ippStsSizeErr          - numIters is less or equal zero
#      ippStsNoErr            - otherwise
#  Note: for Multi-Rate filtering
#          length pSrc = numIters*downFactor
#          length pDst = numIters*upFactor
#          for inplace functions max this values
#

for ( julia_fun, ippf_prefix )  in  [ (   :filt,  "ippsFIR"  ) ]
    for ( Tt, Ts, ippf_suffix ) in FIRFilterTypes
    
        ippfsr  = string( ippf_prefix, ippf_suffix )
        
        @eval begin
            function $(julia_fun)( fltr::FIRFilterSR{$Tt, $Ts}, signal::Array{$Ts, 1} )
                out_vector = similar( signal )
                sig_len = length( signal )
                out_len = length( out_vector )
                state_ptr = fltr.state[1]
                @ippscall( $ippfsr,  (  Ptr{$Ts},   Ptr{$Ts},   IppInt,     Ptr{Uint8}  ),
                                        signal,     out_vector, sig_len,    state_ptr   )         
                out_vector                                                                                                               
            end
        end
    end
end