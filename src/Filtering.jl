module WindowType
    const bartleTt       = 0 
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
        FIRFilter,
        filt,
        FIRInit


        

# These are common to most of the ipps FIR Filter functions 
# I haven't included functions with 
FIRFilterTypes =    [   ( :Float32,             :Float32,               "_32f"       ),
                        ( :(Complex{Float32}),  :(Complex{Float32}),    "_32fc"      ),
                        ( :Float64,             :Float64,               "_64f"       ),
                        ( :(Complex{Float64}),  :(Complex{Float64}),    "_64fc"      ),
                        ( :Float32,             :Int16,                 "32f_16s"    ),
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




################################################################################
#               ____ ____ ____ ____ ____ _    ____ ___ _ ____ _  _             #
#               |    |  | |__/ |__/ |___ |    |__|  |  | |  | |\ |             #
#               |___ |__| |  \ |  \ |___ |___ |  |  |  | |__| | \|             #
################################################################################                                           

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
                    ny > 1 && nx1 > 1 && nx2 > 1 || throw() # TODO: Check constrainTx                    
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
                    ny > 1 && nx1 > 1 && nx2 > 1 || throw() # TODO: Check constrainTx                    
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




################################################################################
#     ____ _ ____    ____ _ _    ___ ____ ____    ____ ___ ____ ___ ____       #
#     |___ | |__/    |___ | |     |  |___ |__/    [__   |  |__|  |  |___       #
#     |    | |  \    |    | |___  |  |___ |  \    ___]  |  |  |  |  |___       #
################################################################################  
#   Name:         ippsFIRRequiredStateSize, ippsFIRMRRequiredStateSize,
#                 ippsFIRInit, ippsFIRMRInit
#   Purpose:      ippsFIRRequiredStateSize - calculates the size of the FIR State
#                                                                    structure;
#                 ippsFIRInit - initialize FIR state - set taps and delay line
#                 using external memory buffer;
#   Parameters:
#       pTaps       - pointer to the filter coefficients;
#       tapsLen     - number of coefficients;
#       pDlyLine    - pointer to the delay line values, can be NULL;
#       ppState     - pointer to the FIR state created or NULL;
#       upFactor    - multi-rate up factor;
#       upPhase     - multi-rate up phase;
#       downFactor  - multi-rate down factor;
#       downPhase   - multi-rate down phase;
#       pStateSize  - pointer where to store the calculated FIR State structure
#                                                              size (in bytes);
#    Return:
#       status      - status value returned, its value are
#          ippStsNullPtrErr       - pointer(s) to the data is NULL
#          ippStsFIRLenErr        - tapsLen <= 0
#          ippStsFIRMRFactorErr   - factor <= 0
#          ippStsFIRMRPhaseErr    - phase < 0 || factor <= phase
#          ippStsNoErr            - otherwise

# IppStatus ippsFIRMRRequiredStateSize_64fc(int tapsLen, int upFactor, int downFactor, int* pBufferSize);

for ( julia_fun, ippf_prefix1, ippf_prefix2 )  in  [ (   :FIRRequiredStateSize,  "ippsFIR",  "GetStateSize"  ) ]
    for ( T_taps, T_sig, ippf_suffix ) in FIRFilterTypes
            
        # Single rate version
        ippfsr = string( ippf_prefix1, ippf_prefix2, ippf_suffix )
        @eval begin
            function $(julia_fun)( ::Type{$T_taps}, ::Type{$T_sig}, tapslen::Integer  )
                buffersize = IppInt[0]  # the function will fill this value
                tapslen > 0 || throw( ) # todo: check this
                @ippscall( $ippfsr, ( IppInt,     Ptr{IppInt} ),
                                      tapslen,    buffersize    )

                return buffersize[1]
            end
        end
        
        # Multirate version
        ippfmr = string( ippf_prefix1, "MR", ippf_prefix2, ippf_suffix )
        @eval begin
            function $(julia_fun)( ::Type{$T_taps}, ::Type{$T_sig}, tapslen::Integer, upFactor::Integer, downFactor::Integer )
                buffersize = IppInt[0]  # the function will fill this value
                tapslen > 0 || throw( ) # todo: check this
                @ippscall( $ippfmr, ( IppInt,   IppInt,    IppInt,      Ptr{IppInt} ),
                                      tapslen,  upFactor,  downFactor,  buffersize    )

                return buffersize[1]
            end
            $(julia_fun)( ::Type{$T_taps}, ::Type{$T_sig}, tapslen::Integer, resamp_ratio::Rational ) = $(julia_fun)( $T_taps, $T_sig, tapslen, num(resamp_ratio), den(resamp_ratio) )
        end        
    end
end



################################################################################
#      ____ _ ____    ____ _ _    ___ ____ ____    ___ _   _ ___  ____         #
#      |___ | |__/    |___ | |     |  |___ |__/     |   \_/  |__] |___         #
#      |    | |  \    |    | |___  |  |___ |  \     |    |   |    |___         #
################################################################################  


# IPP assings a pointer to state
# This points to an invisuble to us struct in buffer
# Intel's documentation is weak in this area, but we allocate the buffer, they
# instantsiat a stuct inside that buffer, and use state to point the location of
# that state

# Both the singe-rate and multirate states are the same, but are in different
# types to help with dispatching

# single rate state
type FIRSRState
    pointer::Vector{Ptr{Uint8}} 
    buffer::Vector{Uint8} 
end

# multirate state
type FIRMRState
    pointer::Vector{Ptr{Uint8}} 
    buffer::Vector{Uint8} 
end

type FIRFilter{Tt, Tx, Ts}
    taps::Vector{Tt}
    delayLine::Vector{Tx}
    state::Ts
    upFactor::Int
    downFactor::Int
    upPhase::Int
    downPhase::Int    
end

# Single rate constructor
function FIRFilter{Tx, Tt}( ::Type{Tx}, taps::Vector{Tt} )
    tapsLen   = length( taps )
    pointer   = Array( Ptr{Uint8}, 1)
    bufLen    = FIRRequiredStateSize( Tt, Tx, tapsLen )
    delayLine = zeros( Tx, tapsLen )
    buffer    = zeros( Uint8, bufLen )                        
    state     = FIRSRState( pointer, buffer )
    FIRInit( taps, delayLine, pointer, buffer )        
    FIRFilter( taps, delayLine, state, 1, 1, 0, 0 )
end    

# Multirate cunstructor
function FIRFilter{ Tx, Tt }( ::Type{Tx}, taps::Array{Tt, 1}, upFactor, downFactor, upPhase = 0, downPhase = 0 )
    tapsLen   = length( taps )
    pointer   = Array( Ptr{Uint8}, 1)
    bufLen    = FIRRequiredStateSize( Tt, Tx, tapsLen )
    delayLine = zeros( Tx, tapsLen )
    buffer    = zeros( Uint8, bufLen )                        
    state     = FIRMRState( pointer, buffer )
    FIRInit( taps, delayLine, pointer, buffer, upFactor, downFactor, upPhase, downPhase )        
    FIRFilter( taps, delayLine, state, upFactor, downFactor, upPhase, downPhase )        
end  



################################################################################
#           ____ _ ____    ____ _ _    ___ ____ ____    _ _  _ _ ___           #
#           |___ | |__/    |___ | |     |  |___ |__/    | |\ | |  |            #
#           |    | |  \    |    | |___  |  |___ |  \    | | \| |  |            #
################################################################################  
#  Documentation is horrible, see this link:
#                hTtps://software.intel.com/en-us/forums/topic/307712
#
#  Name:         ippsFIRRequiredStateSize, ippsFIRMRGeTxtateSize,
#                ippsFIRInit, ippsFIRMRInit
#  Purpose:      ippsFIRRequiredStateSize - calculates the size of the FIR State
#                                                                   structure;
#                ippsFIRInit - initialize FIR state - set taps and delay line
#                using external memory buffer;
#  Parameters:
#      pTaps       - pointer to the filter coefficienTx;
#      tapsLen     - number of coefficienTx;
#      pDlyLine    - pointer to the delay line values, can be NULL;
#      ppState     - pointer to the FIR state created or NULL;
#      upFactor    - multi-rate up factor;
#      upPhase     - multi-rate up phase;
#      downFactor  - multi-rate down factor;
#      downPhase   - multi-rate down phase;
#      pStateSize  - pointer where to store the calculated FIR State structure
#                                                             size (in bytes);
#   Return:
#      status      - status value returned, iTx value are
#         ippSTxNullPtrErr       - pointer(s) to the data is NULL
#         ippSTxFIRLenErr        - tapsLen <= 0
#         ippSTxFIRMRFactorErr   - factor <= 0
#         ippSTxFIRMRPhaseErr    - phase < 0 || factor <= phase
#         ippSTxNoErr            - otherwise
# 
for ( julia_fun, ippf_prefix1, ippf_prefix2 )  in  [ (   :FIRInit,  "ippsFIR",  "Init"  ) ], ( Tt, Tx, ippf_suffix ) in FIRFilterTypes

    # Single rate version
    ippfsr  = string( ippf_prefix1, ippf_prefix2, ippf_suffix )        
    @eval begin
        
        function $(julia_fun)(  taps::Array{$Tt, 1},
                                delayLine::Array{$Tx, 1},
                                state::Array{Ptr{Uint8},1},
                                buffer::Array{Uint8, 1}
                             )
            ntaps      = length( taps )
            ndelayline = length( delayLine ) 
            @ippscall( $ippfsr, (   Ptr{Uint8},     Ptr{$Tt},    IppInt,     Ptr{$Tx},      Ptr{Uint8}  ),
                                    pointer(state), taps,           ntaps,      delayLine,         buffer      )         
            nothing                                                                                                               
        end
    
    end
                
    # Multirate version
    ippfmr  = string( ippf_prefix1, "MR", ippf_prefix2, ippf_suffix )
    @eval begin
    
        function $(julia_fun)(  taps::Vector{ $Tt }, delayLine::Vector{ $Tx }, state::Vector{ Ptr{Uint8} }, buffer::Vector{ Uint8 }, upFactor, downFactor, upPhase, downPhase )
            ntaps      = length( taps )
            ndelayline = length( delayLine )
            0 <= upPhase < upFactor && 0 <= downPhase < downFactor || throw()
            @ippscall( $ippfmr, (   Ptr{Uint8},     Ptr{$Tt}, IppInt, IppInt,    IppInt,  IppInt,     IppInt,    Ptr{$Tx},  Ptr{Uint8} ),
                                    pointer(state), taps,     ntaps,  upFactor,  upPhase, downFactor, downPhase, delayLine, buffer     )         
            nothing                                                                                                               
        end
    
    end                    
end

FIRInit( fir::FIRFilter ) = FIRInit( fir.taps, fir.delayLine, fir.state, fir.buffer )





################################################################################
#      ____ _ ____    ____ _ _    ___ ____ ____    ____ _  _ ____ ____         #
#      |___ | |__/    |___ | |     |  |___ |__/    |___  \/  |___ |            #
#      |    | |  \    |    | |___  |  |___ |  \    |___ _/\_ |___ |___         #
################################################################################
# 
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
#      ippSTxContextMatchErr  - wrong state identifier
#      ippSTxNullPtrErr       - pointer(s) to the data is NULL
#      ippSTxSizeErr          - numIters is less or equal zero
#      ippSTxNoErr            - otherwise
#  Note: for Multi-Rate filtering
#          length pSrc = numIters*downFactor
#          length pDst = numIters*upFactor
#               for inplace functions max this values
# 
for ( julia_fun, ippf_prefix )  in  [ (   :filt,  "ippsFIR"  ) ]
    for ( Tt, Tx, ippf_suffix ) in FIRFilterTypes
    
        ippfsr     = string( ippf_prefix, ippf_suffix )
        julia_fun! = symbol(string(julia_fun, '!'))
        
        @eval begin
            function $(julia_fun!)( self::FIRFilter{$Tt, $Tx, FIRSRState},
                                    buffer::Array{$Tx, 1},
                                    signal::Array{$Tx, 1}
                                  )
                sigLen = length( signal )
                outLen = length( buffer )
                statePtr = pointer(self.state.pointer)
                @ippscall( $ippfsr,  (  Ptr{$Tx},  Ptr{$Tx},  IppInt,     Ptr{Uint8}  ),
                                        signal,         buffer,         sigLen,    statePtr   )         
                buffer                                                                                                               
            end            
            $(julia_fun)( self::FIRFilter{$Tt, $Tx, FIRSRState}, signal::Array{$Tx, 1} ) = $(julia_fun!)( self, similar( signal ), signal )
        end
        
        ippfmr = string( ippf_prefix, "MR", ippf_suffix )
        @eval begin            
            function $(julia_fun!)( self::FIRFilter{ $Tt, $Tx, FIRMRState },
                                    buffer::Array{$Tx, 1},
                                    signal::Array{$Tx, 1} )              
                sigLen = length( signal )
                outLen = length( buffer )
                sigLen * self.upFactor == outLen * self.downFactor || throw()
                iterations  = int( sigLen/self.downFactor )
                @ippscall( $ippfsr,  (  Ptr{$Tx},  Ptr{$Tx},  IppInt,     Ptr{Uint8}                  ),
                                        signal,    buffer,    iterations, pointer(self.state.pointer) )         
                buffer                                                                                                               
            end
            
            function $(julia_fun)( self::FIRFilter{$Tt, $Tx, FIRMRState}, signal::Vector{ $Tx } )
                sigLen = length( signal )
                sigLen % self.downFactor == 0 || throw()
                bufLen = int( sigLen * self.upFactor / self.downFactor )
                buffer = similar( signal, bufLen ) 
                $(julia_fun!)( self, buffer, signal )
            end    
        end
        
    end
end


#  Names:      ippsFIRGenLowpass_64f, ippsFIRGenHighpass_64f, ippsFIRGenBandpass_64f
#              ippsFIRGenBandstop_64f
#
#  Purpose:    This function computes the lowpass FIR filter coefficients
#              by windowing of ideal (infinite) filter coefficients segment
#
#  Parameters:
#      rfreq             cut off frequency (0 < rfreq < 0.5)
#
#      taps              pointer to the array which specifies
#                        the filter coefficients;
#
#      tapsLen           the number of taps in taps[] array (tapsLen>=5);
#
#      winType           the ippWindowType switch variable,
#                        which specifies the smoothing window type;
#
#      doNormal          if doNormal=0 the functions calculates
#                        non-normalized sequence of filter coefficients,
#                        in other cases the sequence of coefficients
#                        will be normalized.
#  Return:
#   ippStsNullPtrErr     the null pointer to taps[] array pass to function
#   ippStsSizeErr        the length of coefficient's array is less than five
#   ippStsSizeErr        the low or high frequency isn't satisfy
#                                    the condition 0 < rLowFreq < 0.5
#   ippStsNoErr          otherwise
#
#
#  IPPAPI(IppStatus, ippsFIRGenLowpass_64f, (Ipp64f rfreq, Ipp64f* taps, int tapsLen,
#                                              IppWinType winType, IppBool doNormal))
#  
#  IPPAPI(IppStatus, ippsFIRGenHighpass_64f, (Ipp64f rfreq, Ipp64f* taps, int tapsLen,
#                                               IppWinType winType, IppBool doNormal))
#  
#  IPPAPI(IppStatus, ippsFIRGenBandpass_64f, (Ipp64f rLowFreq, Ipp64f rHighFreq, Ipp64f* taps,
#                                       int tapsLen, IppWinType winType, IppBool doNormal))
#  
#  IPPAPI(IppStatus, ippsFIRGenBandstop_64f, (Ipp64f rLowFreq, Ipp64f rHighFreq, Ipp64f* taps,
#                                       int tapsLen, IppWinType winType, IppBool doNormal))
