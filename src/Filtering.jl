# ippf  variable to hold the the IPP function name
# sr    single rate
# mr    multirate
# Tt    datatype of filter taps (coefficients)
# Tx    datatype of signal

const WIN_BARTLETT = 0
const WIN_BLACKMAN = 1
const WIN_HAMMING  = 2
const WIN_HANN     = 3
const WIN_RECT     = 4
const LOWPASS      = 0
const HIGHPASS     = 1
const BANDPASS     = 2
const BANDSTOP     = 3

export  conv,
        xcorr,
        autocorr,
        autocorrb,
        autocorru,
        FIRFilter,
        filt,
        FIRInit





################################################################################
#               ____ ____ _  _ _  _ ____ _    _  _ ___ _ ____ _  _             #
#               |    |  | |\ | |  | |  | |    |  |  |  | |  | |\ |             #
#               |___ |__| | \|  \/  |__| |___ |__|  |  | |__| | \|             #
################################################################################

conv_types = [ :IPP32f,
               :IPP64f,
               :IPP16s ]

for ( julia_fun, ippf_prefix ) in [ ( :conv, "ippsConv" ) ]

    julia_fun! = symbol(string(julia_fun, '!'))

    for T in conv_types

        ippf = string( ippf_prefix, IPPSuffix( T ) )

        if eval(T) == IPP16s # integer versions require a scaling factor
            ippf *= "_Sfs"
            @eval begin
                function $(julia_fun!)( y::Array{$T}, x1::Array{$T}, x2::Array{$T}; scale = 0  )
                    ny  = length( y )
                    nx1 = length( x1 )
                    nx2 = length( x2 )
                    ny == nx1+nx2-1 || error("Length(y) must equal length(x1)+length(x2)-1" )
                    if ny > 3
                        @ippscall( $ippf, ( Ptr{$T},    IPPInt,     Ptr{$T},    IPPInt,     Ptr{$T},    IPPInt  ),
                                            x1,         nx1,        x2,         nx2,        y,          scale   )
                    end
                    return y
                end

                $(julia_fun)(  x1::Array{$T}, x2::Array{$T}; args... ) = $(julia_fun!)( similar(x1, length(x1) + length(x2) - 1), x1, x2; args... )

            end
        else
            @eval begin
                function $(julia_fun!)( y::Array{$T}, x1::Array{$T}, x2::Array{$T}  )
                    ny  = length( y )
                    nx1 = length( x1 )
                    nx2 = length( x2 )
                    ny == nx1+nx2-1 || error( "length(y) must be length(x1)+length(x2)-1" )
                    if ny > 3
                        @ippscall( $ippf, ( Ptr{$T},    IPPInt,     Ptr{$T},    IPPInt,     Ptr{$T} ),
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

xcorr_types = [ :IPP32f,
                :IPP64f,
                :IPP32fc,
                :IPP64fc,
                :IPP16s  ]

for ( julia_fun, ippf_prefix ) in [ ( :xcorr, "ippsCrossCorr" ) ]

    julia_fun! = symbol(string(julia_fun, '!'))

    for T in xcorr_types

        ippf = string( ippf_prefix, IPPSuffix( T ) )

        if eval(T) == IPP16s # scaled version
            ippf *= "_Sfs"
            @eval begin
                function $(julia_fun!)( y::Array{$T}, x1::Array{$T}, x2::Array{$T}; scale = 0, lowlag = 0)
                    ny  = length( y )
                    nx1 = length( x1 )
                    nx2 = length( x2 )
                    ny > 1 && nx1 > 1 && nx2 > 1 || error()
                    if ny > 3
                        @ippscall( $ippf, ( Ptr{$T},    IPPInt,     Ptr{$T},    IPPInt,     Ptr{$T},    IPPInt, IPPInt,  IPPInt ),
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
                    ny > 1 && nx1 > 1 && nx2 > 1 || error() # TODO: Check constrainTx
                    if ny > 3
                        @ippscall( $ippf, ( Ptr{$T},    IPPInt,     Ptr{$T},    IPPInt,     Ptr{$T},    IPPInt, IPPInt  ),
                                            x1,         nx1,        x2,         nx2,        y,          ny,     lowlag  )
                    end
                    return y
                end

                $(julia_fun)(  x1::Array{$T}, x2::Array{$T}; args... ) = $(julia_fun!)( similar(x1, length(x1) + length(x2) - 1), x1, x2; args... )
            end
        end
    end
end




autocorr_types = [ :IPP16s,
                   :IPP32f,
                   :IPP64f,
                   :IPP16fc,
                   :IPP32fc,
                   :IPP64fc  ]


for ( julia_fun, ippf_prefix )  in  [   ( :autocorr,  "ippsAutoCorr"       ),
                                        ( :autocorrb, "ippsAutoCorr_NormA" ),
                                        ( :autocorru, "ippsAutoCorr_NormB" )   ]

    julia_fun! = symbol(string(julia_fun, '!'))

    for T = autocorr_types

        ippf = string( ippf_prefix, '_', IPPSuffix( T ) )

        if eval(T) == IPP16s
            ippf *= "_Sfs"
            @eval begin
                function $(julia_fun!)( y::Array{$T}, x::Array{$T}; scale = 0 )
                    ny = length( y )
                    nx = length( x )
                    ny == nx || throw( DimensionMismatch("length(y) != length(x) ") )
                    if ny > 0
                        @ippscall( $ippf, ( Ptr{$T},    IPPInt, Ptr{$T},    IPPInt, IPPInt  ),
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
                        @ippscall( $ippf, ( Ptr{$T},    IPPInt, Ptr{$T},    IPPInt  ),
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
#      ____ _ ____    ____ _ _    ___ ____ ____    ___ _   _ ___  ____         #
#      |___ | |__/    |___ | |     |  |___ |__/     |   \_/  |__] |___         #
#      |    | |  \    |    | |___  |  |___ |  \     |    |   |    |___         #
################################################################################
# These are common to most of the ipps FIR Filter functions
# I haven't included functions that require integer scaling.
# They take difference parameters. Hope I haven't painted
# myself into a corner.
FIRFilterTypes =    [   ( :IPP32f,      :IPP32f  ),
                        ( :IPP32fc,     :IPP32fc ),
                        ( :IPP64f,      :IPP64f  ),
                        ( :IPP64fc,     :IPP64fc ),
                        ( :IPP32f,      :IPP16s  ),
                        ( :IPP64f,      :IPP16s  ),
                        ( :IPP64f,      :IPP32f  ),
                        ( :IPP64f,      :Int32   ),
                        ( :IPP64fc,     :IPP16sc ),
                        ( :IPP64fc,     :IPP32sc ),
                        ( :IPP64fc,     :IPP32fc )   ]

# IPP assings a pointer to state
# This points to an invisuble-to-us struct in buffer
# Intel's documentation is weak in this area, but we allocate the buffer, they
# instantsiat a stuct inside that buffer, and use state to point the location of
# that state

# Both the singe-rate and multirate states are the same, but are in different
# types to help with dispatching

# Single rate state
type FIRSRState
    pointer::Ptr{Void}
    buffer::Vector{Uint8}
end

# Multirate state
type FIRMRState
    pointer::Ptr{Void}
    buffer::Vector{Uint8}
end

# Polyphase state
type FIRPPState
    pointer::Ptr{Void}
    buffer::Vector{Uint8}
end

# FIRFilter type, its type signature includes the FIRState type, which is how functions taking
# FIRFilter as a parapter get diuspatched to the correct single-rate or multirate IPP functions
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
    pointer   = convert( Ptr{Uint8}, 0 )
    bufLen    = FIRRequiredStateSize( Tt, Tx, tapsLen )
    delayLine = zeros( Tx, tapsLen )
    buffer    = zeros( Uint8, bufLen )
    state     = FIRSRState( pointer, buffer )
    FIRInit!( state, taps, delayLine )
    FIRFilter( taps, delayLine, state, 1, 1, 0, 0 )
end

# Multirate constructor
function FIRFilter{ Tx, Tt }( ::Type{Tx}, taps::Array{Tt, 1}, upFactor, downFactor, upPhase = 0, downPhase = 0 )
    tapsLen   = length( taps )
    pointer   = convert( Ptr{Uint8}, 0 )
    bufLen    = FIRRequiredStateSize( Tt, Tx, tapsLen, upFactor, downFactor )
    delayLine = zeros( Tx, tapsLen )
    buffer    = zeros( Uint8, bufLen )
    state     = FIRMRState( pointer, buffer )
    FIRInit!( state, taps, delayLine, upFactor, downFactor, upPhase, downPhase )
    FIRFilter( taps, delayLine, state, upFactor, downFactor, upPhase, downPhase )
end

FIRFilter{ Tx, Tt }( ::Type{Tx}, taps::Array{Tt, 1}, resampleRatio::Rational, upPhase = 0, downPhase = 0 ) = FIRFilter( Tx, taps, num(resampleRatio), den(resampleRatio), upPhase, downPhase )




################################################################################
#           ____ _ ____    ____ _ _    ___ ____ ____    _ _  _ _ ___           #
#           |___ | |__/    |___ | |     |  |___ |__/    | |\ | |  |            #
#           |    | |  \    |    | |___  |  |___ |  \    | | \| |  |            #
################################################################################
#   Name:         ippsFIRGetStateSize, ippsFIRMRGetStateSize,
#                 ippsFIRInit, ippsFIRMRInit
#   Purpose:      ippsFIRGetStateSize   - calculates the size of the FIR Statestructure;
#                 ippsFIRInit           - initialize FIR state - set taps and delay line
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
    for ( Tt, Tx ) in FIRFilterTypes

        ippf_suffix = Tt == Tx ? IPPSuffix( Tt ) : IPPSuffix( (Tt, Tx) )

        # Single rate version
        ippfsr = string( ippf_prefix1, ippf_prefix2, ippf_suffix )
        @eval begin
            function $(julia_fun)( ::Type{$Tt}, ::Type{$Tx}, tapslen::Integer  )
                buffersize = IPPInt[0]  # the function will fill this value
                tapslen > 0 || error()
                @ippscall( $ippfsr, ( IPPInt,     Ptr{IPPInt} ),
                                      tapslen,    buffersize    )

                return buffersize[1]
            end
        end

        # Multirate version
        ippfmr = string( ippf_prefix1, "MR", ippf_prefix2, ippf_suffix )
        @eval begin
            function $(julia_fun)( ::Type{$Tt}, ::Type{$Tx}, tapslen::Integer, upFactor::Integer, downFactor::Integer )
                buffersize = IPPInt[0]  # the function will fill this value
                tapslen > 0 || error() # todo: check this
                @ippscall( $ippfmr, ( IPPInt,   IPPInt,    IPPInt,      Ptr{IPPInt} ),
                                      tapslen,  upFactor,  downFactor,  buffersize    )

                return buffersize[1]
            end
            $(julia_fun)( ::Type{$Tt}, ::Type{$Tx}, tapslen::Integer, resamp_ratio::Rational ) = $(julia_fun)( $Tt, $Tx, tapslen, num(resamp_ratio), den(resamp_ratio) )
        end
    end
end

for ( julia_fun, ippf_prefix1, ippf_prefix2 )  in  [ (   :FIRInit,  "ippsFIR",  "Init"  ) ]
    for ( Tt, Tx ) in FIRFilterTypes

        ippf_suffix = Tt == Tx ? IPPSuffix( Tt ) : IPPSuffix( (Tt, Tx) )
        julia_fun! = symbol(string(julia_fun, '!'))

        # Single rate version
        ippfsr     = string( ippf_prefix1, ippf_prefix2, ippf_suffix )
        @eval begin

            function $(julia_fun!)( state::FIRSRState, taps::Vector{ $Tt }, delayLine::Vector{ $Tx } )
                ntaps          = length( taps )
                ndelayline     = length( delayLine )
                statePtr       = Array( Ptr{Void}, 1 )
                @ippscall( $ippfsr, (   Ptr{Void},         Ptr{$Tt},    IPPInt,     Ptr{$Tx},      Ptr{Uint8}   ),
                                        pointer(statePtr), taps,        ntaps,      delayLine,     state.buffer )
                state.pointer = statePtr[1]
                nothing
            end

        end

        # Multirate version
        ippfmr  = string( ippf_prefix1, "MR", ippf_prefix2, ippf_suffix )
        @eval begin

            function $(julia_fun!)(  state::FIRMRState, taps::Vector{ $Tt }, delayLine::Vector{ $Tx }, upFactor, downFactor, upPhase, downPhase )
                ntaps          = length( taps )
                ndelayline     = length( delayLine )
                statePtr       = Array( Ptr{Void}, 1 )
                0 <= upPhase < upFactor && 0 <= downPhase < downFactor || error()
                @ippscall( $ippfmr, (   Ptr{Void},         Ptr{$Tt}, IPPInt, IPPInt,    IPPInt,  IPPInt,     IPPInt,    Ptr{$Tx},  Ptr{Uint8}   ),
                                        pointer(statePtr), taps,     ntaps,  upFactor,  upPhase, downFactor, downPhase, delayLine, state.buffer )
                state.pointer = statePtr[1]
                nothing
            end

        end
    end
end




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
    for ( Tt, Tx ) in FIRFilterTypes

        ippf_suffix = Tt == Tx ? IPPSuffix( Tt ) : IPPSuffix( (Tt, Tx) )
        julia_fun!  = symbol(string(julia_fun, '!'))

        # single rate
        ippfsr = string( ippf_prefix, ippf_suffix )
        @eval begin

            # TODO: make a lower level version of filt so direct form filters dont need to create a fir filter object
            #       they should only need to call getsize and init
            # using a pre-instanciated FIRFilter object, and a preallocated buffer
            function $(julia_fun!)( self::FIRFilter{$Tt, $Tx, FIRSRState}, buffer::Vector{ $Tx }, signal::Vector{ $Tx } )
                sigLen = length( signal )
                bufLen = length( buffer )
                bufLen >= sigLen || error( "output buffer is too small" )
                @ippscall( $ippfsr,  (  Ptr{$Tx},       Ptr{$Tx},       IPPInt,    Ptr{Void}                    ),
                                        signal,         buffer,         sigLen,    self.state.pointer           )
                buffer
            end

            # using a pre-instanciated FIRFilter object
            $(julia_fun)( self::FIRFilter{$Tt, $Tx, FIRSRState}, signal::Vector{$Tx} ) = $(julia_fun!)( self, similar( signal ), signal )

            # direct form, with taps and signal, and a preallocated buffer
            function $(julia_fun!)( buffer::Vector{ $Tx }, taps::Vector{ $Tt }, signal::Vector{ $Tx })
                self = FIRFilter( $Tx, taps )
                filt!( self, buffer, signal )
            end

            # direct form, with taps and signal
            $(julia_fun)( taps::Vector{ $Tt }, signal::Vector{ $Tx } ) = $(julia_fun!)( similar(signal), taps, signal )

        end

        # multirate
        ippfmr = string( ippf_prefix, "MR", ippf_suffix )
        @eval begin

            # using a pre-instanciated FIRFilter object, and a preallocated buffer
            function $(julia_fun!)( self::FIRFilter{ $Tt, $Tx, FIRMRState }, buffer::Vector{ $Tx }, signal::Array{$Tx} )
                sigLen = length( signal )
                outLen = length( buffer )
                sigLen % self.downFactor == 0 || error( "signal length must be an integer multiple of decimation")
                sigLen * self.upFactor == outLen * self.downFactor || error( "buffer length needs to be >= sigal len * interpolation/decimation" )
                iterations  = int( sigLen/self.downFactor )
                @ippscall( $ippfsr,  (  Ptr{$Tx},  Ptr{$Tx},  IPPInt,     Ptr{Void}            ),
                                        signal,    buffer,    iterations, self.state.pointer   )
                buffer
            end

            # using a pre-instanciated FIRFilter object
            function $(julia_fun)( self::FIRFilter{$Tt, $Tx, FIRMRState}, signal::Vector{ $Tx } )
                sigLen = length( signal )
                bufLen = int( sigLen * self.upFactor / self.downFactor )
                buffer = similar( signal, bufLen )
                $(julia_fun!)( self, buffer, signal )
            end

            # direct form, with taps and signal, and a preallocated buffer
            function $(julia_fun!)( buffer::Vector{ $Tx }, taps::Vector{ $Tt }, signal::Vector{ $Tx }, upFactor::Integer, downFactor::Integer, upPhase = 0, downPhase = 0)
                self = FIRFilter( $Tx, taps, upFactor, downFactor, upPhase, downPhase )
                $(julia_fun!)( self, buffer, signal )
            end

            # direct form, with taps and signal
            function $(julia_fun)( taps::Vector{ $Tt }, signal::Vector{ $Tx }, upFactor::Integer, downFactor::Integer, upPhase = 0, downPhase = 0)
                sigLen = length( signal )
                bufLen = int( sigLen * upFactor / downFactor )
                buffer = similar( signal, bufLen )
                $(julia_fun!)( buffer, taps, signal, upFactor, downFactor, upPhase, downPhase )
            end
        end
    end
end



######################################################################################
# ____ _ ____    ___ ____ ___  ____    ____ ____ _  _ ____ ____ ____ ___ _ ____ _  _ #
# |___ | |__/     |  |__| |__] [__     | __ |___ |\ | |___ |__/ |__|  |  | |  | |\ | #
# |    | |  \     |  |  | |    ___]    |__] |___ | \| |___ |  \ |  |  |  | |__| | \| #
######################################################################################
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
#      winType           the ippwindow switch variable,
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
#  IPPAPI(IppStatus, ippsFIRGenLowpass_64f( Ipp64f rfreq,    Ipp64f* taps,      int tapsLen, IppWinType winType, IppBool doNormal))
#  IPPAPI(IppStatus, ippsFIRGenHighpass_64f(Ipp64f rfreq,    Ipp64f* taps,      int tapsLen, IppWinType winType, IppBool doNormal))
#  IPPAPI(IppStatus, ippsFIRGenBandpass_64f(Ipp64f rLowFreq, Ipp64f  rHighFreq, Ipp64f* taps, int tapsLen, IppWinType winType, IppBool doNormal))
#  IPPAPI(IppStatus, ippsFIRGenBandstop_64f(Ipp64f rLowFreq, Ipp64f  rHighFreq, Ipp64f* taps, int tapsLen, IppWinType winType, IppBool doNormal))

for ( julia_fun, ippf_prefix )  in  [   ( :highpass,  "ippsFIRGenHighpass" ),
                                        ( :lowpass,   "ippsFIRGenLowpass"  )  ]
    julia_fun! = symbol(string(julia_fun, '!'))

    for T in [ IPP64f ]

        ippf_suffix = IPPSuffix( T )
        ippf        = string( ippf_prefix, ippf_suffix )

        @eval begin
            function $(julia_fun!){Tt}( buffer::Vector{Tt}, ω::FloatingPoint, window = WIN_HAMMING, normalize = true )
                Ntaps = length ( buffer )
                0.0 < ω < 0.5   || error( "cutoff frequency ω must satisfy:  0.0 < ω < 0.5" )
                Ntaps >= 5      || error( "number of taps must satisfy: Ntaps >= 5")
                innerBuffer = Tt == $T ? buffer : Array( $T, Ntaps )
                @ippscall( $ippf,  (  IPP64f, Ptr{$T},      IPPInt, IPPInt,     Bool        ),
                                      ω,      innerBuffer,  Ntaps,  window, normalize   )
                return Tt == $T ? buffer : [ buffer[i] = convert( Tt, innerBuffer[i] ) for i = 1:Ntaps ]
            end # function

            $(julia_fun){Tt}( ::Type{Tt}, Ntaps::Integer, ω::FloatingPoint, window = WIN_HAMMING, normalize = true ) = $(julia_fun!)( Array(Tt, Ntaps), ω, window, normalize )
        end # eval
    end # type loop
end # function loop

for ( julia_fun, ippf_prefix )  in  [   ( :bandpass,  "ippsFIRGenBandpass"  ),
                                        ( :bandstop,  "ippsFIRGenBandstop"  )  ]

    julia_fun! = symbol(string(julia_fun, '!'))

    for T in [ IPP64f ]

        ippf_suffix = IPPSuffix( T )
        ippf        = string( ippf_prefix, ippf_suffix )

        @eval begin
            function $(julia_fun!){Tt}( buffer::Vector{Tt}, ω1::FloatingPoint, ω2::FloatingPoint, window = WIN_HAMMING, normalize = true )
                Ntaps = length ( buffer )
                0.0 < ω1 < ω2 < 0.5 || error( "cutoff frequencies ω1, ω2  must satisfy:  0.0 < ω1 < ω2 < 0.5" )
                Ntaps >= 5          || error( "number of taps must satisfy: Ntaps >= 5" )
                innerBuffer = Tt == $T ? buffer : Array( $T, Ntaps )
                @ippscall( $ippf,  (  IPP64f, IPP64f, Ptr{$T},      IPPInt, IPPInt,     Bool        ),
                                      ω1,     ω2,     innerBuffer,  Ntaps,  window, normalize   )
                return Tt == $T ? buffer : [ buffer[i] = innerBuffer[i] for i = 1:Ntaps ]
            end # function

            $(julia_fun){Tt}( ::Type{Tt}, Ntaps::Integer, ω1::FloatingPoint, ω2::FloatingPoint, window = WIN_HAMMING, normalize = true ) = $(julia_fun!)( Array(Tt, Ntaps), ω1, ω2, window, normalize )
        end # eval
    end # type loop
end # function loop

