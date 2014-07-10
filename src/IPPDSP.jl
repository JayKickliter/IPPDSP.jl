module IPPDSP

using IPPCore
import IPPCore: IppInt, ippint


# common

typealias IPP16s    Int16
typealias IPP32s    Int32
typealias IPP16f    Float16
typealias IPP32f    Float32
typealias IPP64f    Float64
typealias IPP16sc   Complex{Int16}
typealias IPP32sc   Complex{Int32}
typealias IPP16fc   Complex{Float16}
typealias IPP32fc   Complex{Float32}
typealias IPP64fc   Complex{Float64}

ipp16s( x )  = convert( IPP16s, x ) 
ipp32s( x )  = convert( IPP32s, x ) 
ipp16f( x )  = convert( IPP16f, x ) 
ipp32f( x )  = convert( IPP32f, x ) 
ipp64f( x )  = convert( IPP64f, x ) 
ipp16sc( x ) = convert( IPP16sc, x )
ipp32sc( x ) = convert( IPP32sc, x )
ipp16fc( x ) = convert( IPP16fc, x )
ipp32fc( x ) = convert( IPP32fc, x )
ipp64fc( x ) = convert( IPP64fc, x )




macro ippscall(ippf, argtypes, args...)
    quote
        ret = ccall(($ippf, "libipps"), IppStatus,
            $argtypes, $(args...))
        ret == 0 || error(ippstatus_string(ret))
    end
end

# source files

include("Windows.jl")
include("Filtering.jl")

end # module
