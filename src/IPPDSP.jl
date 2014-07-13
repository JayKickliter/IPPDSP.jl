module IPPDSP

using IPPCore

IPPTypeDict = [ IPP8u   => "8u",
                IPP8s   => "8s",
                IPP16u  => "16u",
                IPP16s  => "16s",
                IPP32u  => "32u",
                IPP32s  => "32s",
                IPP16f  => "16f",
                IPP32f  => "32f",
                IPP64f  => "64f",
                IPP16uc => "16uc",
                IPP16sc => "16sc",
                IPP32uc => "32uc",
                IPP32sc => "32sc",
                IPP16fc => "16fc",
                IPP32fc => "32fc",
                IPP64fc => "64fc" ]

function IPPTypeSignature( IPPType )
    IPPTypeDict[ eval( IPPType ) ]
end

function IPPSuffix( IPPType )
   string( "_", IPPTypeSignature( IPPType ))
end

function IPPSuffix( IPPTypes::Tuple )
    lenTypes = length( IPPTypes )
    lenTypes == 0 && return ""
    lenTypes == 1 && return IPPSuffix( IPPTypes[1] )

    suffix = IPPTypeSignature( IPPTypes[1] )
    for Tsymb in IPPTypes[2:end]
        suffix *= typeof(Tsymb)<:String ? string( "_", Tsymb ) : IPPSuffix( Tsymb )
    end

    suffix
end

# common

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
