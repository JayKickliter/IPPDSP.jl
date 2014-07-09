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
