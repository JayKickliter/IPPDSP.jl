module IPPDSP

using IPPCore
import IPPCore: IppInt, ippint


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
