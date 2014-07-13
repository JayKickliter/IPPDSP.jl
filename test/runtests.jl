using IPPCore
import IPPDSP
using Base.Test

# write your own tests here

for T = IPPDSP.conv_types
    a = convert(Vector{eval(T)}, [ x = 1:10 ] )
    b = convert(Vector{eval(T)}, [ x = 11:21 ] )
    reference = conv( a, b )
    IPPResult = IPPDSP.conv( a, b )
    @test length( reference ) == length( IPPResult )
    for i = 1:length( reference )
        @test isapprox( reference[i], IPPResult[i] )
    end
end

for T = IPPDSP.xcorr_types
    a = convert(Vector{eval(T)}, [ x = 1:10 ] )
    b = convert(Vector{eval(T)}, [ x = 1:10 ] )
    reference = xcorr( a, b )
    IPPResult = IPPDSP.xcorr( a, b, lowlag = -9 )
    @test length( reference ) == length( IPPResult )
    for i = 1:length( reference )
        @test isapprox( reference[i], IPPResult[i] )
    end
end