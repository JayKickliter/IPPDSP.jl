# IPPDSP#

This package provides wrapper functions for Intel's [Integrated Performance Primitives](https://software.intel.com/en-us/intel-ipp). Specifically, **IPPDSP** targets libIPPS's [DSP](http://en.wikipedia.org/wiki/Digital_signal_processing) functions and data structures.

**IPPDSP** originally started as additions to Dahua Lin's **IPPMath**, as they both wrap `libIPPS`. Even though they use the same library, we decided that DSP specific functionality belonged in its own package. Special thanks to Dahua for laying a foundation for me. Until I read his code I was clueless when it came to using Julia's metaprogramming facilities.


## Installation Notes ##

### Julia Dependancies ###


**IPPDSP** depends on [**IPPCore**](https://github.com/lindahua/IPPCore.jl). It provides common functionality for Julia packages that wrap IPP libraries.


### Binary Dependancies ###


IPP is a commercial product. To use **IPPDSP**, you will need to manually install the IPP libraries.  You can download a 30 day trial for free. They are available for Windows, Linux, and OS X. The location of the libraries will need to be in either your OS library search path, or in Julia's global variable `DS_LOAD_PATH`. If Julia can't find `libIPPS`, add the command `push!(DL_LOAD_PATH, "/path/to/IPPLibs")` in `.juliarc.jl` located in your home folder. Here's the contents of my `.juliarc.jl` (OS X): 

```jl
push!(DL_LOAD_PATH, "/opt/intel/ipp/lib")
```

### Note for OS X Users ###


OS X dynamic libraries have their own paths, and the paths of their dependencies, hard-coded in the file. For some reason, IPP dylibs have their set relative to the library folder. So if you were to launch the Julia process from anywhere other than IPP's `lib` folder, **IPPDSP** will not work. There are three workarounds:

* Launch Julia from `/opt/intel/ipp/lib`
* Launch Julia, then run the command `cd("/opt/intel/ipp/lib")`
* Fix the library paths

I wrote a [Julia script](https://gist.github.com/JayKickliter/4c5c70f47be75a20e43e) that will fix the paths of all the IPP dylibs. If you have plans of using IPP for any purpose outside of **IPPDSP**, I recommend backing up `/opt/intel/ipp/lib` to an alternate location. I only say this because I don't know what Intel's tools expect for library paths.



### Package Installation ###

**IPPDSP** and **IPPCore** are currently unregistered. After installing the IPP libraries, run the following Julia command:

```jlcon
julia> Pkg.clone("https://github.com/JayKickliter/IPPCore.jl.git") # I have added some not-yet-merged features to IPPCore
julia> Pkg.clone("https://github.com/JayKickliter/IPPDSP.jl.git")
```
	
## Documentation ##

Notice that most of the functions listed take a limited number of argument types. That is because IPP provides plain C functions that are statically typed. There are a few places I've overridden IPP's static typing. These are functions that are typically only used as setup functions, like generating FIR filter coefficients.

IPP has depreciated in-place functions. Some not-in-place functions will work when passed the same pointer for source and destination buffer, but I haven't tested which ones work correctly yet. So all in-place functions you see here, denoted by a bang (`funcname!`), overwrite a vector provided by you. If you are doing the same operation over and over, and the output size will always be the same, using an in-place function with a pre-allocated buffer can substantially increase performance.

**Note:** the named argument `scale` only applies when the left (or both) argument is an integer vector. ( TODO: cite Intel's explanation of this )

### Convolution ###
 
#### In Place ####

```julia
conv!(  y::Vector{T}, x1::Vector{T}, x2::Vector{T}[, scale = 0 ])
```

**Where `T`:**

* IPP32f
* IPP64f
* IPP16s

#### Out of Place ####

```julia
y = conv( x1::Vector{T}, x2::Vector{T}[; scale = 0 ])
```

**Where `T`:**

* IPP32f
* IPP64f
* IPP16s






### Cross Correlation ###

#### In Place ####

```julia
xcorr!(  y::Vector{T}, x1::Vector{T}, x2::Vector{T}[, scale = 0 ])
```

#### Out of Place ####

```julia
y = xcorr( x1::Vector{T}, x2::Vector{T}[; scale = 0 ])
```

#### Valid Types `T` ####

* IPP32f
* IPP64f
* IPP16s



### Autocorrelation (Standard) ###

#### In Place ####

```julia
autocorr!(  y::Vector{T}, x1::Vector{T}, x2::Vector{T}[, scale = 0 ])
```

#### Out of Place ####

```julia
y = autocorr( x1::Vector{T}, x2::Vector{T}[; scale = 0 ])
```

#### Valid Types `T` ####

* IPP16s
* IPP32f
* IPP64f
* IPP16fc
* IPP32fc
* IPP64fc


### Autocorrelation (Normalized) ###

#### In Place ####

```julia
autocorrb!(  y::Vector{T}, x1::Vector{T}, x2::Vector{T}[, scale = 0 ])
```

#### Out of Place ####

```julia
y = autocorrb( x1::Vector{T}, x2::Vector{T}[; scale = 0 ])
```

#### Valid Types `T` ####

* IPP16s
* IPP32f
* IPP64f
* IPP16fc
* IPP32fc
* IPP64fc


### Autocorrelation (Un-normalized) ###

#### In Place ####

```julia
autocorru!(  y::Vector{T}, x1::Vector{T}, x2::Vector{T}[, scale = 0 ])
```

#### Out of Place ####

```julia
y = autocorru( x1::Vector{T}, x2::Vector{T}[; scale = 0 ])
```

#### Valid Types `T` ####

* IPP16s
* IPP32f
* IPP64f
* IPP16fc
* IPP32fc
* IPP64fc




## Filtering ##

### `FIRFilter` Object ###

Since IPP is plain C, you must specify ahead of time the data type of your signal.

These are all the valid combinations of `( tapsType, signalType )`:

```julia
( IPP32f,   :IPP32f  )
( IPP32fc,  :IPP32fc )
( IPP64f,   :IPP64f  )
( IPP64fc,  :IPP64fc )
( IPP32f,   :IPP16s  )
( IPP64f,   :IPP16s  )
( IPP64f,   :IPP32f  )
( IPP64f,   :Int32   )
( IPP64fc,  :IPP16sc )
( IPP64fc,  :IPP32sc )
( IPP64fc,  :IPP32fc )
```

Both single-rate and multirate `FIRFilter` objects maintains state, allowing for stream processing. As a matter of fact, intel has depreciated one-off FIR filter functions, requiring you to create a filter state. Although I haven't yet, I will implement a stateless fire-and-forget `filt( taps, signal )` function to quickly filter some data (TODO: don't forget this).

Back to stateful filtering, here's a typical scenerio:

1. Create `FIRFilter` object. Let's name it `myFilt`
2. Read `samples` from a stream of indefinite length
3. Call `filteredSamples = filt( myFilt, samples )`
4. Do something with `mySamples`
5. `GOTO 2`

### Single-Rate Constructor ###

Internally, IPP use polyphase deconstruction to efficiently filter with interpolation and decimation.

```julia
myFilter = FIRFilter( ::Type{Tx}, taps::Vector{Ttaps} )
```

#### Arguments ####

* `::Type{Tx}`: data type of the signal you will be filtering
* `taps`: your FIR filter coefficients

### Multirate Constructor ###


Internally, IPP use polyphase deconstruction to more efficiently filter with interpolation and decimation.

```julia
myFilter = FIRFilter( ::Type{Tx}, taps::Vector{Tt}, upFactor, downFactor, upPhase = 0, downPhase = 0 )
```

#### Arguments ####

* `::Type{Tx}`: The data type of the signal you will be filtering
* `taps`: your FIR filter coefficients
* `upfactor`: interpolation ratio
* `downfactor`: decimation ratio
* `upPhase`: (TODO: understand)
* `downPhase:` which of the `1:downfactor` samples keep for each output sample (TODO: verify)

### `filt` Function

#### Single Rate, In Place ####

```julia
buffer = similar( signal )
filt!( myFilt, buffer, signal )
```

#### Single Rate, Out of Place ####

```julia
y = filt( myFilt, signal )
```

#### Multirate, In Place ####

```julia
# TODO: elaborate the constrains of input/output lengths for multirate filters
sigLength = length( signal )
bufLength = int( sigLength * myFilt.upFactor/myFilt.downFactor )
buffer    = Array( Tx, bufLength )
filt!( myFilt, buffer, signal )
```

#### Multirate, Out of Place ####

```julia
y = filt!( myFilt, signal )
```