
#################################################################################

"""
    const _Float64_threshold = Ref{Float64}(1e-15)

Holds threshold for floating point calculations. In most places, any number below this will be
treated zero. E.g., in `cutoff`.
"""
const _Float64_threshold = Ref{Float64}(1e-15)

"""
    Float64_threshold()
    Float64_threshold(threshold::Float)

Returns threshold for floating point calculations. In most places, any number below this will be
treated zero. E.g., in `cutoff`. Default: `1e-15`.

Optionally, change the threshold as per the input `threshold`.
"""
Float64_threshold() = _Float64_threshold[]

function Float64_threshold(threshold::Float64)
    _Float64_threshold[] = threshold
    return _Float64_threshold[]
end
    
#################################################################################

"""
    const _using_threaded_loop = Ref{Bool}(true)

Holds the condition whether to use threaded loop.
"""
const _using_threaded_loop = Ref{Bool}(true)

"""
    using_threaded_loop()

Returns the condition whether to use threaded loop.
"""
using_threaded_loop() = _using_threaded_loop[]

"""
    enable_threaded_loop()

Enable the use of threaded loop.
"""
function enable_threaded_loop()
    _using_threaded_loop[] = true
    return nothing
end

"""
    disable_threaded_loop()

Disable the use of threaded loop.
"""
function disable_threaded_loop()
    _using_threaded_loop[] = false
    return nothing
end

#################################################################################

"""
    macro threaded_loop(code)

Conditional threaded loop by `Threads.@threads`. If `using_threaded_loop()` is true,
performs parallel loop, otherwise performs serial loop.
"""
macro threaded_loop(code)
    return esc(:(
        if using_threaded_loop() && Threads.nthreads() > 1

            blas_nthreads = BLAS.get_num_threads()
            blas_nthreads > 1 && BLAS.set_num_threads(1)

            using_threaded_bsp = ITensors.using_threaded_blocksparse()
            using_threaded_bsp && ITensors.disable_threaded_blocksparse()

            strd_nthreads = ITensors.Strided.get_num_threads()
            strd_nthreads > 1 && ITensors.Strided.set_num_threads(1)
            
            Threads.@threads($code)

            blas_nthreads > 1 && BLAS.set_num_threads(blas_nthreads)
            using_threaded_bsp && ITensors.enable_threaded_blocksparse()
            strd_nthreads > 1 && ITensors.Strided.set_num_threads(strd_nthreads)
        else
            $code
        end
    ))
end

#################################################################################
