struct RandBFVParamSketch <: HEParamSketch
    ring_param::RingParam
    logP::Int64
    logQ::Int64
    dlen::Int64
    ptxt_modulus::Int64
    ispacking::Bool

    RandBFVParamSketch(ring_param::RingParam, logP::Int64, logQ::Int64, ptxt_modulus::Int64; dlen::Int64=0, ispacking::Bool=true)::BFVParamSketch =
        new(ring_param, logP, logQ, dlen, ptxt_modulus, ispacking)
end

(::Type{BFVParameters})(sketch::RandBFVParamSketch)::BFVParameters = begin
    ring_param, logP, logQ, dlen, t, ispacking = sketch.ring_param, sketch.logP, sketch.logQ, sketch.dlen, sketch.ptxt_modulus, sketch.ispacking

    if logP == 0
        P = missing

        Q0 = t^floor(Int64, 62 / log2(t)) % UInt64
        Qlen = ceil(Int64, (logQ - log2(Q0)) / 62)

        Qbits = (logQ - log2(Q0)) / Qlen
        Qprimes = find_prime(ring_param, Qbits, Qlen)
        Q = vcat(Q0, Qprimes)

        if dlen == 0
            dlen = 1
        end
    else
        maxbits = min(62, logP) # Minimise the key-switching error.
        Plen = ceil(Int64, logP / 62)

        Q0 = t^floor(Int64, 62 / log2(t)) % UInt64
        Qlen = ceil(Int64, (logQ - log2(Q0)) / maxbits)

        Qbits = (logQ - log2(Q0)) / Qlen
        Qprimes = find_prime(ring_param, Qbits, Plen + Qlen)
        
        Pbits = logP / Plen
        P = find_prime(ring_param, Pbits, Plen)

        filter!(x -> x ∉ P, Qprimes)
        Q = vcat(Q0, Qprimes[1:Qlen])

        if dlen == 0
            dlen = floor(Int64, logP / Qbits)
        end
    end

    new(ring_param, P, Q, dlen, t, ispacking, false)
end