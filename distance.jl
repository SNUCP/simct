include("randbfv/RandBfv.jl")
using HIENAA, Printf

function deterministic(x::BFV, y::Vector{UInt64}, rlk::RLEV, rtk::Vector{RLEV}, param::BFVParameters)::BFV
    oper = BFVOperator(param)

    ctout = similar(x)
    tmp = similar(x)
    pty = PlainPoly(param.ring_param.N, length(param.Q))

    @time begin
        # (x - y)^2
        encode_to!(pty, y, oper)
        sub_to!(ctout, x, pty, oper)
        mul_to!(ctout, ctout, ctout, rlk, oper)

        # Rotsum, assuming that the rotation group is a subgroup of <5>.
        for i = length(rtk)-1:-1:0
            rotate_to!(tmp, ctout, (1 << i, 0), rtk[i+1], oper)
            add_to!(ctout, ctout, tmp, oper)
        end
    end

    ctout
end

function randomised(x::BFV, y::Vector{UInt64}, rlk::RLEV, rtk::Vector{RLEV}, pk::RLWE, param::BFVParameters)
    # Parameters
    σ = 48.373546489791295
    τ_cmul = 1.992990407684452e13

    # Generate samplers and randomizer.
    oper = BFVOperator(param)
    pkentor = PKEncryptor(pk, 91.13921251821257, 7383.043789707919, oper.operQ)

    dgs_pmul = TwinCDTSampler(σ)
    rand = RandBFVOperator(oper, pkentor)

    ctout = similar(x)
    tmp = similar(x)
    pty = PlainPoly(param.ring_param.N, length(param.Q))

    @time begin
        # (x - y)^2
        encode_to!(pty, y, oper)
        sub_to!(ctout, x, pty, oper)
        mul_to!(ctout, ctout, ctout, rlk, dgs_pmul, τ_cmul, τ_cmul, rand)

        # rotsum
        # rotation group <5>
        for i = length(rtk)-1:-1:0
            rotate_to!(tmp, ctout, (1 << i, 0), rtk[i+1], rand)
            add_to!(ctout, ctout, tmp, rand)
        end

        mask_randomize!(ctout, rand)
    end

    ctout
end

function optimised(x::BFV, y::Vector{UInt64}, rlk::RLEV, rtk::Vector{RLEV}, pk::RLWE, param::BFVParameters)
    # Parameters
    σ = 48.373546489791295
    τ_cmul = 1.992990407684452e13

    # Generate samplers and randomizer.
    oper = BFVOperator(param)
    pkentor = PKEncryptor(pk, 91.13921251821257, 7383.043789707919, oper.operQ)

    dgs_pmul = TwinCDTSampler(σ)
    rand = RandBFVOperator(oper, pkentor)

    ctout = similar(x)
    tmp = similar(x)
    pty = PlainPoly(param.ring_param.N, length(param.Q))
    buff = similar(y)

    @time begin
        evalT, oper = rand.evalT, rand.oper

        # -2xy
        mul_to!(buff, evalT.Q.Q - 2, y, evalT)
        encode_to!(pty, buff, dgs_pmul, rand)
        mul_to!(ctout, x, pty, τ_cmul, rand)

        # tmp = x^2 (deterministic)
        mul_to!(tmp, x, x, rlk, oper)

        # ctout += tmp
        add_to!(ctout, ctout, tmp, rand)

        # ctout += y^2
        encode_to!(pty, y, oper)
        mul_to!(pty, pty, pty, oper.operQ)
        add_to!(ctout, ctout, pty, rand)

        # rotsum
        # rotation group <5>
        mask_randomize!(ctout, rand)
        for i = length(rtk)-1:-1:0
            rotate_to!(tmp, ctout, (1 << i, 0), rtk[i+1], oper)
            add_to!(ctout, ctout, tmp, oper)
        end

        mask_randomize!(ctout, rand)
    end

    ctout
end

function deterministic_flood(x::BFV, y::Vector{UInt64}, rlk::RLEV, rtk::Vector{RLEV}, param::BFVParameters, pk::RLWE)
    # Generate operator and encryptor.
    oper = BFVOperator(param)
    pkentor = PKEncryptor(pk, 91.13921251821257, 7383.043789707919, oper.operQ)

    ctout = similar(x)
    tmp = similar(x)
    pty = PlainPoly(param.ring_param.N, length(param.Q))

    @time begin
        # (x - y)^2
        encode_to!(pty, y, oper)
        sub_to!(ctout, x, pty, oper)
        mul_to!(ctout, ctout, ctout, rlk, oper)

        # Rotsum, assuming that the rotation group is a subgroup of <5>.
        for i = length(rtk)-1:-1:0
            rotate_to!(tmp, ctout, (1 << i, 0), rtk[i+1], oper)
            add_to!(ctout, ctout, tmp, oper)
        end
        flood_to!(ctout, oper, pkentor)
    end

    ctout
end

function flood_to!(x::BFV, oper::BFVOperator, pkentor::PKEncryptor)
    # Parameters
    bits = 154

    evalQ, buff = oper.operQ.evalQ, oper.ct_buff[end][1:length(x.val)]

    # mask randomize
    rlwe_sample_to!(buff.val, pkentor)
    buff.level[] = x.level[]
    add_to!(x, x, buff, oper)

    # error randomize
    rng = pkentor.rgsampler.rng
    bound = big(1) << bits
    error = ModPoly(rand(rng, -bound:bound, oper.operQ.param.N), evalQ)
    add_to!(x.val.b, x.val.b, error, evalQ)
end

function get_meanstdmax(err::Vector{BigInt})
    mean = sum(err) / length(err)
    std = sqrt(sum([(e - mean)^2 for e = err]) / length(err))
    max = maximum(abs.(err))

    mean, std, max
end

function main()
    println("IN PREPARATION...")

    # Scheme parameters
    m = 1 << 13
    P, Q = missing, UInt64[0x0040000000006001, 8404993*2147565569]
    dlen, t, ispacking, islevelled = 1, 8404993, true, false
    packlen = 128

    # Generate the basic structs.
    ring_param = CyclotomicParam(m)
    param = BFVParameters(ring_param, P, Q, dlen, t, ispacking, islevelled)
    scheme = BFVScheme(param)

    # Generate the secret and evaluation keys.
    us = UniformSampler()
    sk = ternary_ringkey(us, ring_param.N)
    set_encryptor!(sk, scheme)
    rlk = relin_keygen(scheme)
    rtk = [rotate_keygen((1 << i, 0), scheme) for i = 0:trailing_zeros(packlen)-1]

    # Generate the public key.
    pk = public_keygen(scheme)

    # Generate the encryptions.
    msg = UInt64[i % 256 for i = 1:packlen]
    packmsg = repeat(msg, (m >> 1) ÷ packlen)   # Sparse packing.

    ptx = encode(packmsg, scheme)
    ctx = encrypt(ptx, scheme)

    y = reverse(packmsg)

    #===========================================#

    # DETERMINISTIC COMPUTATION

    print("DETERMINISTIC COMPUTATION : ")

    ct_deterministic = deterministic(ctx, y, rlk, rtk, param)
    out = get_error(ct_deterministic, scheme)
    mean, std, max = get_meanstdmax(out)

    @printf("Mean : %.3f-bits, Standard Deviation : %.3f-bits, Max : %.3f-bits \n\n", log2(abs(mean)), log2(std), log2(max))

    #===========================================#

    # RANDOMISED COMPUTATION

    print("RANDOMISED COMPUTATION : ")

    ct_randomised = randomised(ctx, y, rlk, rtk, pk, param)
    out = get_error(ct_randomised, scheme)
    mean, std, max = get_meanstdmax(out)

    @printf("Mean : %.3f-bits, Standard Deviation : %.3f-bits, Max : %.3f-bits \n\n", log2(abs(mean)), log2(std), log2(max))

    #===========================================#

    # OPTIMISED COMPUTATION

    print("OPTIMISED COMPUTATION : ")

    ct_optimised = optimised(ctx, y, rlk, rtk, pk, param)
    out = get_error(ct_optimised, scheme)
    mean, std, max = get_meanstdmax(out)

    @printf("Mean : %.3f-bits, Standard Deviation : %.3f-bits, Max : %.3f-bits \n\n", log2(abs(mean)), log2(std), log2(max))
    
    nothing
end

function main2()
    println("IN PREPARATION...")

    # Scheme parameters
    m = 1 << 14
    P, Q = missing, UInt64[0x0000800000020001, 0x000080000008c001, 0x00008000000dc001, 0x00008000000f4001]
    dlen, t, ispacking, islevelled = 1, 8404993, true, false
    packlen = 128

    # Generate the basic structs.
    ring_param = CyclotomicParam(m)
    param = BFVParameters(ring_param, P, Q, dlen, t, ispacking, islevelled)
    scheme = BFVScheme(param)

    # Generate the secret and evaluation keys.
    us = UniformSampler()
    sk = ternary_ringkey(us, ring_param.N)
    set_encryptor!(sk, scheme)
    rlk = relin_keygen(scheme)
    rtk = [rotate_keygen((1 << i, 0), scheme) for i = 0:trailing_zeros(packlen)-1]

    # Generate the public key.
    pk = public_keygen(scheme)
    
    # Generate the encryptions.
    msg = UInt64[i % 256 for i = 1:packlen]
    packmsg = repeat(msg, (m >> 1) ÷ packlen)   # Sparse packing.

    ptx = encode(packmsg, scheme)
    ctx = encrypt(ptx, scheme)

    y = reverse(packmsg)

    #===========================================#

    # DETERMINISTIC COMPUTATION

    print("DETERMINISTIC COMPUTATION : ")

    ct_flood = deterministic_flood(ctx, y, rlk, rtk, param, pk)
    out = get_error(ct_flood, scheme)
    mean, std, max = get_meanstdmax(out)

    @printf("Mean : %.3f-bits, Standard Deviation : %.3f-bits, Max : %.3f-bits \n\n", log2(abs(mean)), log2(std), log2(max))
end

main()
main2()