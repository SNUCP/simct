const Randomiser = Union{Vector{CDTSampler},TwinCDTSampler}

struct RandBFVOperator <: HEOperator
    oper::BFVOperator
    ct_buff::Vector{BFV}
    ptxt_buff::Vector{PlainPoly}
    poly_buff::Vector{Vector{UInt64}}
    buff_rand::Vector{Int64}
    evalT::PolyEvaluatorArb
    us::UniformSampler
    rgs::RGSampler
    entor::PKEncryptor

    @views function RandBFVOperator(oper::BFVOperator, entor::PKEncryptor)
        ct_buff = BFV[BFV(oper.operQ.param.N, length(oper.operQ.evalQ), 0) for _ = 1:2]
        ptxt_buff = PlainPoly[PlainPoly(oper.operQ.param.N, length(oper.operQ.evalQ)) for _ = 1:2]
        poly_buff = [Vector{UInt64}(undef, oper.operQ.param.N) for _ = 1:2]
        buff_rand = Vector{Int64}(undef, oper.operQ.param.N)
        evalT = PolyEvaluatorArb(oper.operQ.auxeval, oper.ptxt_modulus)
        us = UniformSampler()
        rgs = RGSampler()

        new(oper, ct_buff, ptxt_buff, poly_buff, buff_rand, evalT, us, rgs, entor)
    end
end

#======================================================================================#
################################## ENCODE AND DECODE ###################################
#======================================================================================#

function encode(msg::Vector{UInt64}, dgs::Randomiser, rand::RandBFVOperator; level::Integer=typemax(Int64), ispacking::Bool=true, isPQ::Bool=false, scaling_factor::Real=0.0)::PlainPoly
    oper = rand.oper
    Qatlevel = oper.Qatlevel
    Qlen = Qatlevel[end]

    if isPQ
        if ismissing(oper.operQ.evalP)
            throw(ErrorException("The parameter does not support PQ."))
        end
        res = PlainPoly(oper.operQ.param.N, Qlen + length(oper.operQ.evalP))
        encode_to!(res, msg, dgs, rand, isPQ=true, ispacking=ispacking)
    else
        res = PlainPoly(oper.operQ.param.N, Qlen)
        encode_to!(res, msg, dgs, rand, isPQ=false, ispacking=ispacking)
    end

    res
end

function encode_to!(res::PlainPoly, msg::Vector{UInt64}, dgs::Randomiser, rand::RandBFVOperator; level::Integer=typemax(Int64), ispacking::Bool=true, isPQ::Bool=false, scaling_factor::Real=0.0)::Nothing
    oper = rand.oper

    operQ, Qatlevel = oper.operQ, oper.Qatlevel
    t, packer, buff = oper.ptxt_modulus, oper.packer, oper.tensor_buff[1].coeffs[1]

    # Pack the input message.
    if !ismissing(packer) && ispacking
        pack_to!(buff, msg, packer)
    else
        if length(msg) ≠ length(buff)
            throw(DomainError("The length of the plaintext should match the ring size."))
        end
        copy!(buff, msg)
    end

    # Randomise the input message.
    buff_rand = rand.buff_rand
    @. buff_rand = buff
    if typeof(dgs) == TwinCDTSampler
        @inbounds for i = eachindex(buff)
            I = sample((t.Q - buff[i]) / t.Q, dgs)
            buff_rand[i] += I * Int64(t.Q)
        end
    else
        @inbounds for i = eachindex(buff)
            I = sample(dgs[t.Q-buff[i]])
            buff_rand[i] += I * Int64(t.Q)
        end
    end
    
    if isPQ
        if ismissing(operQ.evalP)
            throw(ErrorException("The parameter does not support PQ."))
        end

        Qlen = Qatlevel[end]
        Plen = length(operQ.evalP)
        evalPQ = geteval_at(Plen + Qlen, operQ, isPQ=true)

        resize!(res, Plen + Qlen)
        for i = 1:Plen+Qlen
            Bred_to!(res.val.coeffs[i], buff_rand, evalPQ[i].Q)
        end
        res.val.isntt[] = false
        ntt!(res.val, evalPQ)

        res.auxQ[] = 0
        res.isPQ[] = true
    else
        Qlen = Qatlevel[end]
        evalQ = geteval_at(Qlen, operQ)

        resize!(res, Qlen)
        for i = 1:Qlen
            Bred_to!(res.val.coeffs[i], buff_rand, evalQ[i].Q)
        end
        res.val.isntt[] = false
        ntt!(res.val, evalQ)

        res.auxQ[] = 0
        res.isPQ[] = false
    end

    return nothing
end

decode(x::PlainPoly, rand::RandBFVOperator; ispacking::Bool=true)::Vector{UInt64} = decode(x, rand.oper, ispacking=ispacking)

decode_to!(res::Vector{UInt64}, x::PlainPoly, rand::RandBFVOperator; ispacking::Bool=true)::Nothing = decode_to!(res, x, rand.oper, ispacking=ispacking)

#======================================================================================#
################################## ENCRYPT AND DECRYPT #################################
#======================================================================================#

encrypt(msg::PlainText, entor::Encryptor, oper::RandBFVOperator)::BFV = encrypt(msg, entor, oper.oper)

encrypt_to!(res::BFV, msg::PlainText, entor::Encryptor, oper::RandBFVOperator)::Nothing = encrypt_to!(res, msg, entor, oper.oper)

decrypt(x::BFV, entor::Encryptor, oper::RandBFVOperator)::PlainPoly = decrypt(x, entor, oper.oper)

decrypt_to!(res::PlainPoly, x::BFV, entor::Encryptor, oper::RandBFVOperator)::Nothing = decrypt_to!(res, x, entor, oper.oper)

get_error(ct::BFV, entor::Encryptor, oper::RandBFVOperator)::Vector{BigInt} = get_error(ct, entor, oper.oper)

#======================================================================================#
################################## PLAINTEXT OPERATIONS ################################
#======================================================================================#

change_level(x::PlainText, targetlvl::Integer, rand::RandBFVOperator)::PlainText = change_level(x, targetlvl, rand.oper)

change_level_to!(res::T, x::T, targetlvl::Integer, rand::RandBFVOperator) where {T<:PlainText} = begin
    change_level_to!(res, x, targetlvl, rand.oper)
end::Nothing

#================================================================================================#
################################## CIPHERTEXT OPERATIONS ########################################
#================================================================================================#

drop_level_to!(res::BFV, x::BFV, targetlvl::Integer, rand::RandBFVOperator)::Nothing = drop_level_to!(res, x, targetlvl, rand.oper)

rescale(x::BFV, rand::RandBFVOperator)::BFV = rescale(x, rand.oper)

rescale_to!(res::BFV, x::BFV, rand::RandBFVOperator)::Nothing = begin
    rescale_to!(res, x, rand.oper)
    return nothing
end

neg(x::BFV, rand::RandBFVOperator)::BFV = neg(x, rand.oper)

neg_to!(res::BFV, x::BFV, rand::RandBFVOperator)::Nothing = begin
    neg_to!(res, x, rand.oper)
    return nothing
end

add(x::BFV, y::PlainPoly, rand::RandBFVOperator)::BFV = add(x, y, rand.oper)

add(x::PlainPoly, y::BFV, rand::RandBFVOperator)::BFV = add(x, y, rand.oper)

add(x::PlainConst, y::BFV, rand::RandBFVOperator)::BFV = add(x, y, rand.oper)

add(x::BFV, y::PlainConst, rand::RandBFVOperator)::BFV = add(x, y, rand.oper)

add(x::BFV, y::BFV, rand::RandBFVOperator)::BFV = add(x, y, rand.oper)

add_to!(res::BFV, x::BFV, y::PlainPoly, rand::RandBFVOperator)::Nothing = begin
    add_to!(res, x, y, rand.oper)
    return nothing
end

add_to!(res::BFV, x::PlainPoly, y::BFV, rand::RandBFVOperator)::Nothing = begin
    add_to!(res, x, y, rand.oper)
    return nothing
end

add_to!(res::BFV, x::PlainConst, y::BFV, rand::RandBFVOperator)::Nothing = begin
    add_to!(res, x, y, rand.oper)
    return nothing
end

add_to!(res::BFV, x::BFV, y::PlainConst, rand::RandBFVOperator)::Nothing = begin
    add_to!(res, x, y, rand.oper)
    return nothing
end

add_to!(res::BFV, x::BFV, y::BFV, rand::RandBFVOperator)::Nothing = begin
    add_to!(res, x, y, rand.oper)
    return nothing
end

sub(x::BFV, y::PlainPoly, rand::RandBFVOperator)::BFV = sub(x, y, rand.oper)

sub(x::PlainPoly, y::BFV, rand::RandBFVOperator)::BFV = sub(x, y, rand.oper)

sub(x::BFV, y::PlainConst, rand::RandBFVOperator)::BFV = sub(x, y, rand.oper)

sub(x::PlainConst, y::BFV, rand::RandBFVOperator)::BFV = sub(x, y, rand.oper)

sub(x::BFV, y::BFV, rand::RandBFVOperator)::BFV = sub(x, y, rand.oper)

sub_to!(res::BFV, x::BFV, y::PlainPoly, rand::RandBFVOperator)::Nothing = begin
    sub_to!(res, x, y, rand.oper)
    return nothing
end

sub_to!(res::BFV, x::PlainPoly, y::BFV, rand::RandBFVOperator)::Nothing = begin
    sub_to!(res, x, y, rand.oper)
    return nothing
end

sub_to!(res::BFV, x::BFV, y::PlainConst, rand::RandBFVOperator)::Nothing = begin
    sub_to!(res, x, y, rand.oper)
    return nothing
end

sub_to!(res::BFV, x::PlainConst, y::BFV, rand::RandBFVOperator)::Nothing = begin
    sub_to!(res, y, x, rand.oper)
    return nothing
end

sub_to!(res::BFV, x::BFV, y::BFV, rand::RandBFVOperator)::Nothing = begin
    sub_to!(res, x, y, rand.oper)
    return nothing
end

mul(x::BFV, y::PlainPoly, τ::Float64, rand::RandBFVOperator; islazy::Bool=false)::BFV = begin
    res = similar(y)
    mul_to!(res, x, y, τ, rand, islazy=islazy)
    res
end

function mul_to!(res::BFV, x::BFV, y::PlainPoly, τ::Float64, rand::RandBFVOperator; islazy::Bool=false)::Nothing
    operQ = rand.oper.operQ

    @assert length(res.val) == length(x.val) "The input and output ciphertext length should match."
    @assert x.val.b.isntt[] && x.val.a.isntt[] "The input ciphertext should be in NTT domain."

    evalQ, Qlen = operQ.evalQ, length(x.val)
    mul_to!(res, x, y, rand.oper)

    # Add noise.
    buffQ = operQ.tensor_buff[1][1:Qlen]
    buffQ.isntt[] = false

    @inbounds for j = 1:buffQ.N
        ej = sample(rand.rgs, τ)
        @simd for i = eachindex(evalQ)
            buffQ.coeffs[i][j] = Bred(ej, evalQ[i])
        end
    end
    ntt!(buffQ, evalQ)
    add_to!(res.val.b, res.val.b, buffQ, evalQ)

    return nothing
end

mul(x::PlainPoly, y::BFV, τ::Float64, rand::RandBFVOperator; islazy::Bool=false)::BFV = mul(y, x, τ, rand, islazy=islazy)

mul_to!(res::BFV, x::PlainPoly, y::BFV, τ::Float64, rand::RandBFVOperator; islazy::Bool=false)::Nothing = mul_to!(res, y, x, τ, rand, islazy=islazy)

mul(x::BFV, y::BFV, rlk::RLEV, dgs::Randomiser, τx::Float64, τy::Float64, rand::RandBFVOperator; islazy::Bool=false)::BFV = begin
    res = similar(x)
    mul_to!(res, x, y, rlk, dgs, τx, τy, rand, islazy=islazy)
    res
end

#Cmul
function mul_to!(res::BFV, x::BFV, y::BFV, rlk::RLEV, dgs::Randomiser, τx::Float64, τy::Float64, rand::RandBFVOperator; islazy::Bool=false)::Nothing
    oper, evalT, us = rand.oper, rand.evalT, rand.us
    t = evalT.Q

    # sample r1, r2 
    p1, p2 = rand.ptxt_buff[1], rand.ptxt_buff[2]
    r1, r2 = rand.poly_buff[1], rand.poly_buff[2]
    uniform_random_to!(us, r1, t)
    uniform_random_to!(us, r2, t)
    encode_to!(p1, r1, dgs, rand, isPQ=false, ispacking=false)
    encode_to!(p2, r2, dgs, rand, isPQ=false, ispacking=false)

    # Compute res = (x + r1)(y + r2)
    xlen, ylen = length(x.val), length(y.val)
    tmpx, tmpy = rand.ct_buff[1][1:xlen], rand.ct_buff[2][1:ylen]
    add_to!(tmpx, x, p1, oper)
    add_to!(tmpy, y, p2, oper)

    # Randomize
    buff = rand.oper.ct_buff[end].val[1:length(tmpx)]
    rlwe_sample_to!(buff, rand.entor)
    add_to!(tmpx.val, tmpx.val, buff, oper.operQ)

    buff = rand.oper.ct_buff[end].val[1:length(tmpy)]
    rlwe_sample_to!(buff, rand.entor)
    add_to!(tmpy.val, tmpy.val, buff, oper.operQ)

    mul_to!(res, tmpx, tmpy, rlk, oper)

    # Compute res -= r1 * (y + r2) + r2 * (x + r1)
    mul_to!(tmpy, p1, tmpy, τy, rand)
    sub_to!(res, res, tmpy, oper)

    mul_to!(tmpx, p2, tmpx, τx, rand)
    sub_to!(res, res, tmpx, oper)

    # Compute res += r1 * r2
    mul_to!(p1, p1, p2, oper.operQ)
    add_to!(res, res, p1, oper)

    return nothing
end

decompose_a(x::BFV, rand::RandBFVOperator)::Tensor = decompose_a(x, rand.oper)

decompose_a_to!(res::Tensor, x::BFV, rand::RandBFVOperator)::Nothing = decompose_a_to!(res, x, rand.oper)

mask_randomize(x::BFV, rand::RandBFVOperator)::BFV = begin
    res = similar(x)
    mask_randomize!(res, rand)
    res
end

mask_randomize!(x::BFV, rand::RandBFVOperator)::Nothing = mask_randomize_to!(x, x, rand)

mask_randomize_to!(res::BFV, x::BFV, rand::RandBFVOperator)::Nothing = begin
    buff = rand.ct_buff[2].val[1:length(x.val)]

    rlwe_sample_to!(buff, rand.entor)
    add_to!(res.val, buff, x.val, rand.oper.operQ)

    return nothing
end

keyswitch(x::BFV, ksk::RLEV, rand::RandBFVOperator)::BFV = begin
    res = similar(x)
    keyswitch_to!(res, x, ksk, rand)
    res
end

keyswitch_to!(res::BFV, x::BFV, ksk::RLEV, rand::RandBFVOperator)::Nothing = begin
    mask_randomize_to!(res, x, rand)
    keyswitch_to!(res, res, ksk, rand.oper)
    return nothing
end

rotate(x::BFV, idx::NTuple{N,Int64}, rtk::RLEV, rand::RandBFVOperator) where {N} = begin
    res = similar(x)
    rotate_to!(res, x, idx, rtk, rand)
    res
end::BFV

function rotate_to!(res::BFV, x::BFV, idx::NTuple{N,Int64}, rtk::RLEV, rand::RandBFVOperator) where {N}
    packer = rand.oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    automorphism_to!(res, x, autidx, rtk, rand)
end

hoisted_rotate(adec::Tensor, x::BFV, idx::NTuple{N,Int64}, rtk::RLEV, rand::RandBFVOperator) where {N} = begin
    res = similar(x)
    hoisted_rotate_to!(res, adec, x, idx, rtk, rand)
    res
end::BFV

function hoisted_rotate_to!(res::BFV, adec::Tensor, x::BFV, idx::NTuple{N,Int64}, rtk::RLEV, rand::RandBFVOperator)::Nothing where {N}
    packer = rand.oper.packer
    if ismissing(packer)
        throw(ErrorException("Rotation operation cannot be defined without SIMD packing."))
    end

    cube, cubegen = packer.cube, packer.cubegen
    if length(cube) ≠ N
        throw(DimensionMismatch("The number of indices should match the number of dimensions."))
    end

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    hoisted_automorphism_to!(res, adec, x, autidx, rtk, rand)

    return nothing
end

automorphism(x::BFV, idx::Int64, atk::RLEV, rand::RandBFVOperator) = begin
    res = deepcopy(x)
    automorphism_to!(res, x, idx, atk, rand)
    res
end

automorphism_to!(res::BFV, x::BFV, idx::Int64, atk::RLEV, rand::RandBFVOperator) = begin
    mask_randomize_to!(res, x, rand)
    automorphism_to!(res, res, idx, atk, rand.oper)
end

hoisted_automorphism(adec::Tensor, x::BFV, idx::Int64, atk::RLEV, rand::RandBFVOperator) = begin
    res = deepcopy(x)
    hoisted_automorphism_to!(res, adec, x, idx, atk, rand)
    res
end

hoisted_automorphism_to!(res::BFV, adec::Tensor, x::BFV, idx::Int64, atk::RLEV, rand::RandBFVOperator)::Nothing = begin
    hoisted_automorphism_to!(res, adec, x, idx, atk, rand.oper)
    return nothing
end
