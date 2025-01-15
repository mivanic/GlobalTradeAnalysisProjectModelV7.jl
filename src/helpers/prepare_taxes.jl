function prepare_taxes(; data=data, hData=hData)

    to = hData["makb"] ./ hData["maks"]
    to[hData["maks"].==0] .= 1
    tfe = hData["evfp"] ./ hData["evfb"]
    tfe[hData["evfb"].==0] .= 1
    txs = hData["vfob"] ./ hData["vxsb"]
    tx = NamedArray(ones(size(txs)[1:2]), names(txs)[1:2])
    tms = hData["vmsb"] ./ hData["vcif"]
    tm = NamedArray(ones(size(tms)[1:2]), names(tms)[1:2])
    tfd = hData["vdfp"] ./ hData["vdfb"]
    tfm = hData["vmfp"] ./ hData["vmfb"]
    tpd = hData["vdpp"] ./ hData["vdpb"]
    tpm = hData["vmpp"] ./ hData["vmpb"]
    tgd = hData["vdgp"] ./ hData["vdgb"]
    tgm = hData["vmgp"] ./ hData["vmgb"]
    tid = hData["vdip"] ./ hData["vdib"]
    tim = hData["vmip"] ./ hData["vmib"]
    tinc = hData["evfb"] ./ hData["evos"]
    tinc[isnan.(tinc)] .= 1

    taxes = Dict(
        :to => to,
        :tfe => tfe,
        :txs => txs,
        :tx => tx,
        :tms => tms,
        :tm => tm,
        :tfd => tfd,
        :tfm => tfm,
        :tpd => tpd,
        :tpm => tpm,
        :tgd => tgd,
        :tgm => tgm,
        :tid => tid,
        :tim => tim,
        :tinc => tinc
    )

    return (data = merge(data, Dict(String(k) => taxes[k] for k ∈ keys(taxes))))
end