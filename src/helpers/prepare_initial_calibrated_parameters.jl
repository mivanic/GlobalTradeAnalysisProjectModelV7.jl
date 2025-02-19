function prepare_initial_calibrated_parameters(; data, sets, parameters, hData)

    # Read all parameters
    (; esubt, esubc, esubva, esubd, etraq, esubq, subpar, incpar, esubg, etrae, esubm, esubs, endowflag) = NamedTuple(Dict(Symbol(k) => parameters[k] for k ∈ keys(parameters)))

    # Read sets
    (; comm, reg, marg, endws, endwc, endwm, endwms) = NamedTuple(Dict(Symbol(k) => sets[k] for k ∈ keys(sets)))


    # Read all data
    (; pint, pva, qint, qva, qo, pfa, qfa, pfe, qfe, pfm, qfm, pfd, qfd, qca, ps, pca, qc, y, yp, yg, ppa, qpa, qpd, qpm, ppd, ppm, pgd, qgd, pgm, qgm, pga, qga, pgov, pid, qid, pim, qim, pia, qia, pinv, qinv, globalcgds, qxs, pmds, qms, qtmfsd, qtm, pds, qst, pes, qes, qe, up, pop, kb) = NamedTuple(Dict(Symbol(k) => data[k] for k ∈ keys(data)))

    # Helper function
    function ρ(σ)
        (σ .- 1) ./ σ
    end

    # The main work function
    function ces_parameters(sigma, sigma_expanded, prices, quantities, output)

        θ = prices .* quantities .^ (1 ./ sigma_expanded)
        θₐ = prices .* quantities
        θₑ = quantities

        θ[isinf.(θ)] .= 0
        θ[isnan.(θ)] .= 0


        inner = Int.(collect(size(θ) ./ size(θ)))
        inner[1] = size(θ, 1)

        neworder = 1:length(size(quantities))
        neworder = [neworder[length(neworder)]; neworder[1:(length(neworder)-1)]]

        α = θ ./ repeat(mapslices(sum, θ, dims=1), inner=inner)
        αₐ = θₐ ./ repeat(mapslices(sum, θₐ, dims=1), inner=inner)
        αₑ = quantities ./ permutedims(repeat(output, inner=[Int.(collect(size(output) ./ size(output))); size(θ, 1)]), neworder)   #θₑ ./ repeat(mapslices(sum, θₑ, dims=1), inner=inner)

        α[sigma_expanded.==1] = αₐ[sigma_expanded.==1]
        α[sigma_expanded.==0] = αₑ[sigma_expanded.==0]

        γ_temp = (α .* quantities .^ ρ(sigma_expanded))
        γ_temp[isnan.(γ_temp)] .= 0

        non_zero_αₑ = repeat(mapslices(f -> length(f[f.!=0]), αₑ, dims=1), inner=inner)

        γₑ_temp = quantities ./ αₑ ./ non_zero_αₑ
        γₑ_temp[isnan.(γₑ_temp)] .= 0

        γ = output .* selectdim(mapslices(sum, γ_temp, dims=1), 1, 1) .^ (-1 ./ ρ(sigma))
        γₐ = output .* selectdim(mapslices(prod, quantities .^ αₐ, dims=1), 1, 1) .^ (-1)
        γₑ = output .* selectdim(mapslices(sum, γₑ_temp, dims=1), 1, 1) .^ (-1)

        γ[sigma.==1] = γₐ[sigma.==1]
        γ[sigma.==0] = γₑ[sigma.==0]
        return (α, γ)
    end

    sigma = copy(esubt)
    sigma_expanded = permutedims(cat(sigma, sigma, dims=3), [3, 1, 2])
    prices = permutedims(cat(pint, pva, dims=3), [3, 1, 2])
    quantities = permutedims(cat(qint, qva, dims=3), [3, 1, 2])
    output = qo

    (α_qintva, γ_qintva) = ces_parameters(sigma, sigma_expanded, prices, quantities, output)


    prices = pfa
    quantities = qfa
    sigma = copy(esubc)
    sigma_expanded = permutedims(repeat(sigma, inner=[1, 1, size(pfa, 1)]), [3, 1, 2])
    output = qint

    (α_qfa, γ_qfa) = ces_parameters(sigma, sigma_expanded, prices, quantities, output)

    prices = pfe
    quantities = qfe
    sigma = copy(esubva)
    sigma_expanded = permutedims(repeat(sigma, inner=[1, 1, size(pfe, 1)]), [3, 1, 2])
    output = qva

    (α_qfe, γ_qfe) = ces_parameters(sigma, sigma_expanded, prices, quantities, output)


    prices = permutedims(cat(pfd, pfm, dims=4), [4, 1, 2, 3])
    quantities = permutedims(cat(qfd, qfm, dims=4), [4, 1, 2, 3])
    sigma = permutedims(repeat(esubd, inner=[1, 1, size(pfd, 2)]), [1, 3, 2])
    sigma_expanded = permutedims(repeat(sigma, inner=[1, 1, 1, 2]), [4, 1, 2, 3])
    output = qfa

    (α_qfdqfm, γ_qfdqfm) = ces_parameters(sigma, sigma_expanded, prices, quantities, output)

    prices = ps
    quantities = qca
    sigma = etraq
    sigma_expanded = permutedims(repeat(sigma, inner=[1, 1, size(quantities, 1)]), [3, 1, 2])
    output = qo

    (α_qca, γ_qca) = ces_parameters(sigma, sigma_expanded, prices, quantities, output)


    prices = pca
    quantities = qca
    sigma = 1 ./ esubq
    sigma_expanded = permutedims(repeat(sigma, inner=[1, 1, size(quantities, 1)]), [3, 1, 2])
    output = qc

    (α_pca, γ_pca) = ces_parameters(sigma, sigma_expanded, prices, quantities, output)

    σyp = yp ./ y
    σyg = yg ./ y

    yp_expanded = permutedims(repeat(yp, inner=[1, size(qpd, 1)]), [2, 1])
    β_qpa_help = (qpd .+ qpm) ./ (subpar ./ (yp_expanded .^ subpar))
    β_qpa = β_qpa_help ./ Array(repeat(mapslices(sum, β_qpa_help, dims=1), inner=[size(β_qpa_help, 1), 1]))

    cy = NamedArray(mapslices(sum, qpa .* ppa, dims=1)[1, :], names(qpa, 2))


    m = JuMP.Model(Ipopt.Optimizer)
    @variables(m,
        begin
            1e-6 <= β2[comm]
            1e-6 <= u3
            pop2
            qpa2
            cy2
            subpar2[comm]
            incpar2[comm]
        end
    )

    @constraints(m,
        begin
            c, log.([Vector(qpa2 / pop2); 1]) .== log.(cde(Vector(1 .- subpar2), Vector(β2), Vector(incpar2), u3, Vector(ppa2), cy2 / pop2))
        end
    )
    u2 = NamedArray(ones(length(reg)),reg)

    for r ∈ reg
        set_start_value(u3, 1000)
        set_start_value.(β2, 1)
        fix.(cy2,cy[r])
        fix.(pop2,pop[r])
        fix.(qpa2,qpa[:,r])
        fix.(subpar2,subpar[:,r])
        fix.(incpar2,incpar[:,r])
        optimize!(m)
        u2[r].=value(u3)
        if is_solved_and_feasible(m)
            throw("Utility not found")
        end
    end

    # Household domestic/imported sourcing

    prices = permutedims(cat(ppd, ppm, dims=3), [3, 1, 2])
    quantities = permutedims(cat(qpd, qpm, dims=3), [3, 1, 2])
    sigma = esubd
    sigma_expanded = permutedims(repeat(sigma, inner=[1, 1, 2]), [3, 1, 2])
    output = qpa

    (α_qpdqpm, γ_qpdqpm) = ces_parameters(sigma, sigma_expanded, prices, quantities, output)

    # Government consumption
    prices = pga
    quantities = qga
    sigma = esubg
    sigma_expanded = permutedims(repeat(sigma, inner=[1, size(prices, 1)]), [2, 1])
    output = yg ./ pgov

    (α_qga, γ_qga) = ces_parameters(sigma, sigma_expanded, prices, quantities, output)

    # Government consumption sourcing
    prices = permutedims(cat(pgd, pgm, dims=3), [3, 1, 2])
    quantities = permutedims(cat(qgd, qgm, dims=3), [3, 1, 2])
    sigma = esubd
    sigma_expanded = permutedims(repeat(sigma, inner=[1, 1, 2]), [3, 1, 2])
    output = qga

    (α_qgdqgm, γ_qgdqgm) = ces_parameters(sigma, sigma_expanded, prices, quantities, output)


    # Investment
    prices = pia
    quantities = qia
    sigma = esubg .* 0
    sigma_expanded = permutedims(repeat(sigma, inner=[1, size(prices, 1)]), [2, 1])
    output = qinv
    (α_qia, γ_qia) = ces_parameters(sigma, sigma_expanded, prices, quantities, output)

    γ_qia = NamedArray(γ_qia, (reg))

    # Government consumption sourcing
    prices = permutedims(cat(pid, pim, dims=3), [3, 1, 2])
    quantities = permutedims(cat(qid, qim, dims=3), [3, 1, 2])
    sigma = esubd
    sigma_expanded = permutedims(repeat(sigma, inner=[1, 1, 2]), [3, 1, 2])
    output = qia

    (α_qidqim, γ_qidqim) = ces_parameters(sigma, sigma_expanded, prices, quantities, output)

    # Imports
    prices = permutedims(pmds, [2, 1, 3])
    quantities = permutedims(qxs, [2, 1, 3])
    sigma = esubm
    sigma_expanded = permutedims(repeat(sigma, inner=[1, 1, size(prices, 1)]), [3, 1, 2])
    output = qms
    (α_qxs_, γ_qxs) = ces_parameters(sigma, sigma_expanded, prices, quantities, output)

    α_qxs = permutedims(α_qxs_, [2, 1, 3])


    # Transportation margins
    α_qtmfsd = qtmfsd ./ permutedims(repeat(qxs, inner=[1, 1, 1, size(qtmfsd, 1)]), [4, 1, 2, 3])

    # Margin distribution parameters
    prices = permutedims(pds, [2, 1])[:, names(qst)[1]]
    quantities = permutedims(qst, [2, 1])
    sigma = esubs
    sigma_expanded = permutedims(repeat(sigma, inner=[1, size(prices, 1)]), [2, 1])
    output = qtm
    (α_qst_, γ_qst) = ces_parameters(sigma, sigma_expanded, prices, quantities, output)

    α_qst = permutedims(α_qst_, [2, 1])
    γ_qst = NamedArray(γ_qst, (marg))

    # Factor allocation
    prices = permutedims(pes[endwms, :, :], [2, 1, 3])
    quantities = permutedims(qes[endwms, :, :], [2, 1, 3])
    sigma = etrae[endwms, :]
    sigma_expanded = permutedims(repeat(sigma, inner=[1, 1, size(prices, 1)]), [3, 1, 2])
    output = qe[endwms, :]

    (α_qes2_, γ_qes2) = ces_parameters(sigma, sigma_expanded, prices, quantities, output)

    α_qes2 = permutedims(α_qes2_, [2, 1, 3])


    # Vector(qes["land", :, "eu"])[Vector(α_qes2["land", :, "eu"]).!=0] 
    # Vector(demand_ces(qe["land", "eu"], Vector(pes["land", :, "eu"])[Vector(α_qes2["land", :, "eu"]).!=0], Vector(α_qes2["land", :, "eu"])[Vector(α_qes2["land", :, "eu"]).!=0], etrae["land", "eu"], γ_qes2[3, 2]))

    

    δ = hData["vdep"] ./ hData["vkb"]
    ρ = mapslices(sum, hData["evos"][endwc, :, :], dims=[1, 2])[1, 1, :] ./ hData["vkb"]

    #α_qinv = (qinv .- δ .* kb) ./ fill(globalcgds, size(qinv, 1))
    α_qinv = (hData["save"]) ./ fill(sum(hData["save"]), size(qinv, 1))


    # ϵs
    ϵ_qxs = copy(γ_qxs)    
    ϵ_qxs[] .= 1

    ϵ_qfe = copy(γ_qfe)    
    ϵ_qfe[] .= 1

    ϵ_qes2 = copy(γ_qes2)    
    ϵ_qes2[] .= 1

    ϵ_qfdqfm = copy(γ_qfdqfm)    
    ϵ_qfdqfm[] .= 1

    ϵ_qpdqpm = copy(γ_qpdqpm)    
    ϵ_qpdqpm[] .= 1

    ϵ_qgdqgm = copy(γ_qgdqgm)    
    ϵ_qgdqgm[] .= 1

    ϵ_qidqim = copy(γ_qidqim)    
    ϵ_qidqim[] .= 1


    # Return the new parameter values along with the old ones
    new_parameters = Dict(
        :α_qintva => α_qintva,
        :γ_qintva => γ_qintva,
        :α_qfa => α_qfa,
        :γ_qfa => γ_qfa,
        :α_qfe => α_qfe,
        :γ_qfe => γ_qfe,
        :ϵ_qfe => ϵ_qfe,
        :α_qfdqfm => α_qfdqfm,
        :γ_qfdqfm => γ_qfdqfm,
        :ϵ_qfdqfm => ϵ_qfdqfm,
        :α_qca => α_qca,
        :γ_qca => γ_qca,
        :α_pca => α_pca,
        :γ_pca => γ_pca,
        :σyp => σyp,
        :σyg => NamedArray(σyg, reg),
        :β_qpa => NamedArray(Array(value.(β)), axes(β)),
        :α_qpdqpm => α_qpdqpm,
        :γ_qpdqpm => γ_qpdqpm,
        :ϵ_qpdqpm => ϵ_qpdqpm,
        :α_qga => α_qga,
        :γ_qga => γ_qga,
        :α_qgdqgm => α_qgdqgm,
        :γ_qgdqgm => γ_qgdqgm,
        :ϵ_qgdqgm => ϵ_qgdqgm,
        :α_qia => α_qia,
        :γ_qia => γ_qia,
        :α_qidqim => α_qidqim,
        :γ_qidqim => γ_qidqim,
        :ϵ_qidqim => ϵ_qidqim,
        :α_qxs => α_qxs,
        :γ_qxs => γ_qxs,
        :ϵ_qxs => ϵ_qxs,
        :α_qtmfsd => α_qtmfsd,
        :α_qst => α_qst,
        :γ_qst => γ_qst,
        :α_qes2 => α_qes2[endws, :, :],
        :γ_qes2 => γ_qes2[endws, :],
        :ϵ_qes2 => ϵ_qes2[endws, :],
        :α_qinv => α_qinv,
        :δ => δ,
        :ρ => ρ
    )

    new_data = Dict(:up => NamedArray(u2, reg))

    return (
        parameters=parameters,
        data=merge(
            data,
            Dict(String(k) => new_data[k] for k ∈ keys(new_data)),
            Dict(String(k) => new_parameters[k] for k ∈ keys(new_parameters)))
    )
end