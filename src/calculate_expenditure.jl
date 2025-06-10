"""
        calculate_expenditure(; sets, data0, data1, parameters, max_iter=50, constr_viol_tol=1e-5, bound_push=1e-15)

Calculate the expenditure at post simulation utility (in data1) at original prices (data0)

### Args:

- `sets`: a dictionary of sets
- `data0`: a dictionary of data from baseline solution
- `data1`: a dictionary of data from simulation solution
- `parameters`: a dictionary of parameters

### Optional args:

- `max_iter`: maximum number of iterations
- `constr_viol_tol`: accuracy for constraint satisfaction
- `bound_push`: mandatory move of the starting values from constraint bounds

### Returns:

A vector of expenditures.
"""
function calculate_expenditure(; sets, data0, data1, parameters, max_iter=50, constr_viol_tol=1e-5, bound_push=1e-15)

    # Read  sets
    (; reg, comm) = NamedTuple(Dict(Symbol(k) => sets[k] for k ∈ keys(sets)))

    # Read hard parameters
    (; subpar, incpar) = NamedTuple(Dict(Symbol(k) => parameters[k] for k ∈ keys(parameters)))

    # Set up the model
    model = JuMP.Model(Ipopt.Optimizer)

    # Set up the general constraints
    p_min = 1e-8
    p_max = 1e+8
    q_min = 1e-9
    q_max = 1e+12
    y_min = 1e-8
    y_max = 1e+12

    # All variables used in the model
    @variables(model,
        begin
            # Population
            q_min <= pop[reg] <= q_max

            # Income
            y_min <= y[reg] <= y_max

            q_min <= u[reg] <= q_max
            q_min <= ug[reg] <= q_max
            #q_min <= us[reg] <= q_max

            # Private consumption        
            y_min <= yp[reg] <= y_max
            q_min <= up[reg] <= q_max
            0 <= uelas[reg] <= 10
            0 <= uepriv[reg] <= 10
            p_min <= ppa[comm, reg] <= p_max
            q_min <= qpa[comm, reg] <= q_max

            # Government consumption
            y_min <= yg[reg] <= y_max
            p_min <= pgov[reg] <= p_max
            p_min <= pga[comm, reg] <= p_max
            q_min <= qga[comm, reg] <= q_max

            # Saving
            -q_max <= qsave[reg] <= q_max
            p_min <= psave[reg] <= p_max

            # Parameters
            1e-8 <= σyp[reg] <= 1
            1e-8 <= σyg[reg] <= 1
            1e-8 <= β_qpa[comm, reg]
            1e-8 <= α_qga[comm, reg] <= 1


            1e-2 <= p[reg] <=1e2
            1e-2 <= ppriv[reg] <=1e2
        end
    )



    # All model equations
    @constraints(model,
        begin
            # Utility
            #e_u, u .== up.^σyp .* ug.^σyg .* us.^(1 .-σyp .- σyg)
            e_u, log.(u) .== log.((y ./ p ./ pop))

            e_ug, ug .== yg ./ pop ./ pgov
            #e_us, us .== qsave ./ pop
            e_uepriv[r=reg], uepriv[r] == sum(qpa[:, r] .* ppa[:, r] .* Vector(incpar[:, r])) / yp[r]
            e_uelas[r=reg], uelas[r] == 1 / (σyp[r] / uepriv[r] + σyg[r] + (1 - σyp[r] - σyg[r]))

            # Household Income
            e_yp, log.(yp) .== log.(y .* Vector(σyp) .* uelas ./ uepriv)

            # Household consumption
            e_qpa[r=reg], log.([Vector(qpa[:, r] ./ pop[r]); 1]) .== log.(cde(Vector(1 .- subpar[:, r]), Vector(β_qpa[:, r]), Vector(incpar[:, r]), up[r], Vector(ppa[:, r]), yp[r] / pop[r]))

            # Government Income
            e_yg, log.(yg) .== log.(y .* Vector(σyg) .* uelas)

            # Government expenditure
            e_qga[r=reg], log.(pga[:, r] .* qga[:, r]) .== log.(yg[r] .* Vector(α_qga[:, r])) ##This one
            e_pgov[r=reg], log.(pgov[r] * sum(qga[:, r])) == log.(sum(qga[:, r] .* pga[:, r]))

            # Saving
            e_qsave, log.(y) .== log.(yp .+ yg .+ psave .* qsave .* uelas)

            e_p, log.(p) .== log.(ppriv .* σyp .+ pgov .* σyg .+ psave .* (1 .- σyp .- σyg))
            e_ppriv[r=reg], log(ppriv[r] * sum(qpa[:, r])) == log(sum(ppa[:, r] .* qpa[:, r]))
        end
    )

    # Set starting values
    for k in keys(data0)
        if Symbol(k) ∈ keys(object_dictionary(model))
            if data0[k] isa NamedArray
                set_start_value.(model[Symbol(k)], Array(data0[k]))
            else
                set_start_value.(model[Symbol(k)], data0[k])
            end
        end
    end

    for vv ∈ ["ppa","pga","psave"]
        fix.(Array(model[Symbol(vv)]), data0[vv]; force = true)
    end

    for vv ∈ ["u","σyp","σyg","β_qpa","α_qga","pop"]
        fix.(Array(model[Symbol(vv)]), data1[vv]; force = true)
    end
    set_attribute(model, "max_iter", max_iter)
    set_attribute(model, "constr_viol_tol", constr_viol_tol)
    set_attribute(model, "bound_push", bound_push)

    # Solve
    optimize!(model)

    # Save results
    results = merge(Dict(
            String(k) => begin
                arrayOut = NamedArray(zeros(map(length, v.axes)), v.axes)
                arrayOut[is_valid.(model, v).data] .= value.(Array(v)[is_valid.(model, v).data])
                arrayOut[.!is_valid.(model, v).data] .= NaN
                arrayOut
            end for (k, v) in object_dictionary(model)
            if v isa AbstractArray{VariableRef}
        ), Dict(
            String(k) => begin
                (is_valid(model, v) ? value.(v) : NaN)
            end for (k, v) in object_dictionary(model)
            if v isa VariableRef
        ))

    return results["y"]
end
