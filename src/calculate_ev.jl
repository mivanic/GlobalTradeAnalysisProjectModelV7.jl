"""
    calculate_ev(; sets, data0, data1, parameters, max_iter=50, constr_viol_tol=1e-5, bound_push=1e-15)

Doc string
"""
function calculate_ev(; sets, data0, data1, parameters, max_iter=50, constr_viol_tol=1e-5, bound_push=1e-15)
    printstyled("Notice: calculate_ev was replaced with calculate_expenditure. Please use that.\n", color = :yellow)
    return calculate_expenditure(; sets=sets, data0=data0, data1=data1, parameters=parameters, max_iter=max_iter, constr_viol_tol=constr_viol_tol, bound_push=bound_push)
end
