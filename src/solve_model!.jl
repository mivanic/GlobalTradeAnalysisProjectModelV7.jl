"""
    solve_model!(model_container; max_iter, constr_viol_tol, bound_push)
"""
function solve_model!(model_container; max_iter, constr_viol_tol, bound_push)

    # Specify the solver parameters
    set_attribute(model_container.model, "max_iter", max_iter)
    set_attribute(model_container.model, "constr_viol_tol", constr_viol_tol)
    set_attribute(model_container.model, "bound_push", bound_push)

    # Solve the model
    optimize!(model_container.model)
end
