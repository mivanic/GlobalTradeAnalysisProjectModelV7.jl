"""
    run_model!(model_container; max_iter=50, constr_viol_tol=1e-8, bound_push=1e-15)

Runs the model (model_container.model) by executing these steps: (1) set all starting 
values to those found in the data (model_container.data), (2) fix all variables that
are specified as exogenous (model_container.fixed), (3) solve the model, (4) load all
variable values to the data (model_container.data)

### Args:
- `model_container`: model container

### Optional args:
- `max_iter`: maximum number of iterations
- `constr_viol_tol`: accuracy for constraint satisfaction
- `bound_push`: mandatory move of the starting values from constraint bounds

"""
function run_model!( model_container; max_iter=50, constr_viol_tol=1e-8, bound_push=1e-15)

    # initialize the model
    initialize_model!(model_container)

    # solve_model
    solve_model!(model_container; max_iter=max_iter, constr_viol_tol=constr_viol_tol, bound_push=bound_push)

    # Exttract values
    extract_model!(model_container)

end