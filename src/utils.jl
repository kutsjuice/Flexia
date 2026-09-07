# This file contains different useful function

function set_initial_state_value!(state, sys, body, values)
    dofs = get_body_position_dofs(sys, body)
    state[dofs] .= values
end

function set_pos_state_value!(state, sys, body, values)
    dofs = get_body_generalized_dofs(sys, body)
    state[dofs] .= values
end