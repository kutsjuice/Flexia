function setid!(conn::AbstractConnector2D, index)
    conn.index = index
end

function add!(sys::MBSystem2D, connector::AbstractConnector2D)
    push!(sys.connectors, connector)
    # sys.connectorsnum += 1;
    last_connector_dof = last_lm_dof(sys) + numberofdofs(connector)
    push!(sys.lmdofs, last_connector_dof)
    setid!(connector, length(sys.connectors))
end