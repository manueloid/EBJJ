using EBJJ, QuantumOptics
function Hamiltonian(q::ConstantQuantity, c::Control)
    Ω(t) = control_function(t, c)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    return H
end
function ξN(q::ConstantQuantity, c::Control)
    H = Hamiltonian(q, c)
    ΔJ(t, psi) = (dagger(psi) * (q.Jz)^2 * psi - (dagger(psi) * (q.Jz) * psi)^2) / (c.N / 4) |> real
    return timeevolution.schroedinger_dynamic([0.0:c.T/1000:c.T;], q.ψ0, H; fout=ΔJ)
end
