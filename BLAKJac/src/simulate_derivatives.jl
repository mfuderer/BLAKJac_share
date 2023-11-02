# Derivatives can be computed using manual implementation of automatic differentiation (ad)
# or finite differences

# First we define, then we set defaults for different types of BlochSimulators. For the EPG model, it turns out that finite difference is usually faster. For the Isochromat model, automatic differentiation is faster

# """
#     _simulate_derivatives_finite_difference(resource, sequence, parameters::TissueParameters{T}, Δ = T(10^-4)) where T

# Currently hard-coded to simulate T₁ and T₂ derivatives using finite differences. A B₁ map is assumed to be present
# """



function simulate_derivatives(m, resource, sequence, simulation_parameters)
    _simulate_derivatives_finite_difference(m, resource, sequence, simulation_parameters)
end



function _simulate_derivatives_finite_difference(m, resource, sequence, parameters)

    Δ = 10^-4

    # derivatives w.r.t. T₁
    Δpars = map(parameters) do p
        B₁ = hasB₁(p) ? p.B₁ : 1.0
        BlochSimulators.T₁T₂B₁(p.T₁+Δ, p.T₂, B₁)
    end
    Δm = simulate_echos(resource, sequence, Δpars)
    ∂m∂T₁ = @. (Δm - m)/Δ

    # derivatives w.r.t. T₂
    Δpars = map(parameters) do p
        B₁ = hasB₁(p) ? p.B₁ : 1.0
        BlochSimulators.T₁T₂B₁(p.T₁, p.T₂+Δ, B₁)
    end
    Δm = simulate_echos(resource, sequence, Δpars)
    ∂m∂T₂ = @. (Δm - m)/Δ

    # derivatives w.r.t. B₁
    Δpars = map(parameters) do p
        B₁ = hasB₁(p) ? p.B₁ : 1.0
        BlochSimulators.T₁T₂B₁(p.T₁, p.T₂, B₁+Δ)
    end
    Δm = simulate_echos(resource, sequence, Δpars)
    ∂m∂B₁ = @. (Δm - m)/Δ

    return ∂m∂T₁, ∂m∂T₂, ∂m∂B₁
end

