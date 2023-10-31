# Derivatives can be computed using manual implementation of automatic differentiation (ad)
# or finite differences

# First we define, then we set defaults for different types of BlochSimulators. For the EPG model, it turns out that finite difference is usually faster. For the Isochromat model, automatic differentiation is faster

# """
#     _simulate_derivatives_finite_difference(resource, sequence, parameters::TissueParameters{T}, Δ = T(10^-4)) where T

# Currently hard-coded to simulate T₁ and T₂ derivatives using finite differences. A B₁ map is assumed to be present
# """
# function _simulate_derivatives_finite_difference(resource, sequence, parameters::AbstractArray{<:AbstractTissueParameters{T}}, Δ = T(10^-4)) where T

#     m = simulate_echos(resource, sequence, parameters)

#     # derivatives w.r.t. T₁
#     Δpars = map(parameters) do p
#         B₁ = hasB₁(p) ? p.B₁ : one(T)
#         BlochSimulators.T₁T₂B₁(p.T₁+Δ, p.T₂, B₁)
#     end

#     Δm = simulate_echos(resource, sequence, Δpars)
#     ∂m∂T₁ = @. (Δm - m)/Δ

#     # derivatives w.r.t. T₂
#     Δpars = map(parameters) do p
#         B₁ = hasB₁(p) ? p.B₁ : one(T)
#         BlochSimulators.T₁T₂B₁(p.T₁, p.T₂+Δ, B₁)
#     end

#     Δm = simulate_echos(resource, sequence, Δpars)
#     ∂m∂T₂ = @. (Δm - m)/Δ

#     # turn separate arrays into appropriate array-of-structs
#     ∂echos_d = map((m,∂T₁,∂T₂)-> (m,∂mˣʸ∂T₁T₂(∂T₁,∂T₂)), m, ∂m∂T₁, ∂m∂T₂)

#     return ∂echos_d
# end

function _simulate_derivatives_ad_manual(::CPU1, sequence::BlochSimulator{T}, parameters::AbstractArray{<:AbstractTissueParameters{T}}) where T

    nr_voxels = length(parameters)
    nr_echos = length(sequence.RF_train)
    # intialize array to store magnetization + derivative w.r.t. all nonlinear parameters
    # at echo times for each voxel
    ∂echo_type = Tuple{complex(T), ∂mˣʸ_type(parameters[1]){complex(T)}}
    ∂echos = zeros(∂echo_type, nr_echos, nr_voxels)
    # initialize array of magnetization vectors + vectors to store partial derivatives
    # one set for each spin in the slice direction
    states = initialize_states_derivatives(sequence, parameters[1])
    # loop over voxels
    for i = 1:nr_voxels
        # all functions used in the sequence have methods for states that include partial derivatives
        BlochSimulators.simulate_echos!(view(∂echos,:,i), sequence, states, parameters[i])
    end

    return ∂echos
end

# function _simulate_derivatives_finite_difference(sequence::BlochSimulator{T}, parameters::AbstractArray{<:AbstractTissueParameters{T}}; Δ = T(10^-4)) where {T}

#     nr_voxels = length(parameters)
#     nr_echos = length(sequence.RF_train)
#     # intialize array to store echos for each voxel
#     echos  = zeros(Complex{T}, nr_echos, nr_voxels)
#     # intialize array to store ∂echos for each voxel
#     ∂echos = zeros(Complex{T}, nr_echos, nr_nonlinpars(parameters), nr_voxels)
#     # to store magnetization at echo times with slightly perturbed parameters
#     buffer = zeros(Complex{T}, nr_echos)

#     states = initialize_states(sequence)

#     # loop over voxels
#     for i = 1:nr_voxels
#         # start simulation for voxel and store echos
#         mˣʸ = view(echos,:,i)
#         ∂mˣʸ = view(∂echos,:,:,i)
#         # simulate magnetization at current parameters
#         BlochSimulators.simulate_echos!(mˣʸ, sequence, states, parameters[i])
#         # compute derivatives using finite difference method
#         _finite_difference!(∂mˣʸ, mˣʸ, sequence, buffer, parameters[i], states, Δ)

#     end

#     return [ (echos[r,v], ∂mˣʸ_type(parameters[1])(∂echos[r,:,v]...)) for r in 1:nr_echos(sequence), v in 1:nr_voxels]
# end

# function _finite_difference!(∂mˣʸ, mˣʸ, sequence::BlochSimulator{T}, buffer, p::AbstractTissueParameters{T}, states, Δ) where {T}
#     # does some allocations but not so bad for performance I guess

#     # for all nonlinear parameters do finite difference
#     idx = 1
#     parnames = fieldnames(typeof(p))

#     for par in setdiff(parnames, [:ρ]) # exclude proton density

#         T₁ = p.T₁ + (par == :T₁ ? Δ : 0)
#         T₂ = p.T₂ + (par == :T₂ ? Δ : 0)
#         q = BlochSimulators.T₁T₂(T₁, T₂)

#         if hasB₁(p)
#             B₁ = p.B₁ + (par == :B₁ ? Δ : 0)
#             q = BlochSimulators.add_B₁(q, B₁)
#         end

#         if hasB₀(p)
#             B₀ = p.B₀ + (par == :B₀ ? Δ : 0)
#             q = BlochSimulators.add_B₁(q, B₀)
#         end

#         if :ρ ∈ parnames
#             q = BlochSimulators.add_ρ(q, ρ)
#         end

#         # Simulate with new parameters and compute finite difference
#         BlochSimulators.simulate_echos!(buffer, sequence, states, q)
#         @. ∂mˣʸ[:,idx] = (buffer - mˣʸ) ./ Δ
#         buffer[:] .= zero(Complex{T})
#         idx += 1
#     end
# end


# function _simulate_derivatives_ad_forwarddiff(simulator, parameters::AbstractArray)

#     # parameters should be an array of NamedTuples
#     # needs to be converted to an array of SVectors because ForwardDiff
#     # does not work with NamedTuples

#     # make array of SVectors instead of array of NamedTuples (ForwardDiff doesnt work with NamedTuples)
#     # npars = length(parameters[1])
#     parameters = [SVector{4}(x[1], x[2], x[3], x[4]) for x in parameters]
#     # intialize array to store echos for each voxel
#     nvoxels = length(parameters)
#     nTRs = length(simulator.RF_train)
#     echotype = complex(eltype(parameters[1]))

#     # allocate output for magnetization at echo time as well as
#     # partial derivatives of magnetization at echo time w.r.t. (nonlinear) parameters
#     # echos = zeros(real(echotype), 2*nTRs, nvoxels)
#     ∂echos = zeros(real(echotype), 2*nTRs, NONLINPARS + 1, nvoxels)

#     # jac = (result, p) -> result = ForwardDiff.jacobian!(result,x -> _wrapper(simulator,x), p)

#     # loop over voxels
#     for i = 1:nvoxels
#         # create DiffResult to store echos and derivatives at the same time (syntax is a bit weird but ok)
#         result = DiffResults.DiffResult(view(∂echos,:,1,i), view(∂echos,:,2:NONLINPARS+1,i))
#         # start simulation for voxel and store echos
#         result = ForwardDiff.jacobian!(result,x -> _wrapper(simulator,x), parameters[i]);
#         # jac(result,parameters[i])
#     end

#     # make complex arrays again
#     # echos = reinterpret(echotype,echos)
#     ∂echos = reinterpret(echotype,∂echos)

#     ∂echos = permutedims(∂echos,[2 1 3])

#     return ∂echos
# end

# function _wrapper(simulator, pars)
#     # ForwardDiff works with SVectors but BlochSimulators with NamedTuples:
#     p = (T₁=pars[1], T₂=pars[2], B₁=pars[3], B₀=pars[4])
#     m = BlochSimulators.simulate_echos(simulator, p)
#     # m is complex array, but ForwardDiff only works with real arrays:
#     m = reinterpret(real(eltype(m)),m)
#     return m
# end

# # function simulate_derivatives_ad_forwarddiff(
# #     simulator::BlochSimulators.BlochSimulator{T},
# #     parameters::AbstractVector{NamedTuple{(:T₁, :T₂, :B₁, :B₀, :ρˣ, :ρʸ),NTuple{6,T}}}
# #     ) where {T}

# #     parameters = [(T₁=p[1], T₂=p[2], B₁=p[3], B₀=p[4]) for p in parameters]
# #     simulate_derivatives_ad_forwarddiff(simulator, parameters)
# # end

# Set defaults for different types of BlochSimulators
# default setting though
# function simulate_derivatives(resource, sequence::IsochromatSimulator, parameters)
#     ∂echos = _simulate_derivatives_ad_manual(resource, sequence, parameters)
#     return ∂echos
# end

# function simulate_derivatives(resource, sequence::EPGSimulator, parameters)
#     ∂echos = _simulate_derivatives_finite_difference(resource, sequence, parameters)
#     return ∂echos
# end

function simulate_derivatives(m, resource, sequence, simulation_parameters, fit_parameters)
    _simulate_derivatives_finite_difference(m, resource, sequence, simulation_parameters, fit_parameters)
end

# Distributed version, choice made above propagates through
# Assumes parameters have already been distributed
function simulate_derivatives(::CPUProcesses, sequence, dparameters::DArray)
    # start computing ∂echos on each worker
    futures = [@spawnat p simulate_derivatives(CPU1(), sequence, dparameters[:lp]) for p in workers()]
    # sync, on each core a 2D array is returned so we need to reshape the futures array accordingly
    d∂echos = DArray(reshape(futures, (1,length(futures))))
    return d∂echos
end

# Distributed version for normal array
simulate_derivatives(resource::CPUProcesses, sequence, parameters) = simulate_derivatives(resource, sequence, distribute(parameters))


    # AD on GPU stuff

    # function cuda_kernel_derivatives!(∂echos, sequence, parameters)

    #     voxel = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    #     if voxel <= length(parameters)

    #         states = initialize_states_derivatives(sequence, parameters[voxel])
    #         simulate!(view(∂echos,:,voxel), sequence, states, parameters[voxel])

    #     end
    #     nothing
    # end

    # function simulate_derivatives(::CUDALibs, sequence::BlochSimulator{T}, parameters::AbstractArray{<:AbstractTissueParameters{T}}; threads_per_block=64) where T

    #     nr_voxels = length(parameters)
    #     nr_echos = length(sequence.RF_train) * echos_per_TR(sequence)

    #     # intialize array to store magnetization + derivative w.r.t. all nonlinear parameters
    #     # at echo times for each voxel
    #     ∂echo_type = Tuple{complex(T), ∂mˣʸ_type(parameters[1]){complex(T)}}
    #         ∂echos_d = fill!(CuArray{∂echo_type}(undef, nr_echos, nr_voxels), zero(∂echo_type))

    #     # @time ∂echos_d = CuArrays.zeros(∂echo_type, nr_echos, nr_voxels)
    #     parameters_d = CuArray(parameters)
    #     sequence_d = BlochSimulators.struct_to_gpu(sequence)

    #     nr_blocks = cld(nr_voxels, threads_per_block)

    #     # launch kernels
    #     CUDA.@sync begin
    #         @cuda blocks=nr_blocks threads=threads_per_block cuda_kernel_derivatives!(∂echos_d, sequence_d, parameters_d)
    #     end

    #     return ∂echos_d
    # end

# recon_options["simulation_parameters"]
# recon_options["fit_parameters"] = (:T₁, :T₂, :B₁⁺, :ΔB₀)

# recon_options["fit_parameters"] = (:T₁, :T₂, :B₁⁺, :ΔB₀)


# function _simulate_derivatives_finite_difference(m, resource, sequence, parameters::AbstractArray{<:T₁T₂B₁{T}}, fitpars::T₁T₂{T}, Δ = T(10^-4)) where T


#     # derivatives w.r.t. T₁
function _simulate_derivatives_finite_difference(m, resource, sequence, simulation_parameters::AbstractArray{<:AbstractTissueParameters{N,T}}, nonlinear_derivatives::Type{<:∂mˣʸ}, Δ = T(10^-4)) where {N,T}

    # Accessors.@set on GPU gives some issues so convert to CPU array first
    simulation_parameters = collect(simulation_parameters)

    ∂echos = map(fieldnames(nonlinear_derivatives)) do derivative

        Δsimulation_parameters = map(simulation_parameters) do p
            if derivative == :∂T₁
                return @set p.T₁ += Δ
            elseif derivative == :∂T₂
                return @set p.T₂ += Δ
            elseif derivative == :∂B₁
                return @set p.B₁ += Δ
            elseif derivative == :∂B₀
                return @set p.B₀ += Δ
            else
                return p
                @warn "Parameter not yet implemented"
            end
        end

        # send to GPU
        if resource == CUDALibs()
            Δsimulation_parameters = cu(Δsimulation_parameters)
        end

        # run simulations
        Δm = simulate_echos(resource, sequence, Δsimulation_parameters)

        # compute finite differences
        @. (Δm - m)/Δ
    end

    ∂m = map(nonlinear_derivatives{Complex{T}}, ∂echos...)

    return map(tuple, m, ∂m)
end


function simulate_derivatives_structarray(m, resource, sequence, simulation_parameters::AbstractArray{<:AbstractTissueParameters{N,T}}, nonlinear_derivatives::Type{<:∂mˣʸ}, Δ = T(10^-4)) where {N,T}

    # Accessors.@set on GPU gives some issues so convert to CPU array first
    simulation_parameters = collect(simulation_parameters)
    ∂m = map(fieldnames(nonlinear_derivatives)) do derivative

        Δsimulation_parameters = map(simulation_parameters) do p
            if derivative == :∂T₁
                return @set p.T₁ += Δ
            elseif derivative == :∂T₂
                return @set p.T₂ += Δ
            elseif derivative == :∂B₁
                return @set p.B₁ += Δ
            elseif derivative == :∂B₀
                return @set p.B₀ += Δ
            else
                return p
                @warn "Parameter not yet implemented"
            end
        end

        # send to GPU
        if resource == CUDALibs()
            Δsimulation_parameters = cu(Δsimulation_parameters)
        end

        # run simulations
        Δm = simulate_echos(resource, sequence, Δsimulation_parameters)

        # compute finite differences
        @. (Δm - m)/Δ
    end
    ∂m = StructArray{nonlinear_derivatives}(∂m)
    return ∂m
end

# struct m∂mˣʸ∂T₁∂T₂{T} <: ∂mˣʸ{3,T}
    # m::T
    # ∂T₁::T
    # ∂T₂::T
# end
