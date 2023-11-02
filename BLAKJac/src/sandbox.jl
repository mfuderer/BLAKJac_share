struct kwereks{momo} <: Any
    juba::momo
    jubo::momo
end



struct FISP3D{T<:AbstractFloat, Ns, U<:AbstractVector{Complex{T}}} <: EPGSimulator{T,Ns}
    RF_train::U # flipangles in degrees, no phase
    TR::T # repetition time in seconds
    TE::T # echo time in seconds
    max_state::Val{Ns} # maximum number of states to keep track of in EPG simulation
    TI::T # inversion delay in seconds
    TW::T     # waiting time between repetitions in seconds
    repetitions::Int # number of repetitions
    inversion_prepulse::Bool    # with or without inversion
    wait_spoiling::Bool    # spoiling is assumed during waiting time 
end

