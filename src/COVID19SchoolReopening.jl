module COVID19SchoolReopening

using OrdinaryDiffEq, RecipesBase, DiffEqCallbacks

#################
#   constants   #
#################

const DEFAULT_DAYS = 5
const DEFAULT_PERIOD = 1
const DEFAULT_SHIFT = 0

const DEFAULT_SCHOOL_START = 9
const DEFAULT_SCHOOL_HOURS = 6
const DEFAULT_COHORT_BREAK = 2

const DEFAULT_K_VALUE = 5

const DAYS_PER_MONTH = 30.4167

#########################
#   helper functions    #
#########################

tanh_pulse(s,s₁,s₂,K) = (tanh(K*(s-s₁)) - tanh(K*(s-s₂))) / 2

"""
### Arguments

- `j`: affected age group
- `k`: affected cohort
- `α`: cohort mixing matrix
- `τ`: transmission risk matrix
- `c`: average contact rate matrix
- `u`: state matrix
"""
function force_of_infection(j, k, t, α, β, u)
    _, nclasses, ncohorts = size(u)
    λ = zero(promote_type(eltype(α), eltype(u)))

    for l in 1:ncohorts, i in 1:nclasses
        I_il = u[3,i,l]
        λ += α[k,l] * β[i,j](t) * I_il
    end

    return λ
end

function check_array_distribution(arr, arr_name)
    if any(<(0), arr) || any(>(100), arr)
        throw(error("Need `0 ≤ $(arr_name)[j] ≤ 100`"))
    elseif sum(arr) != 100
        throw(error("Need `sum($(arr_name)) == 100`"))
    end
end

"""
```
define_initial_conditions(total_infected, demography, cohorts, infected)
```

Generate an array `u0` where `u0[:,j,k]` represents compartments S_jk, E_jk, I_jk, R_jk. Note that `sum(u0) ≈ 1`.

### Arguments

- `total_infected:` Percent of individuals infected in overall population (0-100).
- `demography`: Percent of individuals in each age group (0-100).
- `cohorts`: Percent individuals in each cohort (0-100).
- `infected`: Percent infected in each age group (0-100).
"""
function define_initial_conditions(total_infected, demography, cohorts, infected)
    # sanity checks
    if total_infected < 0 && total_infected > 100
        throw(error("Need `0 ≤ total_infected ≤ 100`"))
    end

    for (arr, arr_name) in zip((demography, cohorts, infected), ("demography", "cohorts", "infected"))
        check_array_distribution(arr, arr_name)
    end

    # extract number of age classes and cohorts
    nclasses = length(demography)
    ncohorts = length(cohorts)

    # initialize output
    u = zeros(4, nclasses, ncohorts)

    # scale to [0,1] so we can convert quantities
    t = total_infected / 100
    p = demography / 100
    q = cohorts / 100
    r = infected / 100

    for k in 1:ncohorts, j in 1:nclasses
        N = p[j] * q[k]     # proportion in age class j, cohort k
        I = t * r[j] * q[k] # proportion infected in class j, cohort k
        S = N - I           # remaining mass belongs to susceptible compartment

        u[1,j,k] = S
        u[3,j,k] = I
    end

    return u
end

#####################
#   contact rates   #
#####################

struct TransmissionRate{T} <: Function
    cmax::T         # maximum
    cmin::T         # minimum
    days::Int       # number of days in a school week
    period::Int     # how often the schedule repeats (1 => every 7-day week)
    shift::Int      # period for repeating a school week (0 => every week)
    start::Float64  # starting hour (24 hour format)
    stop::Float64   # ending hour
    K::Float64      # smoothness parameter
    school::Bool    # indicates whether peaks correspond to rate at school
end

"""
```
TransmissionRate(cmax, cmin, start, stop; days::Int = 5, period::Int = 1, K::Int = 1)
```

Construct a function modeling a transmission rate between individuals within a cohort.

### Arguments

- `cmax`: maximum transmission rate
- `cmin`: minimum transmission rate
- `start`: Start of school, assuming units of days.
- `stop`: End of school, assuming units of days.

### Keyword Arguments

- `days`: Number of days in a school week, assuming a 7-day week. Defaults to `days = 5.`
- `period`: how often the school week repeats. Defaults to `period = 1`, meaning the school schedule repeats every 7-day week.
- `shift`: offset that shifts the schedule when working with multiple cohorts. Defaults to `shift = 0`, meaning the period aligns with simulation time.
- `K`: An integer modulating the smoothness of a square-wave approximation. Larger values correspond to better approximations. Defaults to `K = 1`.
- `school`: Indicates whether peaks represent rates at school.
"""
function TransmissionRate(cmax, cmin, start, stop; days::Int=5, period::Int=1, shift::Int=0, K::Real=1.0, school::Bool=true)
    return TransmissionRate(cmax, cmin, days, period, shift, start, stop, Float64(K), school)
end

import Base: show

function Base.show(io::IO, f::TransmissionRate)
    println(io, "Minimum value: $(f.cmin)")
    println(io, "Maximum value: $(f.cmax)")
    println(io, "Number of school days: $(f.days)")
    println(io, "------------------------------")
    println(io, "Period: $(f.period)")
    println(io, "Shift:  $(f.shift)")
    println(io, "------------------------------")
    println(io, "Period start: $(f.start)")
    println(io, "Period stop:  $(f.stop)")
    println(io, "Sharpness:    $(f.K)")
    println(io, "")
end

##### implementation #####

"""
Evaluate transmission rate at time `t`.
"""
function (f::TransmissionRate)(t)
    K = f.K

    # indicator = (v + f.shift) % f.period == 0
    v = t % (7*f.period)    # current day in schedule
    v₁ = 7*f.shift          # start of schedule
    v₂ = v₁ + f.stop        # end of school week

    # smooth approximation to the indicator above
    indicator = tanh_pulse(v,v₁,v₂,K)

    s = t % 7           # current day in a week
    s₁ = f.start-0.5    # start of school week
    s₂ = f.stop         # end of school week

    # smooth approximation to transmission function
    rate = (f.cmax - f.cmin) * tanh_pulse(s,s₁,s₂,K) * indicator + f.cmin
    if !(f.school)
        rate = f.cmax-rate
    end

    return rate
end

#
# visualization: used to check correctness and smoothness
#
@recipe function plot(f::TransmissionRate; weeks = 1, step = 0.25)
    w = weeks*7
    xs = range(0.0, stop = w, step = step)
    ys = f.(xs)

    day_ticks  = 0:1:w
    day_labels = 0:1:length(day_ticks)

    # formatting
    legend := false
    grid  --> false
    xguide --> "time (days)"
    xticks --> (day_ticks, day_labels)

    @series begin
        seriestype := :line
        xs, ys
    end
end

#################
#   ODE model   #
#################

#
#   implements the ODE model
#
function ode_model(du, u, p, t)
    (α, σ, γ, β, ncohorts, nclasses) = p

    for k in 1:ncohorts, j in 1:nclasses
        # evaluate force of infection and save its value
        λ = force_of_infection(j, k, t, α, β[k], u)

        # pull out compartment values for readability
        S = u[1,j,k]
        E = u[2,j,k]
        I = u[3,j,k]
        R = u[4,j,k]

        # ODE subsysem for class j and cohort k
        du[1,j,k] = dS = -λ*S
        du[2,j,k] = dE = λ*S - σ[j]*E
        du[3,j,k] = dI = σ[j]*E - γ[j]*I
        du[4,j,k] = dR = γ[j]*I
    end

    return nothing
end

"""
Generates a model instance from the given inputs.
By default, assumes each age class has equal representation in the population and each cohort is seeded with the same fraction of infecteds.

### Notes:

Model parameters can be accessed with `prob.p`.
The object `p` will be a tuple: (α, σ, γ, β, ncohorts, nclasses).
Initial conditions are given by `prob.u0`.

### Arguments

- `nclasses`: number of age classes.
- `ncohorts`: number of cohorts.
- `tspan`: timespan for the simulation.
- `beta`: baseline transmission rates.
- `params`: model parameters other than beta.

### Keyword Arguments

- `mult`: represents multiplier used in modeling increased school contacts in `β[1,1]`.
- `demographics`: used to specify different age structure.
- `cohorts`: used to specify distribution of populations across cohorts.
- `total_infected`: used to specify total infecteds in population.
- `infected`: used to specify distribution of infecteds across age classes.
- `tstart`: used to control beginning of peak in transmission function.
- `tstop`: used to control end of peak in tranmission function.
- `K`: sharpness parameter in approximating indicator functions.
"""
function make_model_instance(nclasses, ncohorts, tspan, beta, params;
    mult = 1.0,
    demographics = nothing,
    cohorts = nothing,
    total_infected = 1.0,
    infected = nothing,
    tstart=0.75,     # start for school week, days
    tstop=4.25,     # end for school week, days
    K=3.0,)         # smoothness parameter
    #
    # transmission rates
    #
    matT = Matrix{TransmissionRate{Float64}}
    β = [matT(undef, nclasses, nclasses) for k in 1:ncohorts]

    for k in 1:ncohorts
        # child → child
        β[k][1,1] = TransmissionRate(mult * beta[1,1], beta[1,1], tstart, tstop,
            period=ncohorts, shift=k-1, K=K)
        # child → adult
        β[k][1,2] = TransmissionRate(beta[1,2], beta[1,2], 0.0, 0.0)
        # adult → child
        β[k][2,1] = TransmissionRate(beta[2,1], beta[2,1], 0.0, 0.0)
        # adult → adult
        β[k][2,2] = TransmissionRate(beta[2,2], beta[2,2], 0.0, 0.0)
    end

    # assemble parameters vector
    p = (params..., β, ncohorts, nclasses)

    # intialize demographic structure / age classes
    if demographics === nothing
        q = 100 / nclasses .* ones(nclasses)
    else
        if length(demographics) != nclasses
            error("argument `nclasses` must match length of `demographics`")
        end

        q = 100 .* demographics ./ sum(demographics)
    end

    # initialize cohort sizes
    if cohorts === nothing
        ρ = 100 .* ones(ncohorts) ./ ncohorts
    else
        if length(cohorts) != ncohorts
            error("argument `ncohorts` must match length of `cohorts`")
        end
        ρ = 100 .* cohorts ./ sum(cohorts)
    end

    # initialize infection distribution across age classes
    if infected === nothing
        r = 100 .* ones(nclasses) ./ nclasses
    else
        if length(infected) != nclasses
            error("argument `nclasses` must match length of `infected`")
        end
        r = 100 .* infected ./ sum(infected)
    end

    # define initial conditions, and problem object
    u0 = define_initial_conditions(total_infected, q, ρ, r)
    prob = ODEProblem(ode_model, u0, tspan, p)

    return prob
end

##### Testing Policy #####
struct TestPolicy <: Function
    sensitivity::Float64
    isrotating::Bool
    detected::Vector{Float64}
end

# define affect
function (policy::TestPolicy)(integrator)
    if integrator.opts.userdata[:school_closed]
        return nothing
    end

    # retrieve transmission functions for children
    beta = integrator.p[4]

    # identify cohort currently in school
    cohort = 0

    for k in eachindex(beta)
        period = beta[k][1,1].period
        shift = beta[k][1,1].shift
        v = integrator.t % (7*period)
        if 7*shift+1 ≤ v+1 ≤ 7*(shift+1)
            cohort = k
        end
    end

    # remove identified cases from the infected pool
    if cohort > 0
        # exposed  = policy.sensitivity * integrator.u[2,1,cohort]
        infected = policy.sensitivity * integrator.u[3,1,cohort]
        # total = infected + exposed
        total = infected

        # integrator.u[2,1,cohort] -= exposed
        integrator.u[3,1,cohort] -= infected
        integrator.u[4,1,cohort] += total

        push!(policy.detected, total)
    end

    return nothing
end

function make_testing_policy_callback(test_times, detected;
    sensitivity::Real = 1.0,
    isrotating::Bool  = true,)
    # define the testing policy needed for checking threshold
    test_policy = TestPolicy(sensitivity, isrotating, detected)
    cb_test = PresetTimeCallback(test_times, test_policy)

    return cb_test
end

##### Prevalance Thresholds #####
struct ThresholdPolicy <: Function
    sensitivity::Float64
    isrotating::Bool
    detected::Vector{Float64}
    threshold::Float64
    window::Int
end

# define condition based on 3-arg signature
function (policy::ThresholdPolicy)(u, t, integrator)
    I, J, K = size(u)

    # account for size of children subpop.
    q1 = zero(eltype(u))
    q2 = zero(eltype(u))

    if policy.isrotating
        K = size(u, 3)  # check all cohorts
    else
        K = 1           # check only the in-person cohort
    end

    for i in 1:I, k in 1:K
        q1 += u[i,1,k]
        q2 += u[i,2,k]
    end

    # extract infection data within window
    n = length(policy.detected)
    offset = policy.window - 1
    x = @views sum(policy.detected[max(1, n-offset):n])

    return (x / q1 - policy.threshold) > 0
end

# define affect based on 1-arg signature
function (policy::ThresholdPolicy)(integrator)
    (α, σ, γ, β, ncohorts, nclasses) = integrator.p

    if policy.isrotating
        # iterate over all cohorts
        for k in eachindex(β)
            β11 = β[k][1,1]
            β[k][1,1] = TransmissionRate(β11.cmin, β11.cmin, β11.start, β11.stop)
        end
    else
        # change only for in-person cohort
        β11 = β[1][1,1]
        β[1][1,1] = TransmissionRate(β11.cmin, β11.cmin, β11.start, β11.stop)

        integrator.p = (α, σ, γ, β, ncohorts, nclasses)
    end

    # make sure change is reflected in simulation
    integrator.p = (α, σ, γ, β, ncohorts, nclasses)

    # make sure we don't trigger additional testing callbacks
    integrator.opts.userdata[:school_closed] = true

    return nothing
end

# define function to construct correct callback set
function make_threshold_policy_callbacks(test_times, detected;
    sensitivity::Real = 1.0,
    isrotating::Bool  = true,
    threshold::Real   = 0.05,
    window::Integer   = 14,
)
    # define the testing policy needed for checking threshold
    test_policy = TestPolicy(sensitivity, isrotating, detected)
    cb_test = PresetTimeCallback(test_times, test_policy)

    # define the threshold policy
    threshold_policy = ThresholdPolicy(sensitivity, isrotating, detected, threshold, window)
    cb_threshold = DiscreteCallback(threshold_policy, threshold_policy)

    # combine into CallbackSet - testing affect should trigger first!
    cbs = CallbackSet(cb_test, cb_threshold)

    # return all three objects
    return cb_test, cb_threshold, cbs
end

export define_initial_conditions, force_of_infection, tanh_pulse
export TransmissionRate
export make_model_instance, vary_within_group_transmission, vary_between_group_transmission
export TestingPolicy, make_testing_policy_callback
export ThresholdPolicy, make_threshold_policy_callbacks

end
