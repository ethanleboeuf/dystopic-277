
using StaticArrays, LightGraphs, MetaGraphs,
LazySets, Random, POMDPs, Combinatorics, Parameters
include("ndgrid.jl")
#=
S :: vector of agents and pos of Adversary
A :: a rearrangement of the agents

O ::
    an adversarial position,
    noisy proportional to distance,
    every unobserved node get's a prob

R ::
    seeing Adversary (closer is better),
    minimize uncertainty,
    maximize overall overage,
    maximize agent utilization
=#

const Vec2 = SVector{2, Int64}
const AGENT_ACTIONS = [Vec2(i, j) for i in -1:1 for j in -1:1]
############################################
#########  Environment -related  ###########
############################################
struct Environment
    size::NTuple{2, Int}
    obstacles::Vector{Hyperrectangle}
    visibility::Vector{Vector{Tuple}}
    G::MetaGraph

    function Environment(size::NTuple{2, Int}, obstacles::Vector{Hyperrectangle}, R::Int64)
        G = create_graph(size...)
        remove_obstacles(G, obstacles)
        visibility = limited_visible(G, inf_visibile(G, obstacles))
        return new(size, obstacles, visibility, G)
    end
end

Base.size(env::Environment) = env.size
graph(env::Environment) = env.G
obstacles(env::Environment) = env.obstacles


Base.rand(G::AbstractGraph, rng::AbstractRNG = Random.GLOBAL_RNG) = rand(1:nv(G))
function rand_pos(G::AbstractGraph, rng::AbstractRNG = Random.GLOBAL_RNG)
    n = rand(rng, G)
    Vec2(get_pos_tup(G, n))
end
rand_pos(G::AbstractGraph) = rand_pos(Random.GLOBAL_RNG, G)

function rand_pos(G::AbstractGraph, n::Integer, rng::AbstractRNG)
    points = Set{Vec2}()
    while length(points) < n
        push!(points, rand_pos(G, rng))
    end
    return points
end


struct Region
    center::Tuple # maybe this one isn't necessary
    interior::Vector{Tuple}
end


####################################################
#################### S, A, O #######################
####################################################
struct Agent
    pos::Vec2
    partition::Region
end

partition(ag::Agent) = ag.partition # not sure how much this will be used.
pos(ag::Agent)       = ag.pos

struct SNState
    patrollers::Vector{Agent}
    adversaries::Vector{Vec2}
end

struct Alarm
    guess::Vec2
end
SNAction = Union{Vector{Vec2}, Alarm}

struct SNObs
    adversaries::Vector{Vec2}
end

agents(s::SNState)      = (s.patrollers, s.adversaries)
patrollers(s::SNState)  = s.patrollers
adversaries(s::SNState) = s.adversaries

####################################################
################ Adversary related #################
####################################################

abstract type AdversaryModel end

struct RandomActions <: AdversaryModel end
struct MaximizeDistance <: AdversaryModel end
struct MaximizeSafety <: AdversaryModel end
# struct KnowsBelief{APB<:AbstractParticleBelief} <: AdversaryModel
#     belief::APB
# end
# There was another one wasn't there?

function heuristic_update(::MaximizeDistance, adv::Agent)

end

####################################################
############# POMDP type + functions ###############
####################################################
@with_kw mutable struct SkyNet <: POMDP{SNState, SNAction, SNObs}
    env::Environment                    = Environment((20, 20), Hyperrectangle(low = [4.0, 4.0], high = [8.0, 8.0]), R)
    partition_points::Vector{Vec2}      = rand_pos(graph(env), n_agents) # maybe don't need
    partitions::Vector{Region}          = voronoi_cells(graph(env), obstacles(env), partition_points)
    max_adversaries::Int8               = 1      # maybe don't need
    intrusion_prob::Float64             = 1.0    # maybe don't need
    R::Int64                            = 4
    adversarial_model::AdversaryModel   = RandomActions
    n_agents::Int8                      = 4
    n_adversaries::Int8                 = 1      # maybe don't need
end

environment(p::SkyNet)       = p.env
graph(p::SkyNet)             = p.env.G
obstacles(p::SkyNet)         = p.env.obstacles
adversarial_model(p::SkyNet) = p.adversarial_model

"""
    isterminal checks if the nodes are outside the environment,
    which is used as a flag for being terminal. Should only happen
    if they were intentionally set that way at the end of generate_s,
    if they chose to sound the alarm.
"""
isterminal(p::SkyNet, s::SNState) = any(pat->pos(pat) > size(environment(p)), patrollers(s))

# TODO make it so that the alarm has to guess the right
# position, otherwise periodic alarms always win,
function actions(p::SkyNet, s::SNState)
    agentwise_action_space = [actions(p, pat) for pat in patrollers(s)]
    A = hcat.(vec.(ndgrid(agentwise_action_space...))...)
    push!(A, :alarm)
    return A
end

# maybe use neighbors instead of dir.
function actions(p::SkyNet, s::Agent)
    A = eltype(action_type(SkyNet))[]
    for dir in AGENT_ACTIONS
        new_pos = dir + pos(s)
        if new_pos ∈ partition(s) && all(new_pos .∉ obstacles(p))
            push!(A, dir)
        end
    end
    return A
end






function generate_s(p::SkyNet, s::SNState, a::SNAction, rng::AbstractRNG)

    if isterminal(s)
        error("State is terminal in generate_s. Why is this happening? \ns = $s, \n a = $a")
        # return s
    end
    if a isa Alarm
        out_of_bounds_pats = patrollers(s) .+ Vec2(size(environment(p)))
        return SNState(out_of_bounds_pats, adversaries(s))
    end

    patrollers, adversaries = agents(s)
    updated_patrollers = patrollers .+ a

    # update adversary's position somehow:
    new_adv = heuristic_update(adversarial_model(p), adversaries)

    return SNState(patrollers, adversaries)
end


# NOTE or come up with observation function
function generate_o(p::SkyNet, s::SNState, rng::AbstractRNG)
    noise = rand(p.d, rng)


end


function reward(p::SkyNet, s::SNState, a::SNAction, sp::SNState)
    # subtract some points for each adversary in the scene.
    r = -20.0 * length(adversaries(s))
    if a isa Alarm
        # Alarms are very costly.
        r -= 1000.0
        # If we guess the right location, we don't pay a price!
        # If we're close we pay half price
        dist = norm.([a.guess] .- adversaries(s))
        if any(dist .== 0)
            r += 1000.0
        elseif any(dist < sqrt(2)) # sqrt is diagonal distance
            r += 500.0
        end
    else
        travelled = norm.(a)  # gives, 0, 1, or sqrt(2) depeding on direction of travel
        r -= 5*sum(travelled)
    end
    return r
end