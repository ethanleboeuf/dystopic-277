
using StaticArrays, LightGraphs, MetaGraphs, LazySets, Random, POMDPs, Combinatorics
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
const AGENT_ACTIONS = [Vec2(i, j) for i in -1:1 for j in -1:1 if !(i==j==0)]
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

        inf_vis = inf_visibile(G, obstacles)
        visibility = limited_visible(G, inf_vis)

        return new(size, obstacles, visibility, G)
    end
end

Base.size(env::Environment) = env.size
graph(env::Environment) = env.G
obstacles(env::Environment) = env.obstacles


Base.rand(G::AbstractGraph, rng::AbstractRNG = Random.GLOBAL_RNG) = rand(rng, 1:nv(G))
function rand_pos(G::AbstractGraph, rng)
    n = rand(G, rng)
    get_pos_tup(G, n)
end

struct Region
    center::Tuple # maybe this one isn't necessary
    interior::Vector{Tuple}
end


####################################################
#################### S, A, O #######################
####################################################
mutable struct Agent
    pos::Vec2
    partition::Region
end

partition(ag::Agent) = ag.partition
position(ag::Agent)  = ag.pos
pos(ag::Agent)       = ag.pos

struct SNState
    patrollers::Vector{Agent}
    adversaries::Vector{Vec2}
end
SNAction = Union{Vector{Vec2}, Symbol}

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

struct Random <: AdversaryModel end
struct MaximizeDistance <: AdversaryModel end
struct KnowsBelief{APB<:AbstractParticleBelief} <: AdversaryModel
    belief::APB
end
# There was another one wasn't there?

function heuristic_update(::MaximizeDistance, adv::Agent)

end

####################################################
############# POMDP type + functions ###############
####################################################
mutable struct SkyNet <: POMDP{SNState, SNAction, SNObs}
    env::Environment
    partition_points::Vector{Vector{Tuple}} # maybe don't need
    partitions::Vector{Region}              # maybe don't need
    max_adversaries::Int8                   # maybe don't need
    intrusion_prob::Float64                 # maybe don't need
    R::Int64
    adversarial_model::AdversaryModel

    SkyNet() = new()
end

environment(p::SkyNet) = p.env
graph(p::SkyNet) = p.env.G
obstacles(p::SkyNet) = p.env.obstacles
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
        error("State is terminal in generate_s. Why is this happening? $s")
        # return s
    end
    if a == :alarm
        out_of_bounds_pats = patrollers(s) .+ Vec2(size(environment(p)))
        return SNState(out_of_bounds_pats, adversaries(s))
    end

    patrollers, adversaries = agents(s)
    updated_patrollers = patrollers .+ a

    # update adversary's position somehow:
    new_adv = heuristic_update(adversarial_model(p), adversaries)

    return SNState(patrollers, adversaries)
end


function generate_o(p::SkyNet, s::SNState, rng::AbstractRNG)
    noise = rand(p.d, rng)


end


function reward(p::SkyNet, s::SNState, a::SNAction, sp::SNState)

end