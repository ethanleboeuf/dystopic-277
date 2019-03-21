
using StaticArrays, LightGraphs, MetaGraphs,
LazySets, Random, POMDPs, Combinatorics,
Parameters, POMDPModelTools, Distributions, StaticArrays
include("ndgrid.jl")

import POMDPs: actions, states, discount, isterminal, obs_weight, generate_s, generate_o, reward
#=
S :: vector of agents and pos of Adversary
A :: a rearrangement of the agents or an Alarm

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
const SENSOR_MODEL = Dict(:o1s1 => BoolDistribution(0.8),
                          :o1s0 => BoolDistribution(0.2))

############################################
#########  Environment -related  ###########
############################################
struct Environment
    size::NTuple{2, Int}
    obstacles::Vector{<:Hyperrectangle}
    visibility::Vector{Vector{Tuple}}
    G::MetaGraph
    pos_to_ind::Dict{Vec2, Int64}

    function Environment(size::NTuple{2, Int}, obstacles::Vector{<:Hyperrectangle}, R::Int64)
        G = create_graph(size...)
        remove_obstacles(G, obstacles)
        visibility = limited_visible(G, obstacles, R, inf_visible(G, obstacles))
        node_mapping = Dict(get_pos(G, i) => i for i in 1:nv(G))
        return new(size, obstacles, visibility, G, node_mapping)
    end
end

Base.size(env::Environment) = env.size
graph(env::Environment) = env.G
obstacles(env::Environment) = env.obstacles
visibility(env::Environment) = env.visibility
visibility(env::Environment, pos::Vec2) = env.visibility[env.pos_to_ind[pos]]

get_pos(G::MetaGraph, i::Int64) = Vec2(get_pos_tup(G, i))


Base.rand(G::AbstractGraph, rng::AbstractRNG = Random.GLOBAL_RNG) = rand(1:nv(G))
function rand_pos(G::AbstractGraph, rng::AbstractRNG = Random.GLOBAL_RNG)
    i = rand(G, rng)
    Vec2(get_pos_tup(G, i))
end
rand_pos(G::AbstractGraph) = rand_pos(G, Random.GLOBAL_RNG)

function rand_pos(G::AbstractGraph, n::Integer, rng::AbstractRNG = Random.GLOBAL_RNG)
    points = Set{Vec2}()
    while length(points) < n
        push!(points, rand_pos(G, rng))
    end
    return collect(points)
end


struct Region
    # center::Tuple # maybe this one isn't necessary
    interior::Vector{Tuple}
end

in_obstacle(env::Environment, x) = in_obstacle(env::Environment, Vector{Float64}(x))
in_obstacle(env::Environment, x::Vector{Float64}) = any(x ∈ o for o in obstacles(env))

####################################################
#################### S, A, O #######################
####################################################
struct Agent
    pos::Vec2
    partition::Region
end

partition(ag::Agent) = ag.partition # not sure how much this will be used.
pos(ag::Agent)       = ag.pos
Agent(ag::Agent, pos::Vec2) = Agent(pos, ag.partition)
Base.:+(ag::Agent, a::Vec2) = Agent(ag, ag.pos+a)

struct SNState
    patrollers::Vector{Agent}
    adversaries::Vector{Vec2}
end

struct Alarm
    guess::Vec2
end

SNAction = Union{Vector{Vec2}, Alarm}
SNObservation = BitArray

agents(s::SNState)      = (s.patrollers, s.adversaries)
patrollers(s::SNState)  = s.patrollers
adversaries(s::SNState) = s.adversaries
# occupied_probability(o::SNObservation, pos::Vec2) = o.occupied[pos]

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


####################################################
############# POMDP type + functions ###############
####################################################
@with_kw mutable struct SkyNet <: POMDP{SNState, SNAction, SNObservation}
    size::NTuple{2, Int64}              = (20, 20)
    R::Int64                            = 4
    n_agents::Int8                      = 4
    discount_factor::Float64            = 0.9
    env::Environment                    = Environment(size, [Hyperrectangle(low = [4.0, 4.0], high = [8.0, 8.0])], R)
    partition_points::Vector{Vec2}      = rand_pos(graph(env), n_agents) # maybe don't need
    partitions::Vector{Region}          = Region.(voronoi_cells(graph(env), obstacles(env), partition_points))
    # max_adversaries::Int8               = 1      # maybe don't need
    # n_adversaries::Int8                 = 1      # maybe don't need
    # intrusion_prob::Float64             = 1.0    # maybe don't need
    adversarial_model                   = RandomActions()
    # sensor_distr::Truncated{Normal{Float64},Continuous} = Truncated(Normal(0, -R, R))
end

environment(p::SkyNet)       = p.env
graph(p::SkyNet)             = p.env.G
obstacles(p::SkyNet)         = p.env.obstacles
adversarial_model(p::SkyNet) = p.adversarial_model
discount(p::SkyNet)          = p.discount_factor

"""
    isterminal checks if the nodes are outside the environment,
    which is used as a flag for being terminal. Should only happen
    if they were intentionally set that way at the end of generate_s,
    if they chose to sound the alarm.
"""
isterminal(p::SkyNet, s::SNState) = any(pat->pos(pat) > Vec2(size(environment(p))), patrollers(s))



####################################################
############### More Adversary stuff ###############
####################################################

heuristic_update(p, s, rng) = heuristic_update(adversarial_model(p), p, s, rng)
function heuristic_update(::MaximizeDistance, p::SkyNet, s::SNState, rng)
end
function heuristic_update(::RandomActions, p::SkyNet, s::SNState, rng)
    G = graph(p)
    newadv = similar(adversaries(s))
    for (i, adv) in enumerate(adversaries(s))
        adv2 = adv + rand(rng, AGENT_ACTIONS)
        if all(adv2 .> 0) && all(adv2 .<= Vec2(size(G))) && !in_obstacle(environment(p), adv2)
            newadv[i] = adv2
        else
            newadv[i] = adv
        end
    end
    return newadv
end








function actions(p::SkyNet)
    println("wrong one")
    agentwise_action_space = [AGENT_ACTIONS for i in 1:p.n_agents]
    A = tuple.(vec.(ndgrid(agentwise_action_space...))...)
    return SNAction[collect.(A); alarm_actions(p)]
end

# TODO make it so that the alarm has to guess the right
# position, otherwise periodic alarms always win,
function actions(p::SkyNet, s::SNState)
    agentwise_action_space = [actions(p, pat) for pat in patrollers(s)]
    A = tuple.(vec.(ndgrid(agentwise_action_space...))...)
    return SNAction[collect.(A); alarm_actions(p)]
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

function alarm_actions(p::SkyNet)
    G = graph(p)
    [Alarm(Vec2(get_pos_tup(G, i))) for i in vertices(G)]
end




function generate_s(p::SkyNet, s::SNState, a::SNAction, rng::AbstractRNG)

    if isterminal(p, s)
        error("State is terminal in generate_s. Why is this happening? \ns = $s, \n a = $a")
        # return s
    end
    if a isa Alarm
        out_of_bounds_pats = patrollers(s) .+ [Vec2(size(environment(p)))]
        return SNState(out_of_bounds_pats, adversaries(s))
    end

    pats, advs = agents(s)
    updated_pats = pats .+ a
    # update adversary's position somehow:
    new_adv = heuristic_update(p, s, rng)

    return SNState(updated_pats, new_adv)
end


function POMCPOW.obs_weight(p::SkyNet, s::SNState, a::SNAction, sp::SNState, o::SNObservation)
    advs = adversaries(s)
    prob = 1.0
    for oᵢ in o
        if oᵢ ∈ advs
            pt = SENSOR_MODEL[:o1s1].p
        else
            pt = 1 - SENSOR_MODEL[:o1s0].p
        end

        pt = oᵢ ? pt : 1-pt

        prob *= pt
    end
    return prob
end

function generate_o(p::SkyNet, s::SNState, a::SNAction, sp::SNState, rng::AbstractRNG)

    patrollers, adversaries = agents(sp)

    println(pos.(patrollers))
    if isterminal(p, first(patrollers))
        println("Terminal!!")
        return BitArray(true)
    end

    o = BitArray(undef, 0)

    for pat in patrollers
        visibles = visibility(environment(p), pos(pat))
        for v in visibles
            isin = Vec2(v) ∈ adversaries ? :o1s1 : :o1s0
            push!(o, rand(rng, SENSOR_MODEL[isin]))
        end
    end
    return o
end

# function observation(p::SkyNet, sp::SNState, rng::AbstractRNG)
#     patrollers, adversaries = agents(sp)

#     pat = first(patrollers)
#     adv = first(adversaries)

#     for node in sensor_region(pat)
#         prob = rand(rng, sensor_distr(p))
#         if adv == node
#         else
#         end
#     end
# end


function reward(p::SkyNet, s::SNState, a::SNAction, sp::SNState)
    # subtract some points for each adversary in the scene.
    r = -200.0 * length(adversaries(s))
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
