



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

Vec2 = SVector{2, Int64}

mutable struct SkyNet <: POMDP{SNState, SNAction, SNObs}
    env::Environment
    partition_points::Vector{Vector{Tuple}}
    partitions::Vector{Region}
    max_adversaries::Int8
    intrusion_prob::Float64
end

struct Environment
    size::NTuple{2, Int}
    obstacles::Vector{Hyperrectangle}
    visibility::Vector{Vector{Tuple}}
    G::MetaGraph
end

graph(env::Environment) = env.G
graph(p::SkyNet) = p.env.G
obstacles(env::Environment) = env.obstacles
obstacles(p::SkyNet) = p.env.obstacles


"""
Sample a random position from a metagraph
"""
rand(G::AbstractGraph, rng::AbstractRNG = Random.GLOBAL_RNG) = rand(rng, 1:nv(G))



abstract struct Agent end

mutable struct Patroller <: Agent
    pos::NTuple{2, Int}
    R::Int64
    partition::Region
    state::Symbol # consider making enum
end

struct Adversary <: Agent
    pos::NTuple{2, Int}
    R::Int64
end

struct Region
    center::Tuple # maybe this one isn't necessary
    interior::Vector{Tuple}
end

struct SNState
    patrollers::Vector{Vec2}
    adversaries::Vector{Vec2}
end

agents(s::SNState) = (s.patrollers, s.adversaries)

# Just an alias:
SNAction = Vector{Vec2}

function actions(p::SkyNet, s::SNState)
    obstacles = p.env.obstacles

    A = hcat(actions(p, pat) for pat in patrollers(s))'
end

function actions(p::SkyNet, s::Vec2)
    # consider making global
    dirs = [Vec2(i, j) for i in -1:1 for j in -1:1 if !(i==j==0)]

    A = eltype(action_type(SkyNet))[]
    for dir in dirs
        n_p = dir + s
        if n_p ∈ s.partition && n_p ∈ p.env.G && all(n_p .∉ p.env.obstacles)
            push!(A, dir)
        end
    end
    return A
end

struct SNObs
    adversaries::Vector{Vec2}
end





Base.size(env::Environment) = env.size


function compute_visitibility!(env::Environment)

end


function generate_s(p::SkyNet, s::SNState, a::Vector{Vec2}, rng::AbstractRNG)
    patrollers, adversaries = agents(s)

    updated_patrollers = patrollers .+ a

    # With specified probability, create an adversary if there isn't one already.
    if length(adversaries) < p.max_adversaries
        if rand(rng) > p.intrusion_prob
            push!(adversaries, Vec2(get_pos_tup(G, rand(G, rng))))
        end
    else
        # update the adversaries' positions somehow
    end
end


function generate_o(p::SkyNet, s::SNState, rng::AbstractRNG)
    noise = rand(p.d, rng)


end