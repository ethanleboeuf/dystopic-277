using ARDESPOT, ParticleFilters, BeliefUpdaters, POMDPSimulators
using POMCPOW, MCTS, POMDPPolicies

include("pomdps.jl")
include("visibility_polygons.jl")

# This function overrides an ARDESPOT default behavior to
# make state-dependent actions spaces. Native support for
# that funcationality does not exist yet apparently.
actions(p::SkyNet, b::ARDESPOT.ScenarioBelief) = actions(p, b.scenarios[1][2])
function actions(p::SkyNet, b::POMCPOW.POWTreeObsNode)
    @show belief(b)
    actions(p, rand(belief(b)))
end
function MCTS.next_action(gen::RandomActionGenerator, p::SkyNet, s, snode::AbstractStateNode)
    rand(actions(p, rand(s)))
end
function POMDPPolicies.action(policy::RandomPolicy, b::Nothing)
    return rand(policy.rng, actions(policy.problem))
end




n_parts = 2000

if !isdefined(Main, :p)

    @eval begin
        println("Building SkyNet")
        n_parts_prev = n_parts
        p = SkyNet(R = 3,
                   n_agents = 3,
                   partition_points =[CartesianIndex((2, 3)),
                                       CartesianIndex((5,13)),
                                       CartesianIndex((13,7))])

        println("initializing everything else")
        ag0 = Agent.(p.partition_points, p.partitions)

        b0 = [SNState(ag0, [rand_pos(graph(p))]) for i in 1:n_parts]
        s0 = rand(b0)

        ## MAYBE MAKE OUR OWN
        belief_updater = SIRParticleFilter(p, n_parts)
    end
end


K = 600
D = 100
solver = DESPOTSolver(T_max = 10.0,
                      K = K,
                      # max_trials = 10000,
                      tree_in_info = true,
                      D = D,
                      bounds = IndependentBounds(-1e6, 1e6, check_terminal = true),
                      random_source = MemorizingSource(K, D, Random.GLOBAL_RNG, min_reserve=1000))
# solver = POMCPOWSolver() # has hecka keywords
policy = solve(solver, p)

println("simulating...")
history = simulate(HistoryRecorder(max_steps = 2), p, policy, belief_updater, b0, s0)