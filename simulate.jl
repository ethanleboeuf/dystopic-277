using ARDESPOT, ParticleFilters, BeliefUpdaters, POMDPSimulators




n_particles = 1000


p = SkyNet()


ag0 = Agent.(p.partition_points, p.partitions)

b0 = [SNState(ag0, [rand_pos(graph(p))]) for i in 1:10]
s0 = rand(b0)

## MAYBE MAKE OUR OWN
belief_updater = SIRParticleFilter(p, n_particles)

solver = DESPOTSolver(max_trials = 1) # has hecka keywords
policy = solve(solver, p)

history = simulate(HistoryRecorder(max_steps = 10), p, policy, belief_updater, b0, s0)