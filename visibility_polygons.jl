using Plots, Colors, PyCall, Compose, LazySets, Polyhedra
include("adjacency_generator.jl")
gr()

# Used to convert obstacles to grid
function set2grid(set::Hyperrectangle)
    l = Tuple(Int.(low(set) .- 1))
    u = Tuple(Int.(high(set)))

    set_points = base_grid(u .- l) .+ [CartesianIndex(l)]
end

"""
    segment_blocked(::LineSegment, ::LazySet)

Returns true if the line and set intersect in a non-trivial way, i.e. at more than one point.
Naturally, this does not work on the intersection of two lines (since it is always at a point)
"""
function segment_blocked(line, obs)
    # any(iszero, radius_hyperrectangle(obs)) && error("The obstacle is infinitesimally thin, which will
    #                                                  break this function. obs radius: $(radius_hyperrectangle(obs))")

    if line.q ∈ obs || line.p ∈ obs
        return true
    end
    inter = intersection(line, obs)
    if isempty(inter) || length(vertices_list(inter)) <= 1
        return false
    else
        return true
    end
end

get_pos(G::MetaGraph, i) = collect(Float64.(Tuple(G[i, :pos])))



"""
    inf_visible(G,obs)

Takes in MetaGraph and array of obstacles. obs needs to have array of the Hyperrectangles
for the obstacles
"""
function inf_visible(G, obstacles)
    N = nv(G)
    # initialize a visibility vector for each point in the graph
    visibility = [CartesianIndex{2}[] for i in 1:N]
    for i in 1:N
        i_pos = get_pos(G, i)
        for p in (i+1):N
            p_pos = get_pos(G, p)

            vect = LineSegment(i_pos, p_pos)
            # @show vect
            isobstructed = any(segment_blocked(vect, o) for o in obstacles)
            if !isobstructed
                # If there is no obstruction,
                # then i can see p and p can see i
                push!(visibility[i], G[p, :pos])
                push!(visibility[p], G[i, :pos])
            end
        end
    end
    return visibility
end



"""
TOMERRRRRRRR

With a lot of nodes this takes a while to run (2.5 minutes when doing a 20x20)
Need to figure out how to give inputted nodes and output a vector with lots of
empty values and just the visible cells I care about in their correct positions.

Note: this would only be run once per env/obs pairing so it wouldn't be that big
of a deal.
"""
# Used to get all of the visible cells from all nodes assuming limited views
function limited_visible(G, obs, R, V_x)
    L_V_x = []
    for i in 1:nv(G)
        L_V_x_i = limited_visibility_node(G, obs, i, R, V_x)
        push!(L_V_x, collect(L_V_x_i))
    end
    return L_V_x
end

# Used to get all of the visible cells from a single node assuming limited views
function limited_visibility_node(G, obs, node, R, V_x)
    L_V_x_i = []
    if typeof(node) != CartesianIndex{2}
        V_x_i = V_x[node]
        node = get_prop(G, node, :pos)
    else
        i = G[node, :pos]
        V_x_i = V_x[node]
    end

    [push!(L_V_x_i, pos) for pos in V_x_i if (norm(Tuple(pos - node)) <= R) ]
    return L_V_x_i
end



# Takes in visible nodes (infinite or limited) and outputs all of the visible cells
# from given subset of nodes
function total_coverage(G, V_x, nodes, obs)
    T_V_x = []

    if typeof(nodes[1]) == CartesianIndex{2}
        node_idx = [G[node, :pos] for node in nodes]
    elseif typeof(nodes[1]) == Int64
        node_idx = [node for node in nodes]
    end

    for idx in node_idx
        push!(T_V_x, vcat(V_x[idx]))
    end
    T_V_x = unique(vcat(T_V_x...))

    grid_obs = []
    for ob in obs
        push!(grid_obs, vcat(set2grid(ob)...))
    end

    grid_obs = unique(vcat(grid_obs...))
    T_env = nv(G) - length(grid_obs)
    return T_V_x, T_env
end

# Plots visible cells given a Total Coverage visible set (works for one node too)
function plot_visible(G, T_V_x, nodes, obs, R = [])

    finalplot = graph2gridplot(G,obs)
    xlims0 = xlims(finalplot)
    ylims0 = ylims(finalplot)

    if typeof(nodes[1]) == CartesianIndex{2}
        nodes = [G[node, :pos] for node in nodes]
    end

    for pos in T_V_x
        scatter!(finalplot, Tuple(pos), color="lightgreen", markersize=5)
    end
    for node in nodes
        scatter!(finalplot, Tuple(G[node,:pos]), color="blue", markersize=5)
        if !isempty(R)
            plot!(finalplot, Ball2(Float64.(collect(Tuple(G[node,:pos]))), Float64.(R)),
                    color=:lightgreen, 1e-3, aspectratio=1, opacity=0.5)
        end
    end
    plot!(finalplot, xlims=xlims0, ylims=ylims0)
    png(finalplot, "test.png")

    return finalplot
end

# Example workflow. From graph and obstacles to plot. THIS MIGHT TAKE A WHILE
# RUN DON'T SAY I DIDN'T WARN YOU.
function main()
    R = 4
    obs = [Hyperrectangle(low=[3, 3], high=[6, 8]),Hyperrectangle(low=[4, 9], high=[8, 10]),
            Hyperrectangle(low=[15,15], high=[16,19]),Hyperrectangle(low=[16,3],high=[18,10]),
            Hyperrectangle(low=[3,16],high=[10,18])]
    @time G = create_graph(20, 20)
    @time G = remove_obstacles(G,obs)
    nodes = [CartesianIndex((2, 3)),CartesianIndex((9, 2)),CartesianIndex((5,13)),
            CartesianIndex((14,11)),CartesianIndex((13,7))]
    # @time V_x = inf_visible(G, obs)
    @time L_V_x = limited_visible(G, obs, R, V_x)
    T_L_V_x, n_view = total_coverage(G, L_V_x, nodes, obs)
    println(length(T_L_V_x)/n_view*100)
    fig = plot_visible(G, T_L_V_x, nodes, obs, R)
end

# I dont' think this is needed
#
# function plot_visible_node(G, V_x_i, obs=[], nodenum=[], R = [])
#     finalplot = graph2gridplot(G,obs)
#     xlims0 = xlims(finalplot)
#     ylims0 = ylims(finalplot)
#     if !isempty(R)
#         plot!(finalplot, Ball2(Float64.(collect(Tuple(G[nodenum,:pos]))), Float64.(R)),
#                 color="green", 1e-3, aspectratio=1, opacity=0.5)
#     end
#     if !isempty(nodenum)
#         scatter!(finalplot, Tuple(G[nodenum,:pos]), color="blue", markersize=20)
#     end
#
#     for pos in V_x_i
#         scatter!(finalplot, Tuple(pos), color="red", markersize=20)
#     end
#     plot!(finalplot, xlims=xlims0, ylims=ylims0)
#     finalplot
#     png(finalplot, "test.png")
#     return finalplot
# end
