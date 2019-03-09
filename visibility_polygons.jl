using Plots, Colors, PyCall, Compose, LazySets, Polyhedra
include("adjacency_generator.jl")
gr()

# Used to convert obstacles to grid
function set2grid(set::Hyperrectangle)
    l = Int.(low(set))
    u = Int.(high(set))
    grid = base_grid((u-l .+ 1)...)
    for point in grid
        grid[point] = point + CartesianIndex(Tuple(l-[1,1]))
    end
    return grid
end

"""
    inf_visible(G,obs)

Takes in MetaGraph and array of obstacles. obs needs to have array of the Hyperrectangles
for the obstacles
"""
# Used to get all of the visible cells from every position assuming infinite views
function inf_visible(G,obs)
    V_x = []
    for i in 1:nv(G)
        V_x_i = inf_visible_node(G, obs, i)
        push!(V_x, collect(V_x_i))
    end
    return V_x
end

# Used to get all of the visible cells from a single node assuming infinite views
function inf_visible_node(G, obs, node)
    V_x_i = []
    # Tomer, I'm not sure how to not need these next ten lines. I'm trying to
    # allow for both cartesian index and LI to be inputted but need to do a check.
    # I end up repeating this code in some fashion a bunch throughout this .jl
    if typeof(node) == CartesianIndex{2}
        i_pos = collect(Float64.(Tuple(node)))
        i = G[node, :pos]
    elseif typeof(node) == Int64
        i_pos = collect(Float64.(Tuple(get_prop(G, node, :pos))))
        i = node
    end

    i_to_p = [(LineSegment(i_pos, collect(Float64.(Tuple(get_prop(G, p, :pos))))), p) for p in 1:nv(G) if !(i==p)]
    for (vect, p) in i_to_p
        inter_ob = 0
        num_ob_inter = 0
        p_pos = collect(Float64.(Tuple(get_prop(G, p, :pos))))
        for (k,ob) in enumerate(obs)
            if !is_intersection_empty(vect, ob) && length(vertices_list(intersection(vect, ob))) != 1
                inter_ob = 1
                break
            elseif !is_intersection_empty(Singleton(p_pos), ob)
                inter_ob = 1
                break
            end
            if !is_intersection_empty(vect, ob) && length(vertices_list(intersection(vect, ob))) == 1
                num_ob_inter += 1
            end
        end
        if inter_ob == 0 && num_ob_inter < 2
            push!(V_x_i,get_prop(G, p, :pos))
        end
    end
    return V_x_i
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
function limited_visible(G, obs, R, nodes=[])
    L_V_x = []
    for i in 1:nv(G)
        L_V_x_i = limited_visibility_node(G, obs, i, R)
        push!(L_V_x, collect(L_V_x_i))
    end
    return L_V_x
end

# Used to get all of the visible cells from a single node assuming limited views
function limited_visibility_node(G, obs, node, R)
    V_x = inf_visible_node(G, obs, node)
    L_V_x_i = []

    if typeof(node) != CartesianIndex{2}
        node = get_prop(G, node, :pos)
    end
    [push!(L_V_x_i, pos) for pos in V_x if (norm(Tuple(pos - node)) <= R) ]
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
        scatter!(finalplot, Tuple(pos), color="red", markersize=5)
    end
    for node in nodes
        scatter!(finalplot, Tuple(G[node,:pos]), color="blue", markersize=5)
        if !isempty(R)
            plot!(finalplot, Ball2(Float64.(collect(Tuple(G[node,:pos]))), Float64.(R)),
                    color=:auto, 1e-3, aspectratio=1, opacity=0.5)
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
    nodes = [CartesianIndex((2, 3)),CartesianIndex((9, 2)),CartesianIndex((5,15)),
            CartesianIndex((14,11))]
    @time L_V_x = limited_visible(G, obs, R)
    T_L_V_x, n_view = total_coverage(G, L_V_x, nodes, obs)
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
