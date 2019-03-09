using Plots, Colors, PyCall, Compose, LazySets
include("adjacency_generator.jl")


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
    inf_visibile(G,obs)

Takes in MetaGraph and array of obstacles. obs needs to have array of the Hyperrectangles
for the obstacles
"""
function inf_visibile(G,obs)
    V_x = []
    for i in 1:nv(G)
        V_x_i = inf_visible_node(G, obs, i)
        push!(V_x, collect(V_x_i))
    end
    return V_x
end

function inf_visible_node(G, obs, node)
    V_x_i = []
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


function plot_visible_node(G, V_x_i, obs=[], nodenum=[])
    finalplot = graph2gridplot(G,obs)
    if !isempty(nodenum)
        scatter!(finalplot, Tuple(G[nodenum,:pos]), color="blue", markersize=20)
    end
    for pos in V_x_i
        scatter!(finalplot, Tuple(pos), color="red", markersize=20)
    end
    finalplot
    png(finalplot, "test.png")
end
