using Plots, Colors, PyCall, Compose,LazySets, Polyhedra
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
    ob_union = UnionSetArray(obs)
    for i in 1:nv(G)
        V_x_i = []
        i_pos = Float64.([get_prop(G, i, :pos)[1], get_prop(G, i, :pos)[2]])
        i_to_p = [(LineSegment(i_pos, Float64.([get_prop(G, p, :pos)[1], get_prop(G, p, :pos)[2]])), p) for p in 1:nv(G) if !(i==p)]
        for (vect,p) in i_to_p
            for (k,ob) in enumerate(obs)
                if is_intersection_empty(vect, ob)
                    push!(V_x_i,get_prop(G, p, :pos))
                elseif length(vertices_list(intersection(vect,ob)))==1
                    # println(i=>p)
                    push!(V_x_i,get_prop(G, p, :pos))
                else
                    # println(i=>p=>"didn't work")
                    if !isempty(V_x_i) && k > 1
                        # println(V_x_i=>k)
                        for ii = 1:k-1
                            pop!(V_x_i)
                        end
                    end
                    break
                end
                # if is_intersection_empty(vect, ob_union)
                #     push!(V_x_i,get_prop(G, p, :pos))
                # else
                #     println(intersection(vect,ob_union))
                #     if length(vertices_list(intersection(vect,ob_union)))==1
                #         push!(V_x_i,get_prop(G, p, :pos))
                #     end
                # end
            end
        end
        push!(V_x,collect(V_x_i))
    end
    return V_x
end

function plot_visible_node(G,V_x_i)
    finalplot = graph2gridplot(G)
    for pos in V_x_i
        scatter!(finalplot,[pos[1]],[pos[2]], color = "red", markersize = 20)
    end
    finalplot
    png(finalplot,"test.png")
end
#
# """
#     plot_visible(all_nodes,node)
#
#
# """
# function plot_visible(all_nodes,node)
#     all_nodes = base_adjacency_grid(7,7)
#     V_xall = all_visible([],2) #FOR TESTING ONLY. Will just give A later once
#     # println(V_xall)
#     pyplot()
#
#     dst_pos = findall(x-> x==node,all_nodes)
#     finalplot = plot()
#     x0 = all_nodes[1][1]
#     xf = all_nodes[end][1]
#     y0 = all_nodes[1][2]
#     yf = all_nodes[end][2]
#     plot!(Shape([x0,x0,x0,x0]+[0,xf-x0,xf-x0,0],[y0,y0,y0,y0]+[0,0,yf-y0,yf-y0])
#     ,opacity=0.2)
#
#     vert = []
#     visnodes,out = visibile(node,all_nodes,2)
#     for vis in visnodes
#         push!(vert,[vis[1],vis[2]])
#         scatter!([vis[1]],[vis[2]])
#     end
#     plot!(create_polygon(out),opacity=0.5)
#     ndx = node[1]*length(vert)
#     ndy = node[2]*ones(4)
#     # plot!(Shape([x0,x0,x0,x0]+[0,xf-x0,xf-x0,0],[y0,y0,y0,y0]+[0,0,yf-y0,yf-y0]))
#     png(finalplot,"test.png")
#     return out
# end
#
# """
#     create_polygon(outer_nodes)
#
#
# """
# function create_polygon(outer_nodes)
#     xcoord = []
#     ycoord = []
#     for node in outer_nodes
#         push!(xcoord,node[1])
#         push!(ycoord,node[2])
#     end
#     shape = Shape(xcoord,ycoord)
#     return shape
# end
