using Plots, Colors, PyCall, Compose,LazySets
pyplot()
include("adjacency_generator.jl")


function set2grid(set::Hyperrectangle)
    l = low(set)
    u = high(set)
    grid = base_grid((u-l .+ 1)...)
    # grid = map(x-> x.-[1,1],grid)
    return grid
end


function grid_minus_obstacles(Env::Hyperrectangle,obs)
    grid_env = set2grid(Env)
    grid_obs = []
    for ob in obs
        push!(grid_obs,set2grid(ob))
    end
end

"""
    visible2(node,all_nodes,R)

Takes in the current node, all nodes, and the sensor radius (integer) and outputs the
visible nodes from the current node, V_x
"""
function visibile2(node,all_nodes,R)
    all_nodes = base_adjacency_grid(5,5) #FOR TESTING ONLY. Will just give A later once
    # obstacles are integrated
    pyplot()
    Env = HyperRectangle([3 3],[3 3])
    finalplot = plot(Env)
    gui()
    V_x = []
    outer_x = []
    dirs_card = [[1,0],[0,1],[-1,0],[0,-1]]
    png(finalplot,"test_visible.png")
    return V_x,outer_x
end



"""
    visible(node,all_nodes,R)

Takes in the current node, all nodes, and the sensor radius (integer) and outputs the
visible nodes from the current node, V_x
"""
function visibile(node,all_nodes,R)
    all_nodes = base_adjacency_grid(5,5) #FOR TESTING ONLY. Will just give A later once
    # obstacles are integrated
    V_x = []
    outer_x = []
    dirs_card = [[1,0],[0,1],[-1,0],[0,-1]]
    for dir in dirs_card
        if length(V_x)>0
            if outer_x == []
                push!(outer_x,V_x[end])
            elseif outer_x[end] != V_x[end]
                push!(outer_x,V_x[end])
            end
        end
        for ii = 1:R
            poss_visible_node = [node[1]+dir[1]*ii,node[2]+dir[2]*ii]
            if poss_visible_node in all_nodes
                push!(V_x,poss_visible_node)
            else
                break
            end
        end
    end
    push!(outer_x,V_x[end])

    dir_iter = [[0,1],[0,-1]]
    for jj = -floor(R/2):floor(R/2)
        if outer_x[end] != V_x[end]
            push!(outer_x,V_x[end])
        end
        for dir in dir_iter
            if outer_x[end] != V_x[end]
                push!(outer_x,V_x[end])
            end
            for ii = 1:R
                poss_visible_node = [node[1]+jj+dir[1]*ii,node[2]+dir[2]*ii]
                if poss_visible_node in all_nodes
                    push!(V_x,poss_visible_node)
                else
                    break
                end
            end
        end
    end
    if outer_x[end] != V_x[end]
        push!(outer_x,V_x[end])
    end
    dir_iter = [[1,0],[-1,0]]
    for jj = -floor(R/2):floor(R/2)
        if outer_x[end] != V_x[end]
            push!(outer_x,V_x[end])
        end
        for dir in dir_iter
            if outer_x[end] != V_x[end]
                push!(outer_x,V_x[end])
            end
            for ii = 1:R
                poss_visible_node = [node[1]+dir[1]*ii,node[2]+jj+dir[2]*ii]
                if poss_visible_node in all_nodes
                    push!(V_x,poss_visible_node)
                else
                    break
                end
            end
        end
        if outer_x[end] != V_x[end]
            push!(outer_x,V_x[end])
        end
    end
    return V_x,outer_x
end

"""
    all_visible(A,R)

Takes in adjacency_matrix represnting the grid and outputs all of the sets of
visible cells
"""
function all_visible(A,R)
    all_nodes = base_adjacency_grid(5,5) #FOR TESTING ONLY. Will just give A later once
    n = length(all_nodes)
    V_xall = [] #wasn't sure how to initialize the vector of arrays
    for node in all_nodes
        vis,out = visibile(node,all_nodes,R)
        push!(V_xall,vis)
    end
    return V_xall
end

"""
    plot_visible(all_nodes,node)


"""
function plot_visible(all_nodes,node)
    all_nodes = base_adjacency_grid(7,7)
    V_xall = all_visible([],2) #FOR TESTING ONLY. Will just give A later once
    # println(V_xall)
    pyplot()

    dst_pos = findall(x-> x==node,all_nodes)
    finalplot = plot()
    x0 = all_nodes[1][1]
    xf = all_nodes[end][1]
    y0 = all_nodes[1][2]
    yf = all_nodes[end][2]
    plot!(Shape([x0,x0,x0,x0]+[0,xf-x0,xf-x0,0],[y0,y0,y0,y0]+[0,0,yf-y0,yf-y0])
    ,opacity=0.2)

    vert = []
    visnodes,out = visibile(node,all_nodes,2)
    for vis in visnodes
        push!(vert,[vis[1],vis[2]])
        scatter!([vis[1]],[vis[2]])
    end
    plot!(create_polygon(out),opacity=0.5)
    ndx = node[1]*length(vert)
    ndy = node[2]*ones(4)
    # plot!(Shape([x0,x0,x0,x0]+[0,xf-x0,xf-x0,0],[y0,y0,y0,y0]+[0,0,yf-y0,yf-y0]))
    png(finalplot,"test.png")
    return out
end

"""
    create_polygon(outer_nodes)


"""
function create_polygon(outer_nodes)
    xcoord = []
    ycoord = []
    for node in outer_nodes
        push!(xcoord,node[1])
        push!(ycoord,node[2])
    end
    shape = Shape(xcoord,ycoord)
    return shape
end
