using LightGraphs, MetaGraphs, Compose, LinearAlgebra
using DataFrames, GraphPlot, Colors
include("weights_for_graph.jl")
"""
    base_adjacency_grid

Creates a basic grid of m x n size (m and n must be ints)
"""
function base_adjacency_grid(m,n)
    all_nodes = []
    for ii = 1:m
        for jj = 1:n
            push!(all_nodes,[ii,jj])
        end
    end
    return all_nodes
end
"""
    find_neighbors

Finds neighbors for given node given full environment
"""
function find_neighbors(node,all_nodes)
    dirs = [[1,0],[0,1],[-1,0],[0,-1],[1,1],[1,-1],[-1,1],[-1,-1]]
    result = []
    for dir in dirs
        poss_neighbor = [node[1]+dir[1],node[2]+dir[2]]
        if poss_neighbor in all_nodes
            push!(result,poss_neighbor)
        end
    end
    return result
end

"""
    create_adjacency

Given size m and n, will create a grid adjacency matrix
"""
function create_adjacency(m,n)
    all_nodes = base_adjacency_grid(m,n)
    g = Graph(m*n)
    for node in enumerate(all_nodes)
        neighborset = find_neighbors(node[2],all_nodes)
        for neighbori in neighborset
            dst_pos = findall(x-> x==neighbori,all_nodes)
            add_edge!(g,node[1],dst_pos[1])
        end
    end
    A = adjacency_matrix(g)
    return A
end

"""
    create_weighted_graph

Reads an adjacency matrix and outputs the MetaGraph associated with the optimal
weights.
"""
function create_weighted_graph(A)
    Aw = find_weights(A)
    mg = MetaGraph(A)

    for edge in edges(mg)
        plc = [src(edge),dst(edge)]
        set_prop!(mg,edge,:weight,Aw[plc[1],plc[2]])
    end

    return Aw,mg
end

"""
    read_adjacency

Reads a weighted sparse adjacency matrix and outputs the meta graph.
"""
function read_adjacency(As)

    rows = rowvals(As)
    vals = nonzeros(As)
    m, n = size(As)
    srcs = Vector{Int}(undef, n)
    dsts = Vector{Int}(undef, n)
    wgts = Vector{Float64}(undef, n)
    for i = 1:n
       for j in nzrange(As, i)
          row = rows[j]
          val = vals[j]
          # perform sparse wizardry...
       end
    end
    println()
    for (i, line) in enumerate(filename)
        println(line)
        spl = split(line, ",")
        srcs[i] = parse(Int64, spl[1])
        dsts[i] = parse(Int64, spl[2])
        wgts[i] = parse(Float64, spl[3])
    end
    return construct_graph(srcs, dsts, wgts)
end

"""
    construct_graph(src, dst, wgt)

Takes in three vectors (sources, destinations, weights) and constructs a MetaGraph out of them.
"""
function construct_graph(srcs, dsts, wgts)
    # the maximal value in sources and destinations is the number of nodes
    n = maximum([srcs; dsts])
    G = MetaGraph(SimpleGraph(n))
    defaultweight!(G, 0.0)
    # For each src=>dst, wgt, create a weighted edge for it
    for (s,d,w) in zip(srcs, dsts, wgts)
        # ignore zero weight edges, since the default is 0 anyway
        if w > 0
            add_edge!(G, s, d)
            set_prop!(G, s, d, :weight, w)
        end
    end
    return G
end

"""
    plot_graph(G, fn)

Save the graph `G` as an SVG in the local directory with the filename `fn`.
Uses the default `spring_layout` layout.
"""
function plot_graph(G, fn)
    mn, mx = extrema(weights(G))
    line_colors = map(edges(G)) do ε
        src, dst = Tuple(ε)
        weighted_color_mean(weights(G)[src, dst], colorant"white", colorant"black")
    end

    nodelabel = 1:nv(G)
    p = gplot(G, edgestrokec = line_colors, nodelabel = nodelabel)
    draw(SVG(fn, 10cm, 10cm), p)
end
