using LightGraphs, MetaGraphs, Compose, LinearAlgebra
using DataFrames, GraphPlot, Colors,Combinatorics
include("weights_for_graph.jl")
gr()
"""
    base_adjacency_grid

Creates a basic grid of m x n size (m and n must be ints)
"""
function base_grid(m, n)
    all_nodes = collect(CartesianIndices((m, n)))
    return all_nodes
end
"""
    find_neighbors

Finds neighbors for given node given full environment
"""
function find_neighbors(node, all_nodes)
    dirs = [CartesianIndex(i, j) for i in -1:1 for j in -1:1 if !(i==j==0)]
    N = []
    for dir in dirs
        n_p = dir+node
        if n_p in all_nodes
            push!(N, n_p)
        end
    end
    return N
end

"""
    create_graph

Given size m and n, will create a grid adjacency matrix
"""
function create_graph(m, n)
    V = CartesianIndices((m, n))
    G = MetaGraph(m * n)
    LI = LinearIndices(V)
    for i in V
        Nᵢ = find_neighbors(i, V)
        for j in Nᵢ
            add_edge!(G, LI[i], LI[j])
            set_prop!(G, LI[i], :pos, i)
        end
    end
    set_indexing_prop!(G, :pos)
    return G
end

function remove_obstacles(G, obs)
    obs_coord = []
    for ob in obs
        push!(obs_coord, set2grid(ob))
    end
    for coord in obs_coord
        corners = [coord[1, 1], coord[1, end], coord[end, 1], coord[end, end]]
        borders = [coord[1, :], coord[:, 1], coord[end, :], coord[:, end]]
        for pos in coord
            for neighbori in collect(neighbors(G, G[pos, :pos]))
                rem_edge!(G, G[pos, :pos], neighbori)
            end
        end
    end
    return G
end

function graph2gridplot(G, obs=[])
    finalplot = scatter(legend=false)
    for i in 1:nv(G)
        scatter!(Tuple(G[i, :pos]), color="black")
        for j in neighbors(G, i)
            xpos0, ypos0 = Tuple(G[i, :pos])
            xn, yn = Tuple(G[j, :pos])
            plot!([xpos0, xn], [ypos0, yn], color="black")
        end
    end
    [plot!(finalplot, ob) for ob in obs if !isempty(obs)]
    png(finalplot, "test.png")
    finalplot
end
"""
    create_weighted_graph

Reads an adjacency matrix and outputs the MetaGraph associated with the optimal
weights.
"""
function create_weighted_graph(G)
    A = adjacency_matrix(G)
    Aw = find_weights(adjacency_matrix(G))
    for edge in edges(G)
        plcx, plcy = [src(edge), dst(edge)]
        set_prop!(G, edge, :weight, Aw[plcx, plcy])
    end
    return G
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
        weighted_color_mean(weights(G)[src, dst], colorant"black", colorant"black")
    end

    nodelabel = 1:nv(G)
    p = gplot(G, edgestrokec = line_colors, nodelabel = nodelabel)
    draw(SVG(fn, 10cm, 10cm), p)
end
