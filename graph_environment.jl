using LightGraphs, MetaGraphs, Compose
using DataFrames, GraphPlot, Colors

"""
    read_graph

Reads a comma delimitted text file where each line represents: `source, destination, weight` for each edge.
"""
function read_graph(filename)
    f = open(filename)
    s = readlines(f)
    srcs = Vector{Int}(undef, length(s))
    dsts = Vector{Int}(undef, length(s))
    wgts = Vector{Float64}(undef, length(s))
    for (i, line) in enumerate(s)
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


"""
    main

I don't believe in the concept of "main" functions, but whatever.
"""
function main()
    G = read_graph("weightsSmall.csv")
    plot_graph(G, "test-graph.svg")
end