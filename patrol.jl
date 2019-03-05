using Plots, Colors, PyCall, Compose
pyplot()
include("graph_environment.jl")
# This code simulates a patrol of a grid where expected commute time
# averaged over all pairs of nodes is minimized.

"""
    find_path(time,startNode,mg)

Takes in amount of time to simulate, t, the initial node, node0, and the
MetaGraph, mg, and outputs the calculated travel path
"""
function find_path(time_end,startNode,mg)

    M = mg
    n = nv(M)
    k = sqrt(n)
    # path in node form
    Path = zeros(time_end,1);
    Path[1] = startNode;
    for t = 2:time_end
        x = rand();
        num_neighbors = length(neighbors(M,Int(Path[t-1])))
        prob_choices = zeros(num_neighbors,1)
        for neighbor in enumerate(neighbors(M,Int(Path[t-1])))
            prob_choices[neighbor[1]] = get_prop(M,Int(Path[t-1]),neighbor[2],:weight)
            if x < sum(prob_choices[1:neighbor[1]])
                Path[t] = neighbor[2]
                break
            end
        end
    end

    # path in x and y coordinates
    travel = zeros(time_end,2);
    for i = 1:length(Path)
        if mod(Path[i],k) != 0
            travel[i,1] = k*(Path[i]/k - floor(Path[i]/k));
            travel[i,2] = ceil(Path[i]/k)
        else
            travel[i,1]  = k;
            travel[i,2] = ceil(Path[i]/k)
        end
    end
    return travel
end
"""
    plot_path(travel,fn,mg)

Plots the travel path given (also needs the overall MetaGraph) and saves it to
the given file location
"""
function plot_path(travel,fn,mg)
    pyplot()
    M = mg
    n = nv(M)
    k = sqrt(n)
    time_end = size(travel,1)
    finalplot = plot(legend=false)
    title!("Patrolling of Robots")
    xlabel!("x position")
    ylabel!("y position")
    plot!(collect(1:(k-1)/99:k), ones(100,1), color = :black)
    plot!(collect(1:(k-1)/99:k), ones(100,1)*k, color = :black)
    plot!(ones(100,1), collect(1:(k-1)/99:k), color = :black)
    plot!(ones(100,1)*k, collect(1:(k-1)/99:k), color = :black)
    for ii = 1:time_end-1
        line_color = weighted_color_mean(ii/time_end, colorant"black", colorant"red")
        plot!(travel[ii:ii+1,1], travel[ii:ii+1,2],color = line_color)
    end
    # annotate!(k+5,0,text("this is #5", 16, :red, :center))
    png(finalplot,fn)
    return finalplot
end

"""
    main

I don't believe in the concept of "main" functions, but whatever.
"""
function main()
    time_end = 1000
    startNode = 32
    A = []
    path = find_path(time_end,startNode,A)
    plot_path(path,"ethantest.png",A)
end
