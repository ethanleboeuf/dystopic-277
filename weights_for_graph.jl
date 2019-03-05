using Convex, SCS, LightGraphs, LinearAlgebra

"""
    find_weights(A)

Takes in Adjacency matrix and gives weights of edges, P
"""
function find_weights(A)
    # the maximal value in sources and destinations is the number of nodes
    n = maximum([size(A,1); size(A,2)])
    idx = findall(x-> (x==0),A)
    # P = Semidefinite(n)
    P = Variable(n,n)

    constraints = [P*ones(n,1) == ones(n,1),P>=0,P==P']

    idx = findall(x -> x == 0, A)
    for ii = 1:length(idx)
        xpos = idx[ii][1]
        ypos = idx[ii][2]
        constraint3 = (P[xpos,ypos]==0)
        push!(constraints,constraint3)
    end
    problem = minimize(opnorm(P - (1/n)*ones(n,n)),constraints)
    solve!(problem, SCSSolver(verbose=false,max_iters=10^4,eps=1e-9))

    Pval = copy(P.value)
    Pval[(Pval) .< 1e-4] .= 0
    Pval = Symmetric(Pval)
    return Pval
end
