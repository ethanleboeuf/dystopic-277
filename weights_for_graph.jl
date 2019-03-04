using Convex, SCS, LightGraphs, LinearAlgebra

"""
    find_weights(A)

Takes in Adjacency matrix and gives weights of edges, P
"""
function find_weights(A)
    # the maximal value in sources and destinations is the number of nodes
    n = maximum([size(A,1); size(A,2)])
    idx = findall(x-> (x==0),A)
    P = Semidefinite(n,n)
    constraint = P*ones(n,1) == ones(n,1)
    problem = minimize(opnorm(P - (1/n)*ones(n,n)))
    problem.constraints += constraint
    # Need to constrain the values of P to be 0 at the same indices where A = 0
    # This seems to be where the problem becomes "Infeasible". Not sure how to
    # fix it.
    
    for ii = 1:length(idx)
        xpos = idx[ii][1]
        ypos = idx[ii][2]
        constraint3 = P[xpos,ypos]==0
        problem.constraints += constraint3
    end
    solve!(problem, SCSSolver())
    Pval = evaluate(P)
    return P
end
