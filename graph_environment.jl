using LightGraphs, MetaGraphs,Compose, Cairo, Fontconfig
using LinearAlgebra, CSV, DataFrames, GraphPlot, Colors

function main()

# Assuming all nodes have edges between of 1 meter (ignore Diagonal being longer for now)
numCW = 8; # number of discrete nodes widthwise
numCL = 8; # number of discrete nodes lengthwise

numCellsTotal = numCL*numCW
global prob_weights = CSV.read("weightsSmall.csv")
prob_weights = convert(Matrix{Float64}, prob_weights[:,1:3])
for ii = 1:size(prob_weights,1)
    prob_weights[ii,1] = round(Int,prob_weights[ii,1])
    prob_weights[ii,2] = round(Int,prob_weights[ii,2])
end

global mg = SimpleGraph(numCellsTotal)
mg = MetaGraph(mg)
cellC = 1  #initialize counter

# Adding edges to create a gridlike graph/environment
for ii = 1:numCL
    for jj = 1:numCW
        # Vertical Adjacent Cells
        if ii == 1
            add_edge!(mg,cellC,cellC+numCW)
        elseif ii == numCL
            add_edge!(mg,cellC,cellC-numCW)
        else
            add_edge!(mg,cellC,cellC+numCW)
            add_edge!(mg,cellC,cellC-numCW)
        end
        # Horizontal Adjacent Cells
        if jj == 1
            add_edge!(mg,cellC,cellC+1)
        elseif jj == numCW
            add_edge!(mg,cellC,cellC-1)
        else
            add_edge!(mg,cellC,cellC+1)
            add_edge!(mg,cellC,cellC-1)
        end
        # Diagonal Adjacent Cells
        if jj == 1
            add_edge!(mg,cellC,cellC+numCW+1)
            if ii > 1
                add_edge!(mg,cellC,cellC-numCW+1)
            end
        elseif jj == numCW
            add_edge!(mg,cellC,cellC+numCW-1)
            if ii > 1
                add_edge!(mg,cellC,cellC-numCW-1)
            end
        else
            add_edge!(mg,cellC,cellC+numCW+1)
            add_edge!(mg,cellC,cellC+numCW-1)
            if ii > 1
                add_edge!(mg,cellC,cellC-numCW-1)
                add_edge!(mg,cellC,cellC-numCW+1)
            end
        end
        cellC = cellC + 1 #iterate to next cell
    end
end
line_color = []
# Assigning weights to the edges
counter = 1;
for edge in edges(mg)
    src_node = src(edge)
    dst_node = dst(edge)
    a = [src_node dst_node]
    row_check = Int[a == [round(Int,prob_weights[ii,1]) round(Int,prob_weights[ii,2])] for ii =1:size(prob_weights,1) ]
    row_idx = findall(x->x==1,row_check)
    set_prop!(mg,edge,:weight,prob_weights[row_idx,3]) # assignging weight
    val = prob_weights[row_idx,3][1]
    # Creating vector of colors for the edges, first color is the one with
    # higher probability on the graph
    push!(line_color,weighted_color_mean(val, colorant"white", colorant"black"))
    counter = counter + 1
end

nodelabel = collect(1:cellC-1)
layout=(args...)->spring_layout(args...; C=.9)
draw(PNG("test-graph.png", 16cm, 16cm), gplot(mg,edgestrokec=line_color,nodelabel=nodelabel))
end

main()
