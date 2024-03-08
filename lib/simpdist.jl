# use the Oscar package:
using Polymake
const pm=Polymake
using LinearAlgebra
using DelimitedFiles

# inequality matrix:
#INPUT: num. of prob. parameters:
#OUTPUT: matrix of inequalities for polymake:
function inequality_matrix(n_prob)
    id = Matrix{Int}(I,n_prob,n_prob); #identity
    zero = zeros(Int8,n_prob);         #column of zeros
    ineq = hcat(zero,id);              #concatenate
    return ineq
end

# normalization constraints: Outcomes in Z2:
#INPUT: n-parties, number of contexts, number of probabilities
#OUTPUT: matrix of norm. constraints.
function normalization_matrix(z,n,n_prob, n_context)
    init = zeros(Int8, n_context, n_prob)
    col = -ones(Int8,n_context)
    norm = hcat(col,init)
    for i in range(1,n_context)
        norm[i,(((z^n)*(i-1))+2):(((z^n)*(i-1))+((z^n)+1))] .= 1
    end
    return norm
end

# NS matrix with outcomes on Zd:
# marginalize to single measurements.
#INPUT: d-outcomes, n-parties, incidence array, number of prob.
#OUTPUT: ns constraint matrix
function ns_matrix(d,n,in_array, n_prob)
    init = zeros(Int8, length(in_array), n_prob+1) # create matrix of zeros
    for i in range(1,length(in_array))           # run through each ns constraint
        ns_eq = in_array[i];                     
        for j in range(1,length(ns_eq))          # run through each element of lst
            if j < (((length(ns_eq))/2)+1)
                init[i,ns_eq[j]+1] = 1;          # first two elements are +1
            else
                init[i,ns_eq[j]+1] = -1;         # second two elements are -1
            end
        end
    end
    return init
end

# n-cycle NS conditions:
#NS conditions: Outcomes in Zd: n contexts:
#INPUT: n-number of contexts, d-outcomes
#OUTPUT: array of n(d-1) NS conditions
function ncycle_array(d,n)

    n_prob = (d^2)*n # number of probabilities
    num_ns = n*(d-1) # number of ns conditions
    
    lst = [];
    
    # context for loop:
    for i in range(0,n-1)
        # bad if statement placement:
        if (i<(n))
            # marginal outcome for loop
            for j in range(0,d-2)
                in_lst = [];
                out_lst = [];
                # outcomes consistent with marginal:
                for k in range(1,d)
                    push!(in_lst,(i*(d^2)+j*d+k));
                    # if last context, use periodic boundary:
                    if (i<(n-1))
                        out = (i+1)*(d^2)+((k-1)*d)+(j+1);
                    else
                        out = ((i+1)*(d^2)+((k-1)*d)+(j+1)) % (n*(d^2));
                    end
                    push!(out_lst,out);
                end
            in_lst = vcat(in_lst,out_lst)
            push!(lst,in_lst);
            end
        end
    end
    return lst
end

#INPUT: incidence array, N-number prob, C-number context
#OUTPUT: polymake polytope object
function polymake_polytope(d,n,in_array,N,C)
    ineq = inequality_matrix(N);
    norm = normalization_matrix(d,n,N, C);
    ns   = ns_matrix(d,n,in_array, N);
    eq = vcat(norm,ns);
    p = pm.polytope.Polytope(INEQUALITIES=ineq, EQUATIONS = eq)
    return p
end

#INPUT: Polymake polytope object:
#OUTPUT: representative vertices
function representative_vertices(polytope)
    pm.polytope.combinatorial_symmetries(polytope)
    pm.give(polytope, "GROUP.REPRESENTATIVE_VERTICES");
    return polytope.GROUP.REPRESENTATIVE_VERTICES
end

# convert to probability tables:
function expect_to_pr(vert_array, proj_array)
    pr_array = [];              # empty array to fill
    for i in 1:size(V_L)[1]     # for each vertex
        v_array = Float64[];    # create vertex with prob.
        for j in 1:60           # for each projector
            # calculate inner product of vertex and projector:
            vj = convert(Float64,(transpose(V_L[i,:])*proj_array[j,:]))/4;
            # add entry to vector
            v_array = push!(v_array,vj);
        end
        # add vector to pr_array:
        pr_array = push!(pr_array,transpose(v_array));
    end
    # return pr_array
    return pr_array
end

# input: scenario: d outcomes, n parties, NS conditions, number of probabilities, number of contexts:
# output: bell polytope:
function bell_polytope(d,n,inarray,N,C)
    # generate NS polytope:
    NS = polymake_polytope(d,n,in_array,N,C);
    # isolate deterministic vertices;
    det_idx = [];           # index for deterministic vertices:
    for i in 1:NS.N_VERTICES
        z = length(findall(y->y ==1,transpose(NS.VERTICES[i,:])));
        if z == (C+1)
            push!(det_idx,i)
        end
    end
    # deterministic vertices:
    DET = NS.VERTICES[det_idx,:];
    # generate bell polytope:
    B = pm.polytope.Polytope(POINTS=DET);
    return B
end

# input: NS polytope,
# output: bell polytope:
function bell_from_NS(NS,D)
    # isolate deterministic vertices;
    det_idx = []; D = length(NS.VERTICES[1,:])-1;          # index for deterministic vertices:
    for i in 1:NS.N_VERTICES
        z = length(findall(y->(y ==1) || (y ==-1),transpose(NS.VERTICES[i,:])));
        if z == (D+1)
            push!(det_idx,i)
        end
    end
    # deterministic vertices:
    DET = NS.VERTICES[det_idx,:];
    # generate bell polytope:
    B = pm.polytope.Polytope(POINTS=DET);
    return B
end

# note: bell polytope when given in expectation coordinates:
# input: scenario: d outcomes, n parties, NS conditions, number of probabilities, number of contexts:
# output: bell polytope:
function bell_polytope_exp(C,D,T)
    # generate NS polytope:
    NS_ineq = ns_ineq_matrix(C,D,T);
    NS = pm.polytope.Polytope(INEQUALITIES = NS_ineq);
    # isolate deterministic vertices;
    det_idx = [];           # index for deterministic vertices:
    for i in 1:NS.N_VERTICES
        z_p = length(findall(y->y ==1,transpose(NS.VERTICES[i,:])));
        z_m = length(findall(y->y ==-1,transpose(NS.VERTICES[i,:])));
        if (z_p+z_m) == (D+1)
            push!(det_idx,i)
        end
    end
    # deterministic vertices:
    DET = NS.VERTICES[det_idx,:];
    # generate bell polytope:
    B = pm.polytope.Polytope(POINTS=DET);
    return B
end

#INPUT: Polymake polytope object:
#OUTPUT: representative facets
function representative_facets(polytope)
    pm.polytope.combinatorial_symmetries(polytope)
    pm.give(polytope, "GROUP.REPRESENTATIVE_FACETS");
    return polytope.GROUP.REPRESENTATIVE_FACETS
end

# input: polytope, Int (0/1): 0/1 based indexing:
# output: array of adjacencies of vertices
function adj_to_array(P,b)
    # create adjacency data:
    A = P.GRAPH.ADJACENCY;
    # convert to string:
    A = string(A);
    # assume we know some features of output: remove extra data:
    rmv = ["pm::graph::Graph<pm::graph::Undirected>","\n"];
    for i in 1:length(rmv)
        A = replace(A,rmv[i] => "");
    end
    # split based on "{":
    A_array = split(A,"{");
    A_array = A_array[2:end];
    # remove right bracket
    for i in 1:size(A_array)[1]
        A_array[i] = replace(A_array[i],"}" => "");
    end
    # create array of neighbors for each vertex: 1-based indexing:
    ADJ = [];
    for i in 1:length(A_array)
        a = [parse(Int64,x)+b for x in split(A_array[i]," ")];
        insert!(a,1,(i-(1-b)));
        push!(ADJ,a);
    end
    # return ADJ:
    return ADJ
end

# input: array of integer indices (1-based); e.g., [1,2,3,...]
# output: array of integer indices (0-based); [0,1,2,...]
function one_to_zero_idx(array)
    for i in 1:length(array)
        array[i] = array[i]-1;
    end
    return array
end

# input: polytope (P); b = Int (0/1): whether to use 0/1 based indexing:
# input: strg: for file name (saves into current directory)
# output: text file.
function save_vertices(P,b,strg)
    # save vertices:
    idx = [(x-b) for x in 1:size(P.VERTICES)[1]];
    V = hcat((idx),P.VERTICES);
    # save to file:
    writedlm(strg*".txt", V)
end


# Note: Assume for bipartite scenarios in Z2. beta = 0 assumed.
# Note: Ineq. given as: 1+(-1)^a*T[i][1]+(-1)^b*T[i][2]+(-1)^c*T[i][3] >= 0.
# input: C: number of contexts (Int): D: dimension of polytope (number of edges) (Int)
# input: T: the coordinate in each triangle (array of arrays); e.g., [[1,2,3],[4,5,6], etc.]
# output: matrix array that encodes NS and normalization.
function ns_ineq_matrix(C,D,T)
    N = C*(2^2);
    ineq = zeros(Int8,N,dim);           #array of zeros
    col = ones(Int8,N);                 #column of ones 
    ineq = hcat(col,ineq);              # merge arrays
    beta_array = [T[i][end] for i in 1:length(T)];
    # account for homogenizing coordinate:
    X = [[T[x][y]+1 for y in 1:(length(T[x])-1)] for x in 1:length(T)];
    for i in 1:C
        beta = beta_array[i];
        for a in 0:1                    #Z2 outcomes
            for b in 0:1                #Z2 outcomes
                x = 4*i-3+2*a+b;        #ordering of probabilities
                c = mod(a+b+beta,2);         # mod2 arithmetic for XOR
                # assign values to array:
                ineq[x,X[i]] = Vector{Int8}([(-1)^a,(-1)^b,(-1)^c]);
            end
        end
    end
    return ineq
end

# input: target vertex (numeric vector array)(row vector), V: array of vertices (numeric matrix array):
# note: vertices as row vectors
# output: index
function find_vertex(target,V)
    #initialize index: with 1-based indexing, 0 is not a vertex:
    idx = 0;        # can choose any other index not in V:
    for i in 1:size(V)[1]
        if target == (V[i,:])
            idx = i
        end
    end
    return idx
end

# Polymake uses 0-based indexing, but GAP and Julia use 1-based indexing:
to_one_based_indexing(n::Number) = n + one(n)
to_zero_based_indexing(n::Number) = (n > zero(n) ? n - one(n) : throw(ArgumentError("Can't use negative index")))

# Not sure what this does (yet)
for f in [:to_one_based_indexing, :to_zero_based_indexing]
    @eval begin
        $f(itr) = $f.(itr)
        $f(s::S) where S<:AbstractSet = S($f.(s))
    end
end

#input: polymake polytope object P
#output: GAP group
function combinatorial_automorphism_group(P)
    G = pm.group.automorphism_group(P.VERTICES_IN_FACETS)
    gens_polymake = G.PERMUTATION_ACTION.GENERATORS # acting on the facets
    gens_julia = Vector{Int64}.(pm.to_one_based_indexing(gens_polymake))
    gens_gap = GAP.Globals.PermList.(GAP.julia_to_gap.(gens_julia))
    return GAP.Globals.Group(gens_gap...)
end


# input: polymake polytope:
# output: GAP group for action on VERTICES:
function vertex_action_GAP(P)
    # embed group attribute in P:
    pm.polytope.combinatorial_symmetries(P);
    pm.give(P, "GROUP");
    # permutation group acting on vertices:
    gens_poly = P.GROUP.VERTICES_ACTION.GENERATORS;
    # convert to julia array (1-based indexing):
    gens_julia = Vector{Int64}.(pm.to_one_based_indexing(gens_poly));
    # create GAP generators:
    gens_gap = GAP.Globals.PermList.(GAP.julia_to_gap.(gens_julia));
    return GAP.Globals.Group(gens_gap...)
end

# input: polymake polytope:
# output: GAP group for action on FACETS:
function facets_action_GAP(P)
    # embed group attribute in P:
    pm.polytope.combinatorial_symmetries(P);
    pm.give(P, "GROUP");
    # permutation group acting on vertices:
    gens_poly = P.GROUP.FACETS_ACTION.GENERATORS;
    # convert to julia array (1-based indexing):
    gens_julia = Vector{Int64}.(pm.to_one_based_indexing(gens_poly));
    # create GAP generators:
    gens_gap = GAP.Globals.PermList.(GAP.julia_to_gap.(gens_julia));
    return GAP.Globals.Group(gens_gap...)
end

#input: polymake polytope P:
#output: array of arrays: each element is index of orbit under action of Aut(P):
function vertex_orbit(P)
    # create GAP group:
    G = vertex_action_GAP(P)
    # use polymake to generate representative vertices of P:
    rep = representative_vertices(P);
    # find vertex identifiers: representative vertices:
    rep_idx = [find_vertex(rep[i,:],P.VERTICES) for i in 1:size(rep)[1]];
    # create array of arrays: each array is an orbit:
    orbit_arry = [Vector{Int64}(GAP.Globals.Orbit(G,i)) for i in rep_idx];

    return orbit_arry
end

#input: polymake polytope P:
#output: array of arrays: each element is index of orbit under action of Aut(P):
function facets_orbit(P)
    # create GAP group:
    G = facets_action_GAP(P)
    # use polymake to generate representative vertices of P:
    rep = representative_facets(P);
    # find vertex identifiers: representative vertices:
    rep_idx = [find_vertex(rep[i,:],P.FACETS) for i in 1:size(rep)[1]];
    # create array of arrays: each array is an orbit:
    orbit_arry = [Vector{Int64}(GAP.Globals.Orbit(G,i)) for i in rep_idx];

    return orbit_arry
end