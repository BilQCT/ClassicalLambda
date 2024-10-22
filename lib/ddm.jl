# Struct definition
struct Polytope
    H_representation::Matrix{Rational{Int64}} # H-representation: Inequalities (facets)
    V_representation::Matrix{Rational{Int64}} # V-representation: Vertices
end

# Constructor function for the Polytope struct
function Polytope(H_rep::Matrix{Rational{Int64}}, V_init::Union{Nothing, Matrix{Rational{Int64}}}=nothing, K::Union{Nothing, Vector{Int}}=nothing)
    if (V_init == nothing) && (K == nothing)
        # Compute initial vertex set and corresponding facets if not provided
        K, V_init = initial_DD(H_rep)
        V_rep = ddm(H_rep, K, V_init,true)
    elseif (V_init != nothing) && (K == nothing)
        V_rep = ddm(H_rep, K, V_init,true)
    else
        # Use the provided initial vertex set and facet indices
        V_rep = ddm(H_rep, K, V_init,true)
    end
    return Polytope(H_rep, V_rep)
end

# Method to print the Polytope
function Base.show(io::IO, polytope::Polytope)
    println(io, "Polytope H-representation (Facets):\n", polytope.H_representation)
    println(io, "\nPolytope V-representation (Vertices):\n", polytope.V_representation)
end

# Accessor methods for H and V representations
function get_H_representation(p::Polytope)
    return p.H_representation
end

function get_V_representation(p::Polytope)
    return p.V_representation
end

using AbstractAlgebra

# Helper functions: initial_DD, ddm, and others
function initial_DD(A)
    dim = size(A)[2]-1
    b = (-A[:,1])
    A = A[:,2:end]
    M = [i for i in 1:size(A)[1]]
    K = []
    rK = rank(A[K,:])
    while rK < dim
        i = rand(M)
        Ki = vcat(K,i)
        rKi = rank(A[Ki,:])
        if rKi > rK
            K = Ki
            rK = rKi
        end
        remove = findall(x->x==i, M)
        deleteat!(M, remove)
    end
    AK = (A[K,:])
    bK = reshape(b[K],:,1)
    SA = matrix_space(QQ, size(AK)[1], size(AK)[2])
    Sb = matrix_space(QQ, size(bK)[1], size(bK)[2])
    QAK = SA(AK)
    QbK = Sb(bK)
    R = vcat(ones(Rational{Int64},1,1), Matrix{Rational{Int64}}(solve(QAK, QbK, side=:right)))
    return Vector{Int64}(K), R
end

using DelimitedFiles

function ddm(A::Matrix{Rational{Int64}}, K::Union{Nothing, Vector{Int}}, R::Matrix{Rational{Int64}}, init::Bool=false, save::Bool=true, filename::String="filename")
    # If no initial vertex set is provided, initialize it using initial_DD
    if !init
        K, R = initial_DD(A)
    end

    # Input size
    m = size(A)[1]
    d = size(A)[2]

    # Create vertex matrix to iterate over
    Riter = copy(R)

    # Remaining inequalities
    if K==nothing
        Aiter = copy(A)
        M = Vector{Int64}([i for i in 1:m])
    else
        Aiter = copy(A[K, :])
        M = Vector{Int64}([i for i in 1:m if i ∉ K])
    end

    for i in M
        # Define the facet inserted at step i
        hi = A[i, :]

        # Initialize vertices generated at step i
        Ri = Matrix{Rational{Int64}}(undef, d, 0)

        # Segment Riter into H+, H0, H-
        rp, rt, rn = regions(hi, Riter)

        # For all vertices in H-
        for j in rn
            # Current vertex
            Rj = Riter[:, j]

            # Inner product
            a = Rational{Int64}(transpose(hi) * Rj)

            # All neighbors of j such that j is in rp
            Njp = adjacency(j, rp, Aiter, Riter)

            # New vertices
            for k in Njp
                Rjp = Riter[:, k]

                # Inner product
                b = Rational{Int64}(transpose(hi) * Rjp)

                if (b - a) != 0
                    p = b // (b - a)
                    if (p > 0)
                        # Generate new point
                        Rjk = p * Rj + (1 - p) * Rjp

                        # Add to representation matrix
                        Ri = hcat(Ri, Rjk)
                    end
                end
            end
        end

        # Save new vertices with every iteration
        if save
            open(filename, "w") do io
                writedlm(io, Riter, ',')
            end
        end

        # Update inequalities
        Aiter = vcat(Aiter, transpose(hi))

        # Remove violating extreme rays
        Keep = [i for i in 1:size(Riter)[2] if i ∉ rn]
        Riter = Riter[:, Keep]

        # Add newly generated rays
        Riter = hcat(Riter, Ri)

        # Online updates
        println("Inequality: ", i, "\n")
        println("Intermediate DD vertices: ", size(Riter), "\n")
    end

    return Riter
end


# Additional helper functions (tight_inequalities, adjacency, regions)
function tight_inequalities(x, inequalities)
    x = Vector{Rational{Int64}}(x)
    A = Matrix{Rational{Int64}}(inequalities)
    b = A * x
    return findall(y -> y == 0, b)
end

function adjacency(index, Np, A, R)
    Adj = []
    dim = size(R)[1]
    v = R[:,index]
    Zv = tight_inequalities(v, A)
    for i in Np
        u = R[:,i]
        Zu = tight_inequalities(u, A)
        Z = collect(intersect(Set(Zv), Set(Zu)))
        if length(Z) < dim - 2
            continue
        else
            if i != index
                a = rank(A[Z, :])
                if a == (dim - 2)
                    push!(Adj, i)
                end
            end
        end
    end
    return Adj
end

function regions(h, R)
    hR = transpose(R) * h
    rp = findall(x -> x > 0, hR)
    rt = findall(x -> x == 0, hR)
    rn = findall(x -> x < 0, hR)
    return rp, rt, rn
end





"""
EXTRA

function initial_DD(A)
    dim = size(A)[2]-1;
    b = (-A[:,1]); A = A[:,2:end];
    M = [i for i in 1:size(A)[1]];
    K = []; rK = rank(A[K,:]);
    while rK < dim
        i = rand(M); Ki = vcat(K,i); 
        rKi = rank(A[Ki,:]);

        # check rank increase:
        if rKi > rK
            K = Ki; rK = rKi
        end

        remove = findall(x->x==i,M);
        deleteat!(M,remove)
    end
    
    # initialize Rational Matrices:
    AK = (A[K,:]); bK = reshape(b[K],:,1);
    SA = matrix_space(QQ, size(AK)[1], size(AK)[2])
    Sb = matrix_space(QQ, size(bK)[1], size(bK)[2])
    
    # convert to AbstractAlgebra objects:
    QAK = SA(AK); QbK = Sb(bK); #println("$QAK \n $QbK")
    R = vcat(ones(Rational{Int64},1,1),Matrix{Rational{Int64}}(solve(QAK, QbK, side = :right)));

    return K,R

end



using DelimitedFiles

function ddm(A,K,R, init = false, save = true, filename = "filename")
    
    # if no initial DD pair:
    if init == false
        # do not include homogenizing coordinate:
        K, R = initial_DD(A);
    end
    
    # input size:
    m = size(A)[1]; d = size(A)[2];
    
    # remaining inequalities:
    M = [i for i in 1:m if i ∉ K];
    
    # create matrices to iterate over:
    Riter = copy(R); Aiter = copy(A[K,:]);
    
    for i in M
        # define facet inserted at step i:
        hi = A[i,:];

        # initialize vertices generated at step i:
        Ri = Matrix{Rational{Int64}}(undef,d,0);

        # segment Rn into H+, H0, H-:
        rp, rt, rn = regions(hi,Riter);
        
        # For all vertices in H-:
        for j in rn

            # current vertex:
            Rj = Riter[:,j];

            # inner product:
            a = Rational{Int64}(transpose(hi)*Rj);

            # all neighbors of j such that j in rp:
            Njp = adjacency(j,rp,Aiter,Riter);


            # new vertices:
            for k in Njp

                # define positive vertex:
                Rjp = Riter[:,k];

                # inner product:
                b = Rational{Int64}(transpose(hi)*Rjp);

                if (b-a) != 0

                    p = b//(b-a);


                    if (p > 0)

                        # generate new point:
                        Rjk = p*Rj+(1-p)*Rjp;

                        # add to representation matrix
                        Ri = hcat(Ri,Rjk);
                        
                    end
                end
            end
        end

        println(Ri)

        
        # save new vertices with every iteration:
        if save == true
            open(filename, "w") do io
                writedlm(io, Riter, ',')
            end
        end
                                

        # update inequalities:
        Aiter = vcat(Aiter,transpose(hi));

        # Remove violating extreme rays:
        Keep = [i for i in 1:size(Riter)[2] if i ∉ rn];
        Riter = Riter[:,Keep];

        # Add newly generated rays:
        Riter = hcat(Riter,Ri);
        
        # online updates:
        println("Inequality: ", i, "\n")
        println("Intermediate DD vertices: ", size(Riter),"\n");

    end
    return Riter
end

####################################################################################
####################################################################################
####################################################################################

# Rank functions:

# input: x (vector), inequalities (array)
# output: int array of index values
function tight_inequalities(x,inequalities)
    x = Vector{Rational{Int64}}(x);
    A = Matrix{Rational{Int64}}(inequalities);
    b = A*x;
    return findall(y->y==0,b)
end
""";;