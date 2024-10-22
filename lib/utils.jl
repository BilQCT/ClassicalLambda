# use the Oscar package:
using Polymake
const pm=Polymake

using LinearAlgebra

using Nemo
const nm = Nemo;

function nemo_to_julia(nM)
    k = size(nM)[1]; d = size(nM)[2];
    M = Matrix{Rational{Int64}}(undef,0,d);
    for i in 1:k
        row = Vector{Rational{Int64}}([]);
        for j in 1:d
            r = convert(Rational{Int64},nM[i,j]);
            row = vcat(row,[r]);
        end
        M = vcat(M,transpose(row));
    end
    return M
end

######################################

# compute rank of matrix of tight inequalities
# input: Vector{Int}: array of integer indices corresponding to tight inequalities
# input: Matrix: inequalities
function rank_of_point(tight_array,inequalities)
    if typeof(tight_array) == String
        return "Not a vertex"
    else
        return rank(Array(inequalities[tight_array,:]))
    end
end

######################################

# input: Vector: inequality (dx1)
# input: Matrix: points (n x d)
# output: vector of indices such that R[i,:]*h = 0.
function tight_at_inequality(h,R)
    R = Array(R); h = Vector(h);
    Rh = R*h;
    return findall(x->x == 0,Rh)
end

######################################

# input: vector: p1, p2
# input: matrix: A
function joint_rank(p1,p2,A)
    idx1 = tight_inequalities(p1,A);
    idx2 = tight_inequalities(p2,A);
    intersection = collect(intersect(Set(idx1),Set(idx2)));
    AZ = Array(A[intersection,:]);
    return rank(AZ)
end

# input: vertex (d x 1), inequalities (m x d)
function violated_inequality(v,A)
    v = Vector(v); A = Array(A);
    Av = A*v;
    return findall(x->x < 0,Av)
end

######################################

# input: vector in R(d+1); matrix of dimension m x (d+1)
function vertices_violating_facet(h,R)
    h = Vector(h); R = Array(R);
    Rh = R*h;
    v_idx = findall(x->x<0,Rh); V = R[v_idx,:];
    return v_idx, V
end

######################################

# input: Vector: h: facet (1xd),
# input: Vector: v_minus: violating vector (d x 1)
# input: Vector: v_plus: satisfying vector (d x 1)
# output: Vector: new vertex (d x 1)
function new_vertex(v_minus,v_plus,h)
    d = size(v_minus)[1];
    if size(h)[1] == d
        a = Rational{Int64}((transpose(h)*v_minus)[1]);
        b = Rational{Int64}((transpose(h)*v_plus)[1]);
    else
        a = Rational{Int64}(((h)*v_minus)[1]);
        b = Rational{Int64}(((h)*v_plus)[1]);
    end
    if (b-a) == 0
        error("No vertex possible 1")
    end
    p = Rational{Int64}(b)//Rational{Int64}(b-a);
    if p < 0
        error("No vertex possible 2")
    elseif p == 0
        error("No new vertex possible 3")
    else
        return p*v_minus+(1-p)*v_plus
    end
end


####################################################################################
####################################################################################
####################################################################################

# Lifting functions


# input: x in R9; A in R^60x16; I - subset of inequalities; s - subset of coordinates
# outout: pm polytope object: P subset R7.
function lift_mp1_local(v,A,S,s)
    # convert to Julia Rational Array
    v = Vector{Rational{Int64}}(v);
    A = Matrix{Rational{Int64}}(A);
    # index for all coordinates:
    x = ([i for i in 1:size(A)[2]]);
    # coordinates indices not in x_set
    x_perp = sort(collect(setdiff(Set(x),Set(s))));
    # subset of inequalities
    A_loc = A[S,:]; Bx = A_loc[:,s]*v;
    # Redefine constant column
    A_loc[:,1] .= A_loc[:,1] + Bx;
    return pm.polytope.Polytope(INEQUALITIES = A_loc[:,x_perp])
end

######################################

# input: nonlocal point x, polytope P with vertices v
# output: polytope tilde_P as convex hull of vertices: tilde_v = (v,x).
function append_projection(x,P)
    R = (P.VERTICES); m = size(R)[1]; d = size(R)[2]; R_prime = transpose(x);
    for i in 1:(m-1)
        R_prime = vcat(R_prime,transpose(x));
    end
    return hcat(R,R_prime)
end

######################################

# input: point x to be lifted
# output: Polytope of feasible lifts
function lift_mp1(v,A,S,s)
    Px = lift_mp1_local(v,A,S,s);
    V = append_projection(v,Px);
    return pm.polytope.Polytope(POINTS = V);
end


####################################################################################
####################################################################################
####################################################################################


# input: julia (numeric) vector (perhaps could be array, didn't check):
# output: array of entries with zero probability
function find_zeros(x)
    Z = findall(y->y==0, x);
    return Z;
end

######################################

# inequality matrix:
#INPUT: num. of prob. parameters:
#OUTPUT: matrix of inequalities for polymake:
function inequality_matrix(N)
    id = Matrix{Int64}(I,N,N); #identity
    zero = zeros(Int64,N);         #column of zeros
    ineq = hcat(zero,id);              #concatenate
    return ineq
end

######################################

# function: take vector x and generating matrix A and representation matrix R
#           and compute possible convex decompositions of x in terms of R
# input: x: vector ((d,)), Matrix (mxd), Matrix (nxd)
# input: V: return the set of vertices used in the decomposition. Default set to false.
# output: pm polytope
function convex_decomposition(x,A,R,V = false)
    # convert (possible) pm matrices into julia matrices:
    A = Matrix{Rational{Int64}}(A); R = transpose(Matrix{Rational{Int64}}(R));
    B = transpose(A*R); x = Vector{Rational{Int64}}(x);
    
    # Find tight inequalities:
    px = A*x; Zx = find_zeros(px);
    
    # find all columns of R with common zeros:
    V_idx = [];
    for i in 1:size(R)[2]; y = R[:,i]; Zy = find_zeros(A*y);
        if (intersect(Set(Zx),Set(Zy))) == Set(Zx);
            push!(V_idx,i);
        end;
    end
    
    # Common tight inequalities reduce the number of vertices that may
    # enter the convex decomposition:
    if isempty(V_idx) == true; V = R;
    else; V = (R[:,V_idx]); end;
    
    # create feasible polytope: all convex decompositions of x in R[V_idx]:
    # i.e., x = Ry: y>=0, sum y = 1:
    
    # nonnegativity conditions:
    I = inequality_matrix(size(V)[2]);
    # combine into single array
    E = hcat(-x,V)[2:end,:];
    # create feasible polytope:
    Q = pm.polytope.Polytope(INEQUALITIES=I, EQUATIONS = E);
    if V == false; return Q;
    else; return V,Q; end
end


####################################################################################
####################################################################################
####################################################################################


# Data-type conversions:


# change polymake to Julia array of arrays, which is easier for searching:
function polymake_to_array(P,T)
    A = [];
    for i in 1:size(P)[1]
        r = P[i,:];
        r_array = Vector{T}([k for k in r]);
        push!(A,r_array);
    end
    return A
end

######################################

function array_to_matrix(A,T)
    m = length(A); d = length(A[1]);
    M = Matrix{T}(undef,0,d);
    for i in 1:m
        r = transpose(A[i]);
        M = vcat(M,r);
    end
    return M
end


######################################


# recursive loop to construct all dit strings:
function generate_dit_strings(d,N,n,s)
    if n < N
        Zd = [(i-1) for i in 1:d] ; dit_strings = [];
        for i in 1:length(s)
            for x in Zd
                #println(i)
                push!(dit_strings,push!(copy(s[i]),x));
            end
        end
        s = dit_strings; n = n+1;
        generate_dit_strings(d,N,n,s)
    else
        return s
    end
end

######################################

# input: outcomes (d) and number of generators (N) (N=dim(simplex))
function all_dit_strings(d,N)
    s = [[i-1] for i in 1:d]; n = 1;
    return generate_dit_strings(d,N,n,s)
end
