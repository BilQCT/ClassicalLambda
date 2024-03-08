# Edge coordinates:

# input: dit-string
# output: array of faces of n-tuple of dit-strings
# using Eq.(12) of SimpCont
function nerve_face_map(s,d)
    ds = [];
    for i in 1:(length(s)+1)
        di_s = [];
        # d0 face:
        if i == 1
            di_s = s[2:end]
        # d1 face:
        elseif i==2
            a = (s[1]+s[2]) % d;
            di_s = vcat([a],[s[j] for j in (i+1):length(s)]);
        # di face: 2 < i < n+1
        elseif (i>2) && (i < (length(s)+1))
            a = [s[j] for j in 1:(i-2)];
            b = [(s[i-1]+s[i])%d];
            if length(vcat(a,b)) == length(s)-1
                di_s = vcat(a,b)
            else
                c = [s[j] for j in (i+1):length(s)];
                di_s = vcat(vcat(a,b),c);
            end
        # dn face:
        else
            di_s = s[1:(end-1)]
        end
        #println(di_s)
        push!(ds,di_s);
    end
    return ds
end




# input: two-bit string (ab) (array of 0/1),
# input: twisting (optional): t = [beta,di]: Twisting and the twisted edge: beta in {0,1}, {d0,d1,d2} = {1,2,3} in julia
# input: By default, T is empty and there is no twisting.
#output: array of length three giving outcome assignment for [b,a+b+beta,a]:

function edge_outcome_assignments(ab,T=[])
    # default input no twisting, otherwise T=[beta,di] determines twisting:
    if isempty(T) == false; beta = T[1]; di = T[2]; else;    beta = 0; di = 1; end;
    
    # edge assignments with no twisting: (simplicial distributions on 2-simplex)
    edge_assignments = [ab[2],((ab[1]+ab[2]) % 2),ab[1]];

    # edge assignments with twisting:
    edge_assignments[di] = ((edge_assignments[di]+beta) % 2);

    return edge_assignments
end






function spaces_to_inequalities_EDGE(X1,X2,T=[])
    # number of edges:
    N = length(X1); d = 2;
    # edges
    E = [X1[i][1] for i in 1:N]; dict = Dict(zip(E,[i for i in 1:N]));
    
    #println(dict)

    # extract twisted edges:
    if isempty(T) == false; twisted_idx = [T[i][1] for i in 1:length(T)]; else; twisted_idx = []; end;
    
    # bit strings
    Z2_2 = [[i,j] for i in [0,1] for j in [0,1]];

    # initialize A matrix for P(A,b)
    A = zeros(Int64, 4*length(X2), N+1);

    # For each 2-simplex
    for i in 1:length(X2)
        # extract beta value and face to be applied to:
        if X2[i][1] in twisted_idx
            # find unique identifier: twist = [beta,di]:
            idx = (findall(x->x==X2[i][1],twisted_idx))[1]; twist = [T[idx][2],T[idx][3]];
        else
            twist = [];
        end
        
        
        # assignment of outcomes to edges: (ab) -> [b,a+b,a]:
        edge_assignments = [edge_outcome_assignments(ab,twist) for ab in Z2_2];
        # assign coefficients of edges:
        edge_coefficients = [[(-1)^s[1],(-1)^s[2],(-1)^s[3]] for s in edge_assignments];
        
        
        for j in 1:4
            # initialize row vector of zeros
            a = zeros(Int8,N+1);
            
            
            # row vector defined by sum of projectors (-1)^a*P1+(-1)^b*P2+(-1)^a+b*P3
            for k in 1:3
                P = zeros(Int64,N+1);
                e = X2[i][2][k];
                # if degenerate modify constant term:
                if e < 0
                    P[1] = P[1]+edge_coefficients[j][k];
                elseif e > 0
                    x = dict[e];
                    P[x+1] = edge_coefficients[j][k];
                end
                a = a + P;
            end
            
            # add constant term:
            a[1] = a[1]+1;
            
            # determine ordering or inequalities:
            row = 4*(i-1)+j;
            # append to matrix
            A[row,:] = a;
            
        end
    end
    
    return A
end




################################################################################
################################################################################
############################ Probability Coordinates ###########################
################################################################################
################################################################################




# Determine all faces (edges) ocurring in each higher simplex (triangle)
# input: two-dimensional simplicial set (applies to Delta_n, Delta_n-1, as well)
# output: array of [face,simplex] for each f in Delta_n-1
function faces_in_simplex(X1,X2)
    F = [];
    for x1 in X1
        # isolate edge identifier:
        tau = x1[1]; f = [];

        # check if tau is a face of each triangle sigma:
        for x2 in X2
            sigma = x2[1]; faces = x2[2]; tau_in_sigma = findall(y->y==tau,faces);

            # if tau is a face, register it:
            if isempty(tau_in_sigma) == false
                for k in tau_in_sigma; push!(f,[k,sigma]) end;
            end
        end
        push!(F,f);
    end
    return F
end




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




# input: outcomes (d) and number of generators (N) (N=dim(simplex))
function all_dit_strings(d,N)
    s = [[i-1] for i in 1:d]; n = 1;
    return generate_dit_strings(d,N,n,s)
end




# function: find all dit-strings of length k+1 with face maps coinciding on target dit-string s_k
# input: outcomes (d) (Z_d); target dit-string (s_k); di (ith face map)
function compatible_face_maps(d,s_k,i)
    k = length(s_k); higher_bitstrings = all_dit_strings(d,k+1);

    # simplices sn of (NZd)_n such that di(sn) = sn_1
    S = [];
    for s in higher_bitstrings
        ds = nerve_face_map(s,d);
        if ds[i] == s_k
            push!(S,s)
        end
    end
    return S
end




function indicator_function(bit_string,faces)
    if bit_string in faces; return 1; else; return 0; end;
end


################################################################################
################################################################################


# call Linear Algebra package:
using LinearAlgebra

# Written with Z2 and two-dimensional spaces in mind, but could be generalized.
# input: Array of nondegenerate simplices whose dimension is given by length of array:
function nonnegative(Xn,d)

    # extract dimension of simplices:
    dim = (length(Xn[1][2])-1);

    # Determine number of probability parameters to marginalize over:
    L = length(Xn); N_prob = d^(dim); N = N_prob*L;

    # create matrix of inequalities:
    Id = Matrix{Int}(I,N,N); Z = (zeros(Int8,N));
    #concatenate
    A = hcat(Z,Id);             
    return A
end






# function: normalization constraints: Outcomes in Zd:
# input: Two-dimensional simplicial set. Can be generalized to higher dimensions.
# output: matrix of norm. constraints.
function normalization_matrix(Xn,d)
    L = length(Xn); simplex_dim = [(length(Xn[i][2])-1) for i in 1:L];
    N_simplex = [d^x for x in simplex_dim]; N = sum(N_simplex);

    # initialize matrix of constraints:
    init = zeros(Int8, L, N); col = -ones(Int8,L);
    NORM = hcat(col,init);

    # set counter:
    m = 0; for i in range(1,L) n=m+N_simplex[i]; NORM[i,(m+2):n+1] .= 1; m=n;  end
    return NORM
end



function collapse_conditions(d,X2,T = [])
    # Determine number of probability parameters to marginalize over:
    L = length(X2); N_prob = d^(2); N = N_prob*L;
    
    # Triangle indices:
    S = [X2[i][1] for i in 1:L]
    
    # Twisted triangles:
    Ti = [T[i][1] for i in 1:length(T)];
    twist_dict = Dict(zip(Ti,[i for i in 1:length(T)]))
    
    # initialize constraint matrix
    C = Array{Int64}(undef, 0, N+1);
    
    # generate all dit-strings of length 2:
    Zd_2 = all_dit_strings(d,2);
    # dictionary: map z in Zd_2 to {1,..,d^2}:
    dict = Dict(zip(Zd_2,[i for i in 1:(d^2)]))
    
    # search for degenerate edges:
    for i in 1:L
        E = X2[i][2]; s = X2[i][1];
        D = findall(x->x < 0,E);
        
        # Twisted faces:
        if s in Ti
            beta = T[twist_dict[s]][2];
            dk = T[twist_dict[s]][3]
        else
            dk = []; beta = 0;
        end
        
        
        # elements of Deg are faces of 2-simplex
        for dd in D
            if dd == dk; z = beta; else; z = 0; end;
            f = compatible_face_maps(d,[z],dd);
            f_perp = [z for z in Zd_2 if z ∉ f];
            for k in f_perp
                n = dict[k]; column = (d^2)*(i-1)+1+n;
                
                # create row
                c = zeros(Rational{Int64},(1,N+1));
                c[column] = 1; C = vcat(C,c)
            end
        end
    end
    return C
end




using Combinatorics

# input: two-dimensional simplicial set: (X1, X2).
# input: twisting: T = [Ti], where Ti = [idx,beta,dk]
# output: matrix of compatibility relationships
function two_dimensional_compatibility(X1, X2, d, T=[])

    # extract dimension of simplices:
    dim1 = (length(X1[1][2])-1); dim2 = (length(X2[1][2])-1);

    # Determine number of probability parameters to marginalize over:
    L = length(X2); N_prob = d^(dim2); N = N_prob*L;
    
    # Array of edges occuring in each triangle: [f_i,sigma_k]
    F = faces_in_simplex(X1,X2);

    #println(F)

    # map simplex identifiers to position index in array (X2):
    dict = Dict(zip([x[1] for x in X2],[i for i in 1:length(X2)]));

    # NS condition whenever f in X1 occurs twice:
    # for all edges that appear more than once, we impose compatibility:
    NS = Array{Int64}(undef, 0, N+1); # initialize NS matrix
    for i in 1:length(F)
        # n-simplices (f[2]) for which the (n-1)-simplex (edge) (f[1]) is a face: 
        if length(F[i]) > 1
            # generate all pairs of simplices glued along same edge:
            F_combinations = combinations(F[i],2);
            for f in F_combinations
                # initialize row vector for constraint:
                init1 = zeros(Int64, 1, N+1); init2 = zeros(Int64, 1, N+1);
                
                # specify all n (n-1) simplices of NZ_d:
                lower_outcomes = all_dit_strings(d,dim1); higher_outcomes = all_dit_strings(d,dim2);

                # extract simplex idx:
                simplex_1 = f[1][2]; simplex_2= f[2][2];
                # extract face maps:
                di_1 = f[1][1]; di_2 = f[2][1];

                # twisted edges:
                idx1 = findall(x->x[1]==simplex_1,T); idx2 = findall(x->x[1]==simplex_2,T);

                if isempty(idx1) == false
                    beta1 = T[idx1[1]][2]; face1 = T[idx1[1]][3];
                else
                    beta1 = 0; face1 = copy(di_1);
                end

                if isempty(idx2) == false
                    beta2 = T[idx2[1]][2]; face2 = T[idx2[1]][3];
                else
                    beta2 = 0; face2 = copy(di_2);
                end


                for s in lower_outcomes
                    #println(s)
                    # Determine parameters that enter into marginalization:
                    if (isempty(idx1)==false) && (face1 == di_1)
                        s_prime = [(a + beta1) % 2 for a in s];
                        s1 = compatible_face_maps(d,s_prime,di_1);
                    else
                        s1 = compatible_face_maps(d,s,di_1);
                    end
                    
                    if (isempty(idx2)==false) && (face2 == di_2)
                        s_prime = [(a + beta2) % 2 for a in s];
                        s2 = compatible_face_maps(d,s_prime,di_2);
                    else
                        s2 = compatible_face_maps(d,s,di_2);
                    end

                    eta1 = [indicator_function(x,s1) for x in higher_outcomes];
                    eta2 = [-indicator_function(x,s2) for x in higher_outcomes];

                    # define map of d^n length array into N x d^n length array:
                    sigma1 = f[1][2]; i1 = dict[sigma1]; sigma2 = f[2][2]; i2 = dict[sigma2];
                    END1 = i1*(N_prob)+1; INT1 = END1-(N_prob)+1;
                    END2 = i2*(N_prob)+1; INT2 = END2-(N_prob)+1;#

                    # define interval where constraint is non-trivial:
                    p1 = [i for i in INT1:END1]; p2 = [i for i in INT2:END2];
                    init1[p1] = eta1; init2[p2] = eta2; init = init1+init2;
                    NS = vcat(NS,init);
                end
            end
        end
    end
    return NS
end



function spaces_to_inequalities_PROB(X1,X2,d,T=[])
    # define matrix of nonnegativity inequalities:
    A = nonnegative(X2,d);
    N = normalization_matrix(X2,d);
    C = two_dimensional_compatibility(X1,X2,d,T);
    Co = collapse_conditions(d,X2,T)
    return A, vcat(N,C,Co)
end




################################################################################
################################################################################
################################ Main function #################################
################################################################################
################################################################################


using Polymake
const pm = Polymake
# input: X = (X1,X2)
function two_dimensional_distributions(X,COORDINATES,d = 2,T=[])
    # input to TwoDim is a simplicial set in the convention of SimpSet. 
    #X = collect(X);
    # extract edges and triangles:
    X1 = X[2]; X2 = X[3];

    if (COORDINATES == "EDGE")
        if (d == 2)
            M = spaces_to_inequalities_EDGE(X1,X2,T);
            return pm.polytope.Polytope(INEQUALITIES=M);
        else
            return "Edge coordinates only valid for d = 2."
        end

    elseif COORDINATES == "PROB"
        A, E = spaces_to_inequalities_PROB(X1,X2,d,T);
        return pm.polytope.Polytope(INEQUALITIES=A, EQUATIONS = E);
        
    else
        return "Not a valid input";
    end
end