############################################################################
############################################################################
##########################      TwoDim      ################################
############################################################################
############################################################################

"""
Background on simplicial distributions:

The measurement scenarios of interest can be encoded into 2-dimensional spaces.
In the simplicial framework these are collections of triangles glued along
edges and vertices.

Pauli measurements will be associated with an edge in this space. Pauli
operators in a stabilizer group make up the edges of a triangle; for instance,
{XI,IX,XX} can be three edges of a triangle. There are four possible stabilizer
groups containing these Pauli operators (up to signs). Each suitable choice of sign
corresponds to a stabilizer state which, given a Hermitian operator A, defines an
inequality through Tr(APi_S) >= 0.
"""



function nerve_face_map(s,d)
    """
    Given an n-tuple of dit-strings, calculate the faces, i.e. the 
    (n-1)-dit-strings associated to the n-1 faces of the simplex.
    Given a joint outcome, or tuple (a,b), this assigns an outcome to
    each of the Pauli measurements that respects the fact that the Paulis
    either multiply to +I or -I.
    
    Parameters:
        - s: Array{Int}: dit-string
        - d: Int: Dimensionality
        
    Returns:
        - Array{Array{Int}}: Array of faces of n-tuple of dit-strings
    
    """
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
        push!(ds,di_s);
    end
    return ds
end

######################################

function edge_outcome_assignments(ab,T=[])
    """
    Compute the outcome assignments for a two-bit string based on given twisting parameters.
    E.g., for Paulis A,B,AB, the tuple (a,b) assigns a to A, b to B, and a+b mod 2 to AB
    if the operators multiply to I and a+b+1 mod 2 to AB if they multiply to -I.
    
    Parameters:
        - ab: Array{Int}: Two-bit string (array of 0/1)
        - T: Array{Int} (Optional): Twisting parameters [beta, di]. Default is empty.
    
    Returns:
        - Array{Int}: Outcome assignments for [b, a+b+beta, a]
    """
    # default input no twisting, otherwise T=[beta,di] determines twisting:
    if isempty(T) == false; beta = T[1]; di = T[2]; else;    beta = 0; di = 1; end;
    
    # edge assignments with no twisting: (simplicial distributions on 2-simplex)
    edge_assignments = [ab[2],((ab[1]+ab[2]) % 2),ab[1]];

    # edge assignments with twisting:
    edge_assignments[di] = ((edge_assignments[di]+beta) % 2);

    return edge_assignments
end

######################################

function spaces_to_inequalities_EDGE(X1,X2,T=[])
    """
    Convert spaces to inequalities based on edges.
    
    Parameters:
        - X1: Array{Int}: Edges
        - X2: Array{Tuple}: Triangles
        - T: Array{Array{Int}} (Optional): Twisting parameters [beta, di]. Default is empty.
    
    Returns:
        - Array{Int}: Matrix of inequalities
    """
    # number of edges:
    N = length(X1); d = 2;
    # edges:
    E = [X1[i][1] for i in 1:N]; dict = Dict(zip(E,[i for i in 1:N]));

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

######################################

using Polymake
const pm = Polymake

function two_dimensional_distributions(X,COORDINATES,d = 2,T=[])
    """
    Compute two-dimensional distributions.
    
    Parameters:
        - X: Tuple: (X1,X2)
        - COORDINATES: Array{Array}: Coordinates
        - d: Int: Dimensionality
        - T: Array{Array{Int}} (Optional): Twisting parameters [beta, di]. Default is empty.
    
    Returns:
        - pm.polytope.Polytope: Polytope object representing two-dimensional distributions
    """

    # input to TwoDim is a simplicial set in the convention of SimpSet. 
    X1 = X[2]; X2 = X[3];

    M = spaces_to_inequalities_EDGE(X1,X2,T);

    return pm.polytope.Polytope(INEQUALITIES=M);
end