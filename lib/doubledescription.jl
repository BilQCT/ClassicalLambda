# use the Oscar package:
using Polymake
const pm=Polymake
using LinearAlgebra

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

using LinearAlgebra

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

# input: Vector: inequality (dx1)
# input: Matrix: points (n x d)
# output: vector of indices such that R[i,:]*h = 0.
function tight_at_inequality(h,R)
    R = Array(R); h = Vector(h);
    Rh = R*h;
    return findall(x->x == 0,Rh)
end

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

# input: vector in R(d+1); matrix of dimension m x (d+1)
function vertices_violating_facet(h,R)
    h = Vector(h); R = Array(R);
    Rh = R*h;
    v_idx = findall(x->x<0,Rh); V = R[v_idx,:];
    return v_idx, V
end

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


function array_to_matrix(A,T)
    m = length(A); d = length(A[1]);
    M = Matrix{T}(undef,0,d);
    for i in 1:m
        r = transpose(A[i]);
        M = vcat(M,r);
    end
    return M
end

# matching orbits of P (current polytope) with orbits of Q (target polytope):
function matching_orbits(P, Q, order = [])
    # determine representative vertices of P in terms of orbits of Q:
    P_vertices = polymake_to_array(representative_vertices(P));
    # all vertices of Q:
    Q_vertices = polymake_to_array(Q.VERTICES);
    # sort all vertices of Q into orbits:
    if isempty(order) == true
        Q_orbits = vertex_orbit(Q);
    else
        Q_orbits = vertex_orbit(Q)[order];
    end
    # put matching vertices into array:
    match_array = []; vertices_array = [[] for i in 1:length(Q_orbits)];
    for i in 1:size(P_vertices)[1]
        v = P_vertices[i]; matches = findall(y->y==v,Q_vertices);
        if isempty(matches) == false
            match_idx = matches[1];
            for j in 1:length(Q_orbits)
                type = findall(y->y==match_idx,Q_orbits[j]);
                if length(type) > 0; push!(match_array,j), push!(vertices_array[j],i)  end;
            end
        end
    end
    return [sort(collect(Set(match_array))),vertices_array]
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

# input: nonlocal point x, polytope P with vertices v
# output: polytope tilde_P as convex hull of vertices: tilde_v = (v,x).
function append_projection(x,P)
    R = (P.VERTICES); m = size(R)[1]; d = size(R)[2]; R_prime = transpose(x);
    for i in 1:(m-1)
        R_prime = vcat(R_prime,transpose(x));
    end
    return hcat(R,R_prime)
end

# input: point x to be lifted
# output: Polytope of feasible lifts
function lift_mp1(v,A,S,s)
    Px = lift_mp1_local(v,A,S,s);
    V = append_projection(v,Px);
    return pm.polytope.Polytope(POINTS = V);
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




####################################################################################
####################################################################################
####################################################################################

# Latex functions

function vec_to_str(v)
    strg_array = [];
    for i in 1:length(v)
        strg = string(v[i])
        if length(strg) == 5
            strg = strg[1:2]*strg[(end-1):end]
        else
            strg = strg[1]*strg[(end-1):end]
        end
        push!(strg_array,strg);
    end
    return strg_array
end

function classical_string(idx,typ)
    stringg = "\\noindent \n 
    \\textbf{Classical vertex of \$\\mathbf{\\Lambda_2}\$: \\# \$\\mathbf{"*string(idx)*"}\$, Type \$\\mathbf{"*string(typ)*"}\$}"
    return stringg
end

function nonclassical_string(idx,typ)
    stringg = "\\noindent \n 
    \\textbf{Non-classical vertex of \$\\mathbf{\\Lambda_2}\$: \\# \$\\mathbf{"*string(idx)*"}\$, Type \$\\mathbf{"*string(typ)*"}\$}"
    return stringg
end

# Objective: create tex code for a table of Pauli coefficients for vertex v:
# INPUT: List of Pauli coefficients:
# OUTPUT: readable stringing/latex code:

function h_tab(v)
    stringg = "\$\$\\begin{array}{cccccccccccccccc} \n
    II & XI & YI & ZI & IX & IY & IZ & XX & XY & XZ & YX & YY & YZ & ZX & ZY & ZZ  \\\\ \n
    "*v[1]*" & "*v[2]*" & "*v[3]*" & "*v[4]*" & "*v[5]*" & "*v[6]*" & \n
    "*v[7]*" & "*v[8]*" & "*v[9]*" & "*v[10]*" & "*v[11]*" &
    "*v[12]*" & "*v[13]*" & "*v[14]*" & "*v[15]*" & "*v[16]*"\n
    "*"\\end{array}\$\$\n"
    return stringg

end


# OBJ: create stringing of probability table for tex
# INPUT: list of probabilities:
# OUTPUT: stringing for latex. 

function pl_tab(p)
    stringg = "\\begin{adjustbox}{width=\\columnwidth,center}
    \\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|l|l|l|l|l|l|}
    \\hline
    s \$\\backslash\$ I   & XI,IX & XI,IY & XI,IZ & YI,IX & YI,IY & YI,IZ & ZI,IX & ZI,IY & ZI,IZ & XY,YX & YZ,ZY & ZX,XZ & XY,ZX & YX,XZ & ZZ,YY \\\\ \\hline \n
    00 & "*p[1]*" & "*p[5]*" & "*p[9]*" & "*p[13]*" & "*p[17]*" & "*p[21]*" & "*p[25]*" & "*p[29]*" & "*p[33]*" & "*p[37]*" & "*p[41]*" & "*p[45]*" & "*p[49]*" & "*p[53]*" & "*p[57] *"\\\\ \\hline \n\
    01 & "*p[2]*" & "*p[6]*" & "*p[10]*" & "*p[14]*" & "*p[18]*" & "*p[22]*" & "*p[26]*" & "*p[30]*" & "*p[34]*" & "*p[38]*" & "*p[42]*" & "*p[46]*" & "*p[50]*" & "*p[54]*" & "*p[58] *"\\\\ \\hline    \n\
    10 & "*p[3]*" & "*p[7]*" & "*p[11]*" & "*p[15]*" & "*p[19]*" & "*p[23]*" & "*p[27]*" & "*p[31]*" & "*p[35]*" & "*p[39]*" & "*p[43]*" & "*p[47]*" & "*p[51]*" & "*p[55]*" & "*p[59] *"\\\\ \\hline     \n\
    11 & "*p[4]*" & "*p[8]*" & "*p[12]*" & "*p[16]*" & "*p[20]*" & "*p[24]*" & "*p[28]*" & "*p[32]*" & "*p[36]*" & "*p[40]*" & "*p[44]*" & "*p[48]*" & "*p[52]*" & "*p[56]*" & "*p[60] *"\\\\ \\hline    \n\
    \\end{tabular}
    \\end{adjustbox}"
    return stringg
end

# nonlocal vertices:
# Objective: create tex code for a table of Pauli coefficients for vertex v:
# INPUT: List of Pauli coefficients:
# OUTPUT: readable stringing/latex code:

function nl_h_tab(v)
    stringg = "\$\$\\begin{array}{cccccccccccccccc} \n
    II & XX & XY & XZ & YX & YY & YZ & ZX & ZY & ZZ  \\\\ \n
    "*v[1]*" & "*v[2]*" & "*v[3]*" & "*v[4]*" & "*v[5]*" & "*v[6]*" & \n
    "*v[7]*" & "*v[8]*" & "*v[9]*" & "*v[10]*"\\end{array}\$\$\n"
    return stringg

end


# OBJ: create stringing of probability table for tex
# INPUT: list of probabilities:
# OUTPUT: stringing for latex. 

function nl_p_tab(p)
    stringg = "\\begin{adjustbox}{width=0.85\\columnwidth,center}
    \\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|l|l|l|l|l|l|}
    \\hline
    s \$\\backslash\$ I & XY,YX & YZ,ZY & ZX,XZ & XY,ZX & YX,XZ & ZZ,YY \\\\ \\hline \n
    00 & "*p[1]*" & "*p[5]*" & "*p[9]*" & "*p[13]*" & "*p[17]*" & "*p[21]*"\\\\ \\hline \n\
    01 & "*p[2]*" & "*p[6]*" & "*p[10]*" & "*p[14]*" & "*p[18]*" & "*p[22]*"\\\\ \\hline    \n\
    10 & "*p[3]*" & "*p[7]*" & "*p[11]*" & "*p[15]*" & "*p[19]*" & "*p[23]*"\\\\ \\hline     \n\
    11 & "*p[4]*" & "*p[8]*" & "*p[12]*" & "*p[16]*" & "*p[20]*" & "*p[24]*"\\\\ \\hline    \n\
    \\end{tabular}
    \\end{adjustbox}"
    return stringg
end

function f_string(v)
    stringg = "\n\
    \\noindent
    \$\\mathbf{\\# "*string(v)*"}\$: \n
    \\vskip .25 cm \n
    "
    return stringg
end

function type_string(v,typ)
    stringg = "\n
    \\noindent
    \$\\mathbf{\\# "*string(v)*" :}\$ \\textbf{Subtype "*string(typ)*"}: \n
    \\vskip .25 cm \n
    "
    return stringg
end

function skip_string(len)
    stringg = """\n\
    \\vskip """*string(len)*""" cm \n\
    """
    return stringg
end

function section(v)
    stringg = """\n\
    \\noindent \n\
    \\section{\\textbf{Neighbors of Canonical Vertex: Type """*string(v)*""":}}\\label{"""*string(v*1)*"""}\n\
    """
    return stringg
end

function section_sub(v,s)
    stringg = """\n\
    \\noindent \n\
    \\section{\\textbf{Neighbors of Canonical Vertex: Type """*string(v)*""": Subtype """*string(s)*\
    """}}\\label{"""*string(v*1)*"""}\n\
    """
    return stringg
end

function sec_stringg(sec_string, lab_string)
    stringg = """\n\
    \\noindent \n\
    \\section{\\textbf{"""*string(sec_string)*""":}}\\label{"""*string(lab_string)*"""}\n\
    """
    return stringg
end

function latex_matrix(A)
    m = size(A)[1]; n = size(A)[2]; s = "";
    for i in 1:m
        for j in 1:n
            if j < n; s = s*str_func(A[i,j])*" & ";
            else; s = s*str_func(A[i,j]); end;
        end
        if i < m; s = s*" \\\\ \n ";
        else; s = s*" \n "; end
    end
    return print(s)
end
