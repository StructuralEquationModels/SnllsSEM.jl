############################################################################
### Types
############################################################################

struct SNLLS{N1, N2, N3, I1, I2, I3, I6, M1, M2, M3, M4, N4, I4, M5, I5, D} <: SemImply

    q_directed::N1
    q_undirected::N2
    size_σ::N3

    A_indices_linear::I1
    A_indices_cartesian::I2
    S_indices::I3
    σ_indices::I6

    A_pre::M1
    I_A::M2
    G::M3
    ∇G::M4

    q_mean::N4
    M_indices::I4
    G_μ::M5
    G_μ_indices::I5

    identifier::D

end

############################################################################
### Constructors
############################################################################

function SNLLS(;
    specification,
    start_val = start_fabin3,
    kwargs...)

    specification = spec_mean

    ram_matrices = RAMMatrices(specification)
    identifier = StructuralEquationModels.identifier(ram_matrices)
    
    A_ind, S_ind, F_ind, M_ind, parameters = ram_matrices.A_ind, 
        ram_matrices.S_ind, ram_matrices.F_ind, ram_matrices.M_ind, ram_matrices.parameters
    
    n_var, n_nod = ram_matrices.size_F
    
    σ_indices = findall(isone, LowerTriangular(ones(n_var, n_var)))
    
    S = zeros(n_nod, n_nod)
    S = Matrix{Any}(S)
    for (i, par) in enumerate(parameters)
        for ind in S_ind[i]
            S[ind] = par
        end
    end
    
    S_indices = get_parameter_indices(parameters, S; index_function = eachindex_lower, linear_indices = true)
    A_indices = A_ind
    
    A_pars, S_pars = get_partition(A_indices, S_indices)
    
    # dimension of undirected and directed parameters
    q = size(start_test_mean)
    q_undirected = length(S_pars)
    q_directed = length(A_pars)
    
    A_indices_linear = A_indices[A_pars]
    A_indices_cartesian = linear2cartesian.(A_indices_linear, [size(S)])
    
    S_indices = linear2cartesian.(S_indices, [size(S)])
    S_indices = S_indices[S_pars]
    
    # A matrix
    A_pre = zeros(n_nod, n_nod)
    S_pre = zeros(n_nod, n_nod)
    !isnothing(M_indices) ? M_pre = zeros(n_nod) : M_pre = nothing
    
    set_RAMConstants!(A_pre, S_pre, M_pre, ram_matrices.constants)
    
    A_pre = check_acyclic(A_pre, n_par, A_ind)
    
    I_A = zeros(n_nod, n_nod)
    
    size_σ = Int(0.5*(n_var^2+n_var))
    
    if !isnothing(M)
    
        @error "meanstructure not reworked"
    
    else
        
        q_mean = nothing
        M_indices = nothing
        G_μ = nothing
        G_μ_indices = nothing
    
        G = zeros(size_σ, q_undirected)
        # TODO: analyze sparsity pattern of G
    
        ∇G = nothing
    
    end

    return SNLLS(
        q_directed,
        q_undirected,
        size_σ,

        A_indices_linear,
        A_indices_cartesian,
        S_indices,
        σ_indices,

        A_pre,
        I_A,
        G,
        ∇G,

        q_mean,
        M_indices,
        G_μ,
        G_μ_indices,

        identifier
    )
end

############################################################################
### functors
############################################################################

function objective!(imply::SNLLS, par, model::AbstractSemSingle) 

    fill_matrix(
        imply.A_pre,
        imply.A_indices_linear,
        par)

    imply.I_A .= inv(I - imply.A_pre)

    fill_G!(
        imply.G,
        imply.q_undirected,
        imply.size_σ,
        imply.S_indices,
        imply.σ_indices,
        imply.I_A)

end

############################################################################
### additional functions
############################################################################

function fill_G!(G, q_undirected, size_σ, S_indices, σ_indices, I_A)

    fill!(G, zero(eltype(G)))

    for s in 1:q_undirected
        for ind in S_indices[s]
            l, k = ind[1], ind[2]
            # rows
            for r in 1:size_σ
                i, j = σ_indices[r][1], σ_indices[r][2]
                G[r, s] += I_A[i, l]*I_A[j, k]
                if l != k
                    G[r, s] += I_A[i, k]*I_A[j, l]
                end
            end
        end
    end

end

function fill_G_μ!(G_μ, q_mean, M_indices, FI_A)

    fill!(G_μ, zero(eltype(G_μ)))

    for i in 1:q_mean
        for j in M_indices[i]
            G_μ[:, i] .+= FI_A[:, j]
        end
    end

end

function fill_∇G!(∇G, q_undirected, size_σ, q_directed, S_indices, σ_indices, A_indices, I_A, q_mean::Nothing, M_indices::Nothing, n_var)

    fill!(∇G, zero(eltype(∇G)))

    for s in 1:q_undirected

        for c ∈ S_indices[s]
            l, k = c[1], c[2]

            for r in 1:size_σ
                i, j = σ_indices[r][1], σ_indices[r][2]
                t = (s-1)*size_σ + r

                for m in 1:q_directed
                    u, v = A_indices[m][1][1], A_indices[m][1][2]
                    ∇G[t, m] += I_A[i, u]*I_A[v, l]*I_A[j, k] + I_A[i, l]*I_A[j, u]*I_A[v, k]
                    if l != k
                        ∇G[t, m] += I_A[i, u]*I_A[v, k]*I_A[j, l] + I_A[i, k]*I_A[j, u]*I_A[v, l]
                    end
                end

            end

        end

    end

end

function fill_∇G!(∇G, q_undirected, size_σ, q_directed, S_indices, σ_indices, A_indices, I_A, q_mean, M_indices, n_var)

    fill!(∇G, zero(eltype(∇G)))

    size_G_rows = size_σ + Int(n_var)

    for s in 1:q_undirected
        
        for c ∈ S_indices[s]
            l, k = c[1], c[2]

            for r in 1:size_σ
                i, j = σ_indices[r][1], σ_indices[r][2]
                t = (s-1)*size_G_rows + r

                for m in 1:q_directed
                    u, v = A_indices[m][1][1], A_indices[m][1][2]
                    ∇G[t, m] += I_A[i, u]*I_A[v, l]*I_A[j, k] + I_A[i, l]*I_A[j, u]*I_A[v, k]
                    if l != k
                        ∇G[t, m] += I_A[i, u]*I_A[v, k]*I_A[j, l] + I_A[i, k]*I_A[j, u]*I_A[v, l]
                    end
                end
            end
        end
    end

    for s in 1:q_mean
        
        for k ∈ M_indices[s]

            for r in 1:Int(n_var)
                t =  (q_undirected + s - 1)*size_G_rows + size_σ + r

                for m in 1:q_directed
                    u, v = A_indices[m][1][1], A_indices[m][1][2]
                    ∇G[t, m] += I_A[r, u]*I_A[v, k]
                end
            end
        end
    end

end

############################################################################
### Pretty Printing
############################################################################

function Base.show(io::IO, struct_inst::SNLLS)
    print_type_name(io, struct_inst)
    print_field_types(io, struct_inst)
end