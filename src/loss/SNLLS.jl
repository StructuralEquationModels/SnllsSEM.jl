##### weighted least squares

############################################################################
### Types
############################################################################

mutable struct SemSNLLS{Vt, St, B} <: SemLossFunction
    V::Vt
    s::St
    sᵀV::B
end

############################################################################
### Constructors
############################################################################

function SemSNLLS(;
    observed,
    V = nothing,
    kwargs...)
    
    ind = findall(!iszero, LowerTriangular(observed.obs_cov))
    s = observed.obs_cov[ind]
    
    # compute V
    if isnothing(V)
        D = duplication_matrix(observed.n_man)
        S = inv(observed.obs_cov)
        S = kron(S, S)
        V = Symmetric(0.5*(D'*S*D))
    end
    
    sᵀV = transpose(s)*V

    return SemSNLLS(
        V, 
        s,
        sᵀV)
end

############################################################################
### functors
############################################################################

function objective!(semsnlls::SemSNLLS, par, model::AbstractSemSingle)

    if !isnothing(model.imply.Gc)
        s = semsnlls.s - model.imply.Gc*model.imply.c
        sᵀV = s'*semsnlls.V
    else
        sᵀV = semsnlls.sᵀV
    end
    
    outer = sᵀV*model.imply.G
    prod = Symmetric(model.imply.G'*semsnlls.V*model.imply.G)
    b = cholesky(prod; check = false)

    if !isposdef(b)
        a = zeros(size(prod)); a[diagind(a)] .= 0.001
        b = cholesky(prod + a; check = false)

        if !isposdef(b)
            return Inf
        else
            a = b\(transpose(outer))
            return -outer*a
        end

    else
        a = b\(transpose(outer))
        return -outer*a
    end
end

############################################################################
### additional functions
############################################################################

############################################################################
### Pretty Printing
############################################################################

function Base.show(io::IO, struct_inst::SemSNLLS)
    print_type_name(io, struct_inst)
    print_field_types(io, struct_inst)
end