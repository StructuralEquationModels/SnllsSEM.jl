using StructuralEquationModels, LinearAlgebra

############################################################################################
### data
############################################################################################

dat = example_data("political_democracy")
# dat_missing = example_data("political_democracy_missing")
solution_lav = example_data("political_democracy_solution")

############################################################################################
### specification - Graph
############################################################################################

observed_vars = [:x1, :x2, :x3, :y1, :y2, :y3, :y4, :y5, :y6, :y7, :y8]
latent_vars = [:ind60, :dem60, :dem65]

# meanstructure
mean_labels = label.([:m1, :m2, :m3, :m4, :m5, :m6, :m7, :m4, :m5, :m6, :m7])

graph = @StenoGraph begin
    # loadings
    ind60 → fixed(1)*x1 + x2 + x3
    dem60 → fixed(1)*y1 + y2 + y3 + y4
    dem65 → fixed(1)*y5 + y6 + y7 + y8
    # latent regressions
    label(:a)*dem60 ← ind60
    dem65 ← dem60
    dem65 ← ind60
    # variances
    _(observed_vars) ↔ _(observed_vars)
    _(latent_vars) ↔ _(latent_vars)
    # covariances
    y1 ↔ y5
    y2 ↔ y4 + y6
    y3 ↔ y7
    y8 ↔ y4 + y6
    # means
    Symbol("1") → _(mean_labels).*_(observed_vars)
end

spec_mean = ParameterTable(
    latent_vars = latent_vars,
    observed_vars = observed_vars,
    graph = graph)

# sort!(spec_mean)

partable_mean = spec_mean

start_test_mean = [fill(0.5, 8); fill(0.05, 3); fill(1.0, 11); fill(0.05, 3); fill(0.05, 13)]

semoptimizer = SemOptimizerOptim

model_ls_sym = Sem(
    specification = spec_mean,
    data = dat,
    imply = RAMSymbolic,
    loss = SemWLS,
    optimizer = semoptimizer,
    meanstructure = true
)

############################################################################################
### SNLLS imply code
############################################################################################

observed = SemObservedData(data = dat, specification = spec_mean)

ind = findall(!iszero, LowerTriangular(observed.obs_cov))
s = observed.obs_cov[ind]

# compute V
if isnothing(V)
        D = duplication_matrix(observed.n_man)
        S = inv(observed.obs_cov)
        S = kron(S, S)
        V = 0.5*(D'*S*D)
end

sᵀV = transpose(s)*V