module SnllsSEM

    using StructuralEquationModels, LinearAlgebra
    import StructuralEquationModels: 
        get_partition, 
        get_parameter_indices, 
        linear2cartesian, 
        cartesian2linear,
        eachindex_lower, 
        set_constants!, 
        fill_matrix,
        duplication_matrix,
        set_RAMConstants!,
        check_acyclic,
        n_par,
        identifier,
        print_type_name,
        print_field_types,
        objective!

    include("imply/SNLLS.jl")
    include("loss/SNLLS.jl")

    export SNLLS, SemSNLLS

end
