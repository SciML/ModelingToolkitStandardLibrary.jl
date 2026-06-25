using SciMLTesting, ModelingToolkitStandardLibrary, Test

run_qa(
    ModelingToolkitStandardLibrary;
    explicit_imports = true,
    aqua_kwargs = (; ambiguities = (; recursive = false)),
    ei_kwargs = (;
        # `unwrap` is re-exported by Symbolics but owned by SymbolicUtils, and is
        # not declared public in Symbolics.
        all_explicit_imports_via_owners = (; ignore = (:unwrap,)),
        # Names explicitly imported but not declared public in their source package:
        #   unwrap     -> Symbolics, ifelse -> IfElse (its main export),
        #   getdefault -> ModelingToolkitBase
        all_explicit_imports_are_public = (; ignore = (:unwrap, :ifelse, :getdefault)),
        # Names accessed qualified but not declared public in their source package:
        #   ifelse -> IfElse, SConst -> Symbolics, depwarn -> Base,
        #   isvariable, t_nounits -> ModelingToolkitBase
        all_qualified_accesses_are_public = (;
            ignore = (:ifelse, :SConst, :depwarn, :isvariable, :t_nounits),
        ),
    ),
    # The component submodules `using ModelingToolkitBase, Symbolics, IfElse` and
    # rely on their exported names/macros (@component, @named, @variables, System,
    # Equation, Flow, ...). Making these explicit is a large, mechanical refactor
    # tracked in https://github.com/SciML/ModelingToolkitStandardLibrary.jl/issues/467.
    ei_broken = (:no_implicit_imports,),
)
