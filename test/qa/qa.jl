using SciMLTesting, ModelingToolkitStandardLibrary, Test

run_qa(
    ModelingToolkitStandardLibrary;
    explicit_imports = true,
    api_docs_kwargs = (; rendered = true),
    aqua_kwargs = (; ambiguities = (; recursive = false)),
    ei_kwargs = (;
        # `ifelse` is IfElse's sole export but it is not declared public in IfElse.
        all_explicit_imports_are_public = (; ignore = (:ifelse,)),
        # Names accessed qualified but not declared public in their source package:
        #   ifelse -> IfElse, SConst -> Symbolics, isvariable -> ModelingToolkitBase
        all_qualified_accesses_are_public = (; ignore = (:ifelse, :SConst, :isvariable)),
    ),
)
