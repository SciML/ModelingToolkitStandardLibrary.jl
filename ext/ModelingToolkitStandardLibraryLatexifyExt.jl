module ModelingToolkitStandardLibraryLatexifyExt

using ModelingToolkitStandardLibrary
using Latexify

@latexrecipe function f(n::AnalysisPoint)
    env --> :equation
    cdot --> false
    index --> :subscript

    return nameof(n)
end

end