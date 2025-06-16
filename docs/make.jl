using Yagami
using Documenter
using DocStringExtensions

DocMeta.setdocmeta!(Yagami, :DocTestSetup, :(using Yagami); recursive=true)

makedocs(;
    checkdocs = :none,
    modules=[Yagami],
    authors="Marco <m1menari@eng.ucsd.edu> and contributors",
    sitename="Yagami.jl",
    format=Documenter.HTML(;
        canonical="https://uriele.github.io/Yagami.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "YagamiCore" => "YagamiCore.md",
        "MaterialProperties" => "MaterialProperties.md",
        "RayTracing" => "RayTracing.md",
        "CurtisGodson" => "CurtisGodson.md"
    ],
)

deploydocs(;
    repo="github.com/uriele/Yagami.jl",
    devbranch="master",
)
