using Documenter, NMRInversions, GLMakie

using DocumenterCitations

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric
)

makedocs(sitename="NMRInversions.jl",
         authors = "Aristarchos Mavridis",

         modules=[
                  NMRInversions,
                  isdefined(Base, :get_extension) ? 
                  Base.get_extension(NMRInversions, :gui_ext) :
                  NMRInversions.gui_ext
                 ],

         pages=["Overview" => "index.md",
                "Theory" => "theory.md",
                "Tutorial" => "tutorial.md",
                "Functions" => "functions.md",
                "Types and Structures" => "types_structs.md",
                "Saving data" => "savefiles.md",
                "References" => "references.md"
               ],
         checkdocs=:none,
         format = Documenter.HTML(
         assets=String["assets/citations.css"],
         ),
         plugins=[bib]
        )

deploydocs(
    repo = "github.com/aris-mav/NMRInversions.jl.git",
)
