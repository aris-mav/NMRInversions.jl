using Documenter, NMRInversions, GLMakie

makedocs(sitename="NMRInversions.jl",
         authors = "Aristarchos Mavridis",

         modules=[
                  NMRInversions,
                  isdefined(Base, :get_extension) ? 
                  Base.get_extension(NMRInversions, :gui_ext) :
                  NMRInversions.gui_ext
                 ],

         pages=["Overview" => "index.md",
                "Tutorial" => "tutorial.md",
                "Functions" => "functions.md",
                "Types and Structures" => "types_structs.md",
                "Saving data" => "savefiles.md",
                #="Theory" => "theory.md",=#
               ],
         checkdocs=:none
        )

deploydocs(
    repo = "github.com/aris-mav/NMRInversions.jl.git",
)
