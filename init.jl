# initialize dependencies:
using Pkg
# create project
current_dir = pwd()

if isdir(joinpath(current_dir,"env")) == false
    Pkg.generate("./env")
end

# activate project
Pkg.activate("./env")

dependencies = ["Polymake","GAP","Combinatorics","Nemo","PrettyTables","DelimitedFiles","AbstractAlgebra"];

for package in dependencies
    Pkg.add(package)
end

Pkg.instantiate()