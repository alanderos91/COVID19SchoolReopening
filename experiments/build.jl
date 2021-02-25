using Weave

# Select .jmd files.
jmd_files = readdir("experiments")
filter!(contains(".jmd"), jmd_files)
filenames = map(file -> splitext(file)[1], jmd_files)

for filename in filenames
    jmd_file = joinpath("experiments", filename * ".jmd")
    ipynb_file = joinpath("experiments", filename * ".ipynb")

    # Weave to HTML for easy viewing.
    weave(jmd_file)

    # Weave to a Jupyter notebook for interactive version.
    convert_doc(jmd_file, ipynb_file)
end
