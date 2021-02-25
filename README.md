# COVID19 School Dynamics

Project name is a work in progress!

### Getting started

##### "Installation"

These steps prevent polluting and breaking packages on your machine.

1. Move the enclosing folder, `COVID19SchoolReopening`, to a convenient location.

2. Start Julia.

3. Type `]activate .` to switch to the environment for this project. If this was successful Julia prompt should now read `(COVID19SchoolReopening) pkg>`.

4. Type `]instantiate` to reproduce the environment defined in `Manifest.toml`.

Now the package code and examples will work whenever the `COVID19SchoolReopening` environment is active.

##### Using the code

**In a REPL:**

1. Start julia

2. Type `]activate .`

3. The usual `using COVID19SchoolReopening` will work.

**In a Jupyter Notebook:**

1. Set the working directory to the enclosing folder, e.g. `path/to/COVID19SchoolReopening` using `cd` in the command line or `cd()` in Julia.

2. Start the notebook server as usual.

The notebooks in this project should just work using this method.
The following workflow should guarantee that everything works:

```julia
using IJulia

# assuming you're currently in the folder COVID19SchoolReopening
notebook(detached=true, dir=".")

# otherwise
notebook(detached=true, dir="path/to/COVID19SchoolReopening")
```

##### File structure

- `src`: Contains code used to generate ODE models.
- `experiments`: Notes with Julia code I used to play around with the models. Everything in here, including the notebooks, is generated from the `*.jmd` files using Weave.
