# Yagami

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://uriele.github.io/Yagami.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://uriele.github.io/Yagami.jl/dev/)
[![Build Status](https://github.com/uriele/Yagami.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/uriele/Yagami.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/uriele/Yagami.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/uriele/Yagami.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)


Package for ray tracing and retrieval for the CAIRT project.

## Structure
The package is structured as follow:\

```mermaid
flowchart TD

subgraph Yagami["📂 Yagami Package"]
  
  subgraph YagamiCore["📦 Yagami Core"]
    CoreDescription["Common functions,<br/>constants, and utilities<br/>shared across modules"]
  end

  subgraph Modules
    direction LR
    MaterialProperties["📚 Material Properties<br/>(Refractive indices,<br/>atmosphere discretization)"]
    RayTracing["🔭 RayTracing<br/>(Finding intersections<br/>with Earth)"]
    CurtisGodson["📐 CurtisGodson<br/>(Integral computations<br/>for direct model)"]
  end

  YagamiCore --> MaterialProperties
  YagamiCore --> RayTracing
  YagamiCore --> CurtisGodson
end
```



Python scripts for the testing the material properties are adapted from the [Refractive Index Database](https://github.com/polyanskiy/refractiveindex.info-database?tab=readme-ov-file) created by [Mikhail Polyanskiy](https://www.bnl.gov/staff/polyanskiy).[^1]




# References

[^1]: Polyanskiy, Mikhail N. ["Refractiveindex. info database of optical constants."](https://www.nature.com/articles/s41597-023-02898-2) Scientific Data 11.1 (2024): 94.


