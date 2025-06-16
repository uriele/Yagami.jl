```@meta
CurrentModule = Yagami
```

# Yagami

Documentation for [Yagami](https://github.com/uriele/Yagami.jl).


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




```@docs 
Yagami
```

```@autodocs
Modules = [Yagami]
```
