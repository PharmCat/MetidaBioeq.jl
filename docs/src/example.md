## Example

```@example beexample
using MetidaBioeq, CSV, DataFrames, CategoricalArrays;

bedf2x2 = CSV.File(joinpath(dirname(pathof(MetidaBioeq)), "..", "test", "csv",  "2x2rds1.csv")) |> DataFrame
transform!(bedf2x2, :Subject => categorical, renamecols = false)
transform!(bedf2x2, :Period => categorical, renamecols = false)
bedf2x2.logVar = log.(bedf2x2.Var)

nothing; # hide
```

### Data

```@example beexample
bedf2x2[1:5, :]
```

### BE object

```@example beexample
    be2 = MetidaBioeq.bioequivalence(bedf2x2, 
    vars = :Var, 
    subject = :Subject, 
    formulation = :Formulation, 
    period = :Period,
    sequence = :Sequence, 
    reference = "R",
    design = "2x2",
    autoseq = true,
    logt = false)
```

### Estimation 

#### GLM

```@example beexample
    beAglm  = MetidaBioeq.result(be2;  estimator = "glm", method = "A")
```

```@example beexample
    beAglm.models[1]
```


#### MixedModels

```@example beexample
    beBmm  = MetidaBioeq.result(be2;  estimator = "mm", method = "B")
```

```@example beexample
    beBmm.models[1]
```

#### Metida

```@example beexample
    beBmet  = MetidaBioeq.result(be2;  estimator = "met", method = "B")
```

```@example beexample
    beBmet.models[1]
```