## Example

```@example beexample
using MetidaBioeq, CSV, DataFrames, CategoricalArrays;

bedf2x2 = CSV.File(joinpath(dirname(pathof(MetidaBioeq)), "..", "test", "csv",  "2x2rds1.csv")) |> DataFrame
transform!(bedf2x2, :Subj => categorical, renamecols = false)
transform!(bedf2x2, :Per => categorical, renamecols = false)
bedf2x2.logVar = log.(bedf2x2.Var)

bedf2x2x4 = CSV.File(joinpath(dirname(pathof(MetidaBioeq)), "..", "test", "csv",  "2x2x4rds1.csv")) |> DataFrame
transform!(bedf2x2x4, :Subject => categorical, renamecols = false)
transform!(bedf2x2x4, :Period => categorical, renamecols = false)

nothing; # hide
```

### Data

```@example beexample
bedf2x2[1:5, :]
```

### BE object

Simple 2x2 study. Data (`var`) not log-transgormed.

```@example beexample
    be1 = MetidaBioeq.bioequivalence(bedf2x2, 
    vars = :Var, 
    subject = :Subj, 
    formulation = :Trt, 
    period = :Per,
    sequence = :Seq, 
    reference = "R",
    design = "2x2",
    autoseq = true,
    logt = false)
```

Replicate design 2x2x4 study. Data (`logVat`) already log-transformed.

```@example beexample
    be2 = MetidaBioeq.bioequivalence(bedf2x2x4, 
    vars = :logVar, 
    subject = :Subject, 
    formulation = :Formulation, 
    period = :Period,
    sequence = :Sequence, 
    reference = "R",
    autoseq = true,
    logt = true)
```

### Estimation 

#### GLM

Estimation witn GLM (simple model).

```@example beexample
    beAglm  = MetidaBioeq.estimate(be1;  estimator = "glm", method = "A")
```

```@example beexample
    beAglm.models[1]
```


#### MixedModels

Estimation witn MixedModels.jl (method B).

```@example beexample
    beBmm  = MetidaBioeq.estimate(be2;  estimator = "mm", method = "B")
```

```@example beexample
    beBmm.models[1]
```

#### Metida

Estimation witn Metida.jl (method C).

```@example beexample
    beBmet  = MetidaBioeq.estimate(be2;  estimator = "met", method = "C")
```

```@example beexample
    beBmet.models[1]
```