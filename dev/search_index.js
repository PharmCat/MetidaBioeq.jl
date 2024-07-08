var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#MetidaBioeq.bioequivalence","page":"API","title":"MetidaBioeq.bioequivalence","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaBioeq.bioequivalence","category":"page"},{"location":"api/#MetidaBioeq.bioequivalence","page":"API","title":"MetidaBioeq.bioequivalence","text":"bioequivalence(data;\n    vars = nothing,\n    subject::Union{String, Symbol},\n    period::Union{String, Symbol, Nothing} = nothing,\n    formulation::Union{String, Symbol},\n    sequence::Union{String, Symbol, Nothing} = nothing,\n    stage::Union{String, Symbol, Nothing} = nothing,\n    reference::Union{String, Symbol, Nothing} = nothing,\n    design::Union{String, Symbol, Nothing} = nothing,\n    io::IO = stdout,\n    seqcheck::Bool = true,\n    designcheck::Bool = true,\n    dropcheck::Bool = true,\n    dropmissingsubj = false,\n    info::Bool = true,\n    warns::Bool = true,\n    autoseq::Bool = false,\n    logt::Bool = true)\n\nvars - variabel's column(s);\nsubject - subject's column;\nperiod - period's column;\nformulation - formulation's column;\nsequence -sequence's column;\nstage - stage's column;\nreference - reference value for formulation column;\ndesign - design: \"parallel\", \"2X2\", \"2X2X2\", \"2X2X4\", ets. (formulations X sequences X periods);\nseqcheck - check sequencs; \ndesigncheck - check design correctness;  \ndropcheck - dropuot check;\ndropmissingsubj - drop subjects with missing data;\ninfo - show information;\nwarns - show warnings;\nautoseq - try to make sequence collumn;\nlogt - if true (default) data is already log-transformed, else log() will be used.\n\n\n\n\n\n","category":"function"},{"location":"api/#MetidaBioeq.estimate","page":"API","title":"MetidaBioeq.estimate","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaBioeq.estimate","category":"page"},{"location":"api/#MetidaBioeq.estimate","page":"API","title":"MetidaBioeq.estimate","text":"estimate(be; estimator = \"auto\", method = \"auto\", supresswarn = false)\n\nmethod - Model settings.\n\nif method == \"auto\" than method A used for \"2X2\" and \"2X2X2\"  designes, method P for \"parallel\" design and method B for any other.\n\nMethods:\n\nA using GLM and model @formula(var ~ formulation + period + sequence + subject)\nB using MixedModels and model @formula(var ~ formulation + period + sequence + (1|subject)) or Metida and model @lmmformula(v ~ formulation + period + sequence, random = 1|subject:SI)\nC using Metida and model @lmmformula(v ~ formulation + period + sequence, random = formulation|subject:CSH, repeated = formulation|subject:DIAG)\nP using GLM and model @formula(var ~ formulation)\n\nestimator - Estimator settings.\n\nif estimator == \"auto\" than GLM used for \"parallel\" design; for \"2X2\" design used GLM if no droputs and MixedModels if dropout == true; for other designes with method C Metida used and MixedModel for other cases.\n\nEstimators:\n\n\"glm\" for GLM (https://juliastats.org/GLM.jl/stable/)\n\"mm\" for MixedModels (https://juliastats.org/MixedModels.jl/stable/)\n\"met\" for Metida (https://pharmcat.github.io/Metida.jl/stable/)\n\nOther autosettings:\n\nIf design is \"parallel\" estimator set as \"glm\" and method as \"P\".\n\nIf design is \"2X2\" and method is \"P\" or \"C\" than if estimator == \"glm\" method set as \"A\" and \"B\" for other estimators. \n\nIf design not \"parallel\" or \"2X2\": \n\nif method not \"A\", \"B\" or \"C\" than set as \"A\" for \"glm\" ann as B for other estimators;\n\nif estimator == \"glm\" and method == \"B\" than estimator set as \"mm\", if estimator == \"glm\" or \"mm\" and method == \"C\" than estimator set as \"met\".\n\nReference:\n\nEMA: GUIDELINE ON THE INVESTIGATION OF BIOEQUIVALENCE\n\nEMA: GUIDELINE ON THE INVESTIGATION OF BIOEQUIVALENCE, Annex I\n\n\n\n\n\n","category":"function"},{"location":"api/#MetidaBioeq.result","page":"API","title":"MetidaBioeq.result","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaBioeq.result","category":"page"},{"location":"api/#MetidaBioeq.result","page":"API","title":"MetidaBioeq.result","text":"result(beres::BEResults)\n\nReturns dataframe with bioequivalence results.\n\n\n\n\n\n","category":"function"},{"location":"api/#MetidaBioeq.makeseq","page":"API","title":"MetidaBioeq.makeseq","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaBioeq.makeseq","category":"page"},{"location":"api/#MetidaBioeq.makeseq","page":"API","title":"MetidaBioeq.makeseq","text":"makeseq(data;\n    subject = :subject,\n    period = :period,\n    formulation = :formulation)\n\nMake sequence vector from data and subject, period, formulation columns.\n\n\n\n\n\n","category":"function"},{"location":"#Bioequivalence","page":"Home","title":"Bioequivalence","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Bioequivalence is a term in pharmacokinetics used to assess the expected in vivo biological equivalence of two proprietary preparations of a drug. If two products are said to be bioequivalent it means that they would be expected to be, for all intents and purposes, the same.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Birkett (2003) defined bioequivalence by stating that, \"two pharmaceutical products are bioequivalent if they are pharmaceutically equivalent and their bioavailabilities (rate and extent of availability) after administration in the same molar dose are similar to such a degree that their effects, with respect to both efficacy and safety, can be expected to be essentially the same. Pharmaceutical equivalence implies the same amount of the same active substance(s), in the same dosage form, for the same route of administration and meeting the same or comparable standards.\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"For The World Health Organization (WHO) \"two pharmaceutical products are bioequivalent if they are pharmaceutically equivalent or pharmaceutical alternatives, and their bioavailabilities, in terms of rate (Cmax and tmax) and extent of absorption (area under the curve), after administration of the same molar dose under the same conditions, are similar to such a degree that their effects can be expected to be essentially the same\".","category":"page"},{"location":"","page":"Home","title":"Home","text":"The United States Food and Drug Administration (FDA) has defined bioequivalence as, \"the absence of a significant difference in the rate and extent to which the active ingredient or active moiety in pharmaceutical equivalents or pharmaceutical alternatives becomes available at the site of drug action when administered at the same molar dose under similar conditions in an appropriately designed study.\"","category":"page"},{"location":"#Requirements:","page":"Home","title":"Requirements:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Julia ≥1.8+\nMixedModels ≥4.11\nGLM ≥1.8.2\nMetida ≥0.15","category":"page"},{"location":"#Models","page":"Home","title":"Models","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Basic models for bioequivalence.","category":"page"},{"location":"#A","page":"Home","title":"A","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Standard fixed-effect model. For parallel design Welch's correction not used.","category":"page"},{"location":"","page":"Home","title":"Home","text":"# Parallel design GLM\nfit(LinearModel, @formula(var ~ formulation), data; \ncontrasts = Dict(formulation => DummyCoding(base = reference)),\n)\n\n\n# Cross-over designs GLM\nfit(LinearModel, @formula(var ~ formulation + period + sequence + subject), data; \ncontrasts = Dict(formulation => DummyCoding(base = reference)),\n)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Validated with:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Schütz, H., Labes, D., & Fuglsang, A. (2014). Reference datasets for 2-treatment, 2-sequence, 2-period bioequivalence studies. The AAPS journal, 16(6), 1292–1297. https://doi.org/10.1208/s12248-014-9661-0\nFuglsang, A., Schütz, H., & Labes, D. (2015). Reference datasets for bioequivalence trials in a two-group parallel design. The AAPS journal, 17(2), 400–404. https://doi.org/10.1208/s12248-014-9704-6","category":"page"},{"location":"#Type-B","page":"Home","title":"Type B","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Random-effect model, no covariance.","category":"page"},{"location":"","page":"Home","title":"Home","text":"# MixedModels\n\nfit(MixedModel, @formula(var ~ formulation + period + sequence + (1|subject)), data; \n    contrasts = Dict(formulation => DummyCoding(base = reference)),\n    REML=true\n)\n\n# Metida\nlmm = LMM(@formula(var~sequence+period+formulation), data;\nrandom = VarEffect(@covstr(1|subject), SI),\ncontrasts = Dict(formulation => DummyCoding(base = reference)),\n)\n","category":"page"},{"location":"","page":"Home","title":"Home","text":"Validated with: ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Schütz, H., Labes, D., Tomashevskiy, M., la Parra, M. G., Shitova, A., & Fuglsang, A. (2020). Reference Datasets for Studies in a Replicate Design Intended for Average Bioequivalence with Expanding Limits. The AAPS journal, 22(2), 44. https://doi.org/10.1208/s12248-020-0427-6","category":"page"},{"location":"#Type-C","page":"Home","title":"Type C","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Random-effect model, with covariance (FDA reference code).","category":"page"},{"location":"","page":"Home","title":"Home","text":"# Metida\nlmm =LMM(@formula(var~sequence+period+formulation), data;\nrandom = VarEffect(@covstr(formulation|subject), CSH),\nrepeated = VarEffect(@covstr(formulation|subject), DIAG),\ncontrasts = Dict(formulation => DummyCoding(base = reference)),\n)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Validated against SPSS model:","category":"page"},{"location":"","page":"Home","title":"Home","text":"MIXED var BY period formulation sequence subject\n  /CRITERIA=CIN(90) MXITER(200) MXSTEP(20) SCORING(2) SINGULAR(0.000000000001) HCONVERGE(0,\n    RELATIVE) LCONVERGE(0.0000000000001, RELATIVE) PCONVERGE(0, RELATIVE)\n  /FIXED=period formulation sequence | SSTYPE(3)\n  /METHOD=REML\n  /RANDOM=formulation | SUBJECT(subject) COVTYPE(CSH)\n  /REPEATED=formulation | SUBJECT(subject*period) COVTYPE(DIAG)\n  /EMMEANS=TABLES(formulation) COMPARE REFCAT(FIRST) ADJ(LSD).","category":"page"},{"location":"","page":"Home","title":"Home","text":"What complies with the FDA code (FDA Guidance for Industry: Statistical Approaches to Establishing Bioequivalence, APPENDIX F):","category":"page"},{"location":"","page":"Home","title":"Home","text":"PROC MIXED;\nCLASSES SEQ SUBJ PER TRT;\nMODEL Y = SEQ PER TRT/ DDFM=SATTERTH;\nRANDOM TRT/TYPE=FA0(2) SUB=SUBJ G;\nREPEATED/GRP=TRT SUB=SUBJ;\nESTIMATE 'T vs. R' TRT 1 -1/CL ALPHA=0.1;","category":"page"},{"location":"","page":"Home","title":"Home","text":"In the Random statement, TYPE=FA0(2) could possibly be replaced by TYPE=CSH. This guidance recommends that TYPE=UN not be used, as it could result in an invalid (i.e., not non- negative definite) estimated covariance matrix","category":"page"},{"location":"","page":"Home","title":"Home","text":"Exported SPSS results.","category":"page"},{"location":"#Reference:","page":"Home","title":"Reference:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"EMA: GUIDELINE ON THE INVESTIGATION OF BIOEQUIVALENCE\nEMA: GUIDELINE ON THE INVESTIGATION OF BIOEQUIVALENCE, Annex I\nFDA: Guidance for Industry Statistical Approaches to Establishing Bioequivalence\nEEU: Regulations conducting BE studies within the framework of the EEU\nChow SC. Bioavailability and Bioequivalence in Drug Development. Wiley Interdiscip Rev Comput Stat. 2014;6(4):304-312. doi: 10.1002/wics.1310. PMID: 25215170; PMCID: PMC4157693.\nLu D, Lee SL, Lionberger RA, Choi S, Adams W, Caramenico HN, Chowdhury BA, Conner DP, Katial R, Limb S, Peters JR, Yu L, Seymour S, Li BV. International Guidelines for Bioequivalence of Locally Acting Orally Inhaled Drug Products: Similarities and Differences. AAPS J. 2015 May;17(3):546-57. doi: 10.1208/s12248-015-9733-9. Epub 2015 Mar 11. PMID: 25758352; PMCID: PMC4406956.\nDavit BM, Chen ML, Conner DP, Haidar SH, Kim S, Lee CH, Lionberger RA, Makhlouf FT, Nwakama PE, Patel DT, Schuirmann DJ, Yu LX. Implementation of a reference-scaled average bioequivalence approach for highly variable generic drug products by the US Food and Drug Administration. AAPS J. 2012 Dec;14(4):915-24. doi: 10.1208/s12248-012-9406-x. Epub 2012 Sep 13. PMID: 22972221; PMCID: PMC3475857.\nKaza, M., Sokolovskyi, A. & Rudzki, P.J. 10th Anniversary of a Two-Stage Design in Bioequivalence. Why Has it Still Not Been Implemented?. Pharm Res 37, 140 (2020). https://doi.org/10.1007/s11095-020-02871-3\nSchütz, H., Labes, D., & Fuglsang, A. (2014). Reference datasets for 2-treatment, 2-sequence, 2-period bioequivalence studies. The AAPS journal, 16(6), 1292–1297. https://doi.org/10.1208/s12248-014-9661-0\nFuglsang A, Schütz H, Labes D. Reference datasets for bioequivalence trials in a two-group parallel design. AAPS J. 2015 Mar;17(2):400-4. doi: 10.1208/s12248-014-9704-6. Epub 2014 Dec 9. PMID: 25488055; PMCID: PMC4365103.\nBioequivalence reference datasets: Schütz, H., Labes, D., Tomashevskiy, M. et al. Reference Datasets for Studies in a Replicate Design Intended for Average Bioequivalence with Expanding Limits. AAPS J 22, 44 (2020). https://doi.org/10.1208/s12248-020-0427-6","category":"page"},{"location":"#See-also:","page":"Home","title":"See also:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"https://github.com/PharmCat/MetidaNCA.jl\nhttps://github.com/PharmCat/Metida.jl\nhttps://juliastats.org/GLM.jl/stable/\nhttps://juliastats.org/MixedModels.jl/stable/\nhttps://bebac.at/\nhttps://pharmcat.net/","category":"page"},{"location":"example/#Example","page":"Example","title":"Example","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"using MetidaBioeq, CSV, DataFrames, CategoricalArrays;\n\nbedf2x2 = CSV.File(joinpath(dirname(pathof(MetidaBioeq)), \"..\", \"test\", \"csv\",  \"2x2rds1.csv\")) |> DataFrame\ntransform!(bedf2x2, :Subj => categorical, renamecols = false)\ntransform!(bedf2x2, :Per => categorical, renamecols = false)\nbedf2x2.logVar = log.(bedf2x2.Var)\n\nbedf2x2x4 = CSV.File(joinpath(dirname(pathof(MetidaBioeq)), \"..\", \"test\", \"csv\",  \"2x2x4rds1.csv\")) |> DataFrame\ntransform!(bedf2x2x4, :Subject => categorical, renamecols = false)\ntransform!(bedf2x2x4, :Period => categorical, renamecols = false)\n\nnothing; # hide","category":"page"},{"location":"example/#Data","page":"Example","title":"Data","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"bedf2x2[1:5, :]","category":"page"},{"location":"example/#BE-object","page":"Example","title":"BE object","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"Simple 2x2 study. Data (var) not log-transgormed.","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"    be1 = MetidaBioeq.bioequivalence(bedf2x2, \n    vars = :Var, \n    subject = :Subj, \n    formulation = :Trt, \n    period = :Per,\n    sequence = :Seq, \n    reference = \"R\",\n    design = \"2x2\",\n    autoseq = true,\n    logt = false)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"Replicate design 2x2x4 study. Data (logVat) already log-transformed.","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"    be2 = MetidaBioeq.bioequivalence(bedf2x2x4, \n    vars = :logVar, \n    subject = :Subject, \n    formulation = :Formulation, \n    period = :Period,\n    sequence = :Sequence, \n    reference = \"R\",\n    autoseq = true,\n    logt = true)","category":"page"},{"location":"example/#Estimation","page":"Example","title":"Estimation","text":"","category":"section"},{"location":"example/#GLM","page":"Example","title":"GLM","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"Estimation witn GLM (simple model).","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"    beAglm  = MetidaBioeq.estimate(be1;  estimator = \"glm\", method = \"A\")","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"    beAglm.models[1]","category":"page"},{"location":"example/#MixedModels","page":"Example","title":"MixedModels","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"Estimation witn MixedModels.jl (method B).","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"    beBmm  = MetidaBioeq.estimate(be2;  estimator = \"mm\", method = \"B\")","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"    beBmm.models[1]","category":"page"},{"location":"example/#Metida","page":"Example","title":"Metida","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"Estimation witn Metida.jl (method C).","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"    beBmet  = MetidaBioeq.estimate(be2;  estimator = \"met\", method = \"C\")","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"    beBmet.models[1]","category":"page"}]
}
