## Bioequivalence

Bioequivalence is a term in pharmacokinetics used to assess the expected in vivo biological equivalence of two proprietary preparations of a drug. If two products are said to be bioequivalent it means that they would be expected to be, for all intents and purposes, the same.

Birkett (2003) defined bioequivalence by stating that, "two pharmaceutical products are bioequivalent if they are pharmaceutically equivalent and their bioavailabilities (rate and extent of availability) after administration in the same molar dose are similar to such a degree that their effects, with respect to both efficacy and safety, can be expected to be essentially the same. Pharmaceutical equivalence implies the same amount of the same active substance(s), in the same dosage form, for the same route of administration and meeting the same or comparable standards."

For The World Health Organization (WHO) "two pharmaceutical products are bioequivalent if they are pharmaceutically equivalent or pharmaceutical alternatives, and their bioavailabilities, in terms of rate (Cmax and tmax) and extent of absorption (area under the curve), after administration of the same molar dose under the same conditions, are similar to such a degree that their effects can be expected to be essentially the same".

The United States Food and Drug Administration (FDA) has defined bioequivalence as, "the absence of a significant difference in the rate and extent to which the active ingredient or active moiety in pharmaceutical equivalents or pharmaceutical alternatives becomes available at the site of drug action when administered at the same molar dose under similar conditions in an appropriately designed study."

### Requirements:

* Julia ≥1.8+
* MixedModels ≥4.11
* GLM ≥1.8.2
* Metida ≥0.14.6

### Models 

Basic models for bioequivalence.

#### A

Standard fixed-effect model.

```julia
# Parallel design GLM
fit(LinearModel, @formula(var ~ formulation), data; 
contrasts = Dict(formulation => DummyCoding(base = reference)),
)


# Cross-over designs GLM
fit(LinearModel, @formula(var ~ formulation + period + sequence + subject), data; 
contrasts = Dict(formulation => DummyCoding(base = reference)),
)
```

#### Type B

Random-effect model, no covariance.

```julia
# MixedModels

fit(MixedModel, @formula(var ~ formulation + period + sequence + (1|subject)), data; 
    contrasts = Dict(formulation => DummyCoding(base = reference)),
    REML=true
)

# Metida
lmm = LMM(@formula(var~sequence+period+formulation), data;
random = VarEffect(@covstr(1|subject), SI),
contrasts = Dict(formulation => DummyCoding(base = reference)),
)

```

#### Type C

Random-effect model, with covariance (FDA reference code).

```julia
# Metida
lmm =LMM(@formula(var~sequence+period+formulation), data;
random = VarEffect(@covstr(formulation|subject), CSH),
repeated = VarEffect(@covstr(formulation|subject), DIAG),
contrasts = Dict(formulation => DummyCoding(base = reference)),
)
```

### Reference:


* EMA: [GUIDELINE ON THE INVESTIGATION OF BIOEQUIVALENCE](https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-investigation-bioequivalence-rev1_en.pdf)
* EMA: [GUIDELINE ON THE INVESTIGATION OF BIOEQUIVALENCE, Annex I](https://www.ema.europa.eu/en/documents/other/31-annex-i-statistical-analysis-methods-compatible-ema-bioequivalence-guideline_en.pdf)
* FDA: [Guidance for Industry Statistical Approaches to Establishing Bioequivalence](https://www.fda.gov/media/70958/download)
* EEU: [Regulations conducting BE studies within the framework of the EEU](https://docs.eaeunion.org/docs/ru-ru/01411942/cncd_21112016_85)
* Chow SC. Bioavailability and Bioequivalence in Drug Development. Wiley Interdiscip Rev Comput Stat. 2014;6(4):304-312. doi: 10.1002/wics.1310. PMID: 25215170; PMCID: PMC4157693.
* Lu D, Lee SL, Lionberger RA, Choi S, Adams W, Caramenico HN, Chowdhury BA, Conner DP, Katial R, Limb S, Peters JR, Yu L, Seymour S, Li BV. International Guidelines for Bioequivalence of Locally Acting Orally Inhaled Drug Products: Similarities and Differences. AAPS J. 2015 May;17(3):546-57. doi: 10.1208/s12248-015-9733-9. Epub 2015 Mar 11. PMID: 25758352; PMCID: PMC4406956.
* Davit BM, Chen ML, Conner DP, Haidar SH, Kim S, Lee CH, Lionberger RA, Makhlouf FT, Nwakama PE, Patel DT, Schuirmann DJ, Yu LX. Implementation of a reference-scaled average bioequivalence approach for highly variable generic drug products by the US Food and Drug Administration. AAPS J. 2012 Dec;14(4):915-24. doi: 10.1208/s12248-012-9406-x. Epub 2012 Sep 13. PMID: 22972221; PMCID: PMC3475857.
* Kaza, M., Sokolovskyi, A. & Rudzki, P.J. 10th Anniversary of a Two-Stage Design in Bioequivalence. Why Has it Still Not Been Implemented?. Pharm Res 37, 140 (2020). https://doi.org/10.1007/s11095-020-02871-3
* Schütz, H., Labes, D., & Fuglsang, A. (2014). Reference datasets for 2-treatment, 2-sequence, 2-period bioequivalence studies. The AAPS journal, 16(6), 1292–1297. https://doi.org/10.1208/s12248-014-9661-0
* Fuglsang A, Schütz H, Labes D. Reference datasets for bioequivalence trials in a two-group parallel design. AAPS J. 2015 Mar;17(2):400-4. doi: 10.1208/s12248-014-9704-6. Epub 2014 Dec 9. PMID: 25488055; PMCID: PMC4365103.
* Bioequivalence reference datasets: Schütz, H., Labes, D., Tomashevskiy, M. et al. Reference Datasets for Studies in a Replicate Design Intended for Average Bioequivalence with Expanding Limits. AAPS J 22, 44 (2020). https://doi.org/10.1208/s12248-020-0427-6



### See also:

* https://github.com/PharmCat/MetidaNCA.jl
* https://github.com/PharmCat/Metida.jl
* https://juliastats.org/GLM.jl/stable/
* https://juliastats.org/MixedModels.jl/stable/

* https://bebac.at/
* https://pharmcat.net/