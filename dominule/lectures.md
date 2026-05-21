
# 1

Molfile

Maximum common substructure

Morgan algorithm
--> canonical atom numbering
--> canonical SMILES generation
--> calculation of ECFP fingerprints

### Formats
- SMILES		(1986)
- SMARTS		--> pro substructure searching (all SMILES legal in SMARTS :)), spíš na vyhledávání
- InChI		= IUPAC international chemical identifier
- InChIKey	= fixed length with 27 characters (hashing, doesn't work backwards)


## Molecular fingerprints ##
Structural keys
Morgan Fingerprints	- ECFP, klasický cirkulární přístup (radius
RDKit fingerprint	= pouze pth-based


to know:
- convert mol to SMILES and back
- Morgan algorithm
- Morgan fingerprints



# 2

Chemical space = all possible molecules (drug-like: 10^33 - 10^60 (mol.weiht <500Da), fragment-like: 10^7 (mo.weight <160Da)

## Databases:
CAS		= chemical abstrakt service = database (human curated) - 290 mil látek
DrugBank	= 15000 látek
ZINC		= vše, co lze koupit - dobrý začátek pro VS
Enamine		= i molekuly, ke kterým by se mělo dát snadno dostat (cca 100dolarů na látku, ale není úplně jistý, že to uvaří) - největší prostor

# VS
0) similarity search - nejjednodušší (lze zpracovat 10^9)
1) pharmacophore search
2) QSAR
3) Docking (max 10^6)
--> 10 - 100 molecules for lab

VS methods:
# Structure based (máme 3d strukturu targetu)	--> docking
# Ligand based		--> similarity search, pharmacophore search, QSAR

X HTS - experimental - max 10^6 (za týden max 10^5)

## SIMILARITY SEARCH
- scaffold = substruktura (třeba ta aktivní) - např. benzen od toluenu, "hledáme aktivní scaffold"
1) fingerprint calculation
2) similarity coeficient = how much similar? - the most known = Tanimoto coefficient (other: simple matching, Tversky, Pearson...)

!fingerprint affects the value of coefficient!

FILTROVÁNÍ:
- radši explicitní vodíky a soli odendat
- více tautomerů - radši se zeptat toho člověka, co s tím
- use InCHIKey to remove duplicates --> list(set(inchikeys))
- use rather SW: Papyrus structure pipeline or cHEMBL
- max 2 chirální centra (pokud je řešíme)
- Rule of 5 (Lipinského pravidlo)
- diversity filters (některé látky samy sobě hodně podobné) --> computational cost reduction
--> až nyní lze provést similarity search


Lead-like
= menší látky, to hledáme a potom optimalizujeme
X
Drug-like
= optimalizovaná drig-like molekula


"hot" molecules = velmi reaktivní - pro reakce (dostupné např v ZINC) - když si generujeme vlastní knihovnu
"promiscuous inhibitor" = aktivní i vůči dalším proteinům v těle - nejsou dostatečně specifické :(
"PAINS" = make problems in the essays (they are fluorescent, they denaturate proteins, ...)


Performance metrics - to evaluate VS
- precision = TP /(TP+FP)
- rest for Zk?
- early recognition problem = jak rychle na začátku roste ROC křivka - chceme co nejrychlejší růst na začátku (AUC ukazuje pořád stejně = nevypovídající)
- enrichment factor (EF)

Data fusion
- similarity fusion = do similarity search with more different fingerprints of the same molecule (higher difference)
- group fusion = more molecules with the same type of fingerprint

Aggregation rank

to know:
!drug-like space is bigger than fragment-like (you can combine fragment
!asymmetric similarity coefficients take into account: only ON bits (common and only X-ON and Y-ON)
?twersky coefficient
! AUC is not suffient alone because it does not encounter the enrichment factor and doesn't show the ROC curve



# 3 

## Target identification
1)  knowledge of which gene cause illness
2)  target identification
3)  validation that target is druggable
4) VS

3 categories of genes:
- disease genes
- disease modifying genes - most important
- druggable genes - subset of /\


phenotypic drug discovery:	you know the drug before you know the target		- lower throughput
X
target-based drug discovery:	you know target and search for drug (more often)	- drugs may act on more targets but faster in general


## Protein targets
- GPCR, enzymes (inhibition), transporters and ion canals

- enzymes inhibitors:
---competitive
---allosteric (non-competitive)
---irreversible = nevratné poškození enzymu = "degraders" (e.g., PROTACs = degrader with E3 ligase binding domain for ubiquitynylation)


## Lead
- sources:	VS, HTS, fragment-based screening, natural products, drug repurposing

## HTS
- primary screening = 30 000 - 100 000 compounds per day
- pipeline:	assay development, compound selection --> primary screening --> hit confirmation --> potency determination --> confirmed hits --> lead series --> drug candidate
- confirmed hits we can use as training set for QSAR
- hit = 3std devs from average (e.g., how much get the bioluminiscence lower)
- not very suitable for cells
- not suitable for heterogeneous assays (e.g., filtering and washing = too much work and problematic for HTS)
- rather FRET than SPA

- which substances to use? - design library with diversity - to have greater range
- e.g., Enamine (commercial), Natural products, HTOS
- databases:	PubChem BioAssays (noisy, requires preprocessing) or ChEMBL (more curated)


### BIOCHEMIC ASSAYS - cell-free
= kinetic assays --> Kd, Kon, Koff
- can distinguish binders and non-binders
- can not distinguish antagonist and agonists

RADIOLIGANDS
- radioligandy:		saturation / competition / kinetic (how fast to bind)
- radioligand s disociační konstantou (Kd) --> we get relevant constant for our ligand ...přímo úměrné vazbě ligandu

SCINTILLATION
- scintillant emits light when excited --> we induce it by other molecule and then detect this luminescence from scintillant
- přímo úměrné vazbě
- Scintillation Proximity Essay = SPA
- problem for HTS


### CELL-BASED ASSAYS

REPORTER GENE - LUMINISCENCE
- we can see when it is activated or not in the cell
- e.g., GFP gene or Luciferase gene in světluška
- associated with our needed promotor
- bioluminescence - PROTACs (chimeras) cause degradation of the protein --> decrease of luminiscence

to know
- induced fit doesn't necessarily induce complex changes in the target
- advantages of biochemical assays for HTS


# --- 4 --- 

PD = pharmacodynamics	- effects on the body (affinity and efficacy)
PK = pharmacokinetics	- effects on the drug itself (ADMET)

## PD
### AFFINITA
- the lower (nM, microM), the better
- Kd, Ki

### POTENCY
= amount of drug needed for response
- the lower (nM, microM), the better
- EC50 (effective concentration), IC50 (inhibitory concentration) --> log
- agonists:	full / partial / reverse

### EFFICACY
= ability of the drug to elicit the response (the maximum value)
- subset of potency

## PK
ADME = absorption, distribution, metabolism (liver), excretion

### Absorption
- Lipinski
- correct balance of water vs fat solubility

### Distribution
- water-soluble tends to be in the blood
- fat-soluble tends to concentrate in fatty tissues
- BBB (blood-brain-barrier)

### Metabolism
- Phase I
    - oxidation - Cytochrome P450 (e.g., isoform 3A4 big binding pocket &rarr; for almost everything)
    - adding polar groups, hydroxylation
    - metabolite can be toxic
    - metabolite from prodrug will be drug - yay!
    - &rarr; active metabolites
- Phase II
  - conjugation 
  - transform active metabolites to not that active

### Excretion
- mostly kidneys

### Tools for modelling AMDE
- pkCSM
- SwissADME
- ...

## Preclinical phase
- therapeutic index
- therapeutic window (between efficacy and toxicity) - the bigger window, the better
- often: off-target effects




to know:
- definition of Kd
- antagonist vs reverse agonist (antagonists block agonists when agonists present)
- know how do affinity and potency move - know the definitions
- definition of prodrugs
- pharmacodynamics - sth coefficient, antagonism, agonism
- pharmacokinetics - ADME (tox)


# --- 5 ---
# Molecular descriptors

### Experimental
* logP = log(non-polar/polar)
* logS = solubility in water - depends on pH (oral bioactivity logS > -4)
  * cancer tissue - different pH

### Calculated
* 0D - molecular weight
* 1D - fingerprints
* 2D - the whole molecule - e.g., in graph
* 3D - 3d conformation

2D descriptors:
* adjacency matrix - how atoms are connected
* topological indices - Wiener index (sth with distance matrix)
* burden matrix - atomic numbers on diagonal 

3d descriptors (3d matrices with eucledian distances):
* radius of gyration - axis, look for distances between axis and atoms (similar to Wiener index) - it says how the molecule is compact
* SAS - solvent accessible surface
* PSA - polar surface area
* TPSA - topological polar surface area (very correlated with PSA)
* etc. in RDKit (200+)

- quantum physics descriptors - package MOPAC (semi-empirical quantum chemistry)
- descriptors for atoms


## Preprocessing
- normalization - minmax scaling
- standardization - classical z-value
- distances
  * manhattan
  * tanimoto
  * mahalanobis (takes into account correlation between features and the variability))
- covariance, correlation

## Diversity selection
...clustering, dissimilarity-based metrics

* chemical diversity measures:
    * pair similarity (e.g., average closest neighbor distance)
    * based on coverage (e.g., average occupancy)


* clustering (e.g., Butina in RDKit)
  * dissimilarity-based metrics - DBCS (heuristic, easier than calculating distance matrix)
    * then, from each cluster, choose the most similar molecules
  * sphere exclusion algorithm


* cell-based methods - terrible number of dimensions - too bad


## Dimensionality reduction
- linear - PCA
- manifold methods
  - t-SNE (just tries to isolate clusters, local)
  - UMAP (more global, captures the whole structure of the space)

slido:
- easiest descriptors to calculate: clogP, TPSA, fraction of sp3 atoms, ECFP
- harder: PSA, logP, radius of gyration
- min-max scaling and standard scaling relies on an initial set of molecules
- correlation --> covariance either positive or negative
- Enamine contains building blocks - we can filter molecules using these blocks


# --- 6 ---
# METHODS OF VIRTUAL SCREENING
1) pharmacophore search
2) QSAR
3) docking
4) HTS
* (molecular dynamics, FEP)

## Pharmacophore search
- interaction features
- compare compound 3D structures in view of preferences for specific molecular interaction
- pharmacophoric features (potential pharmacophore points PPPs):
  - H-bond acceptor/donor
  - anionic/cationic
  - hydrophobic/hydrophilic
  - aromatic (they are hydrophobic)
  - halogenic
- pick only relevant (in practice about 4-5 PPPs):
  - if we already have active molecules, we align them and look for the consensus

1) database creation: we have to sample all ligands to 3d (generate rotamers)
2) database search --> score
3) hit list (--> analysis)


### Bioisosteres
- structural moiety (= funkční skupina) with similar shape and function
- mimic the original molecule while improving some properties (solubility, ...)
- e.g., carbonyl group, benzene group, ...
- SwissBioisostere database
- add halogen to 

### Matched Molecular Pairs (MMPs)
- MMP analysis
- same "key" but variable "value" fragments
- we search for molecules that would have higher activity or etc.
--> identification of activity cliffs
- but tanimoto stays the same

### Scaffold hopping
- for lead optimization
- we search for new chemotypes with the same activity
- e.g., from diazepam -> zopiclone -> zolpidem
- we search for analogs to understand the mechanism of interaction

SCAFFOLD - levels of abstractions:
1) molecule on its own
2) BM scaffold - "the main core"
3) CSK
4) RCS = reduced cyclic skeleton ("I have two rings somehow bound together")

Scaffold analysis SW
- ScaffoldGraph - builds a scaffold tree


### Pharmacophore keys/fingerprints
- take eg three PPPs that are closer together
- possible eg 3point features:
  - DDD, DDA, ..., DAR, ...
  - D = donor, A = acceptor, R = ring, H = hydrophobic, N = negative, P = positive
- fingerprint = vector of bins
  - all combinations of available PPPs of the molecule to the fingerprint
  - it is more sparse than normal fingerprint
  - 37 000 bits for 11 bins ...heavy :(
- good for similarity seach
- 2point are not specific enough
- we need enough conformations!
  - or: you can instead count the number of bounds between the PPPs - better way :)


PPPs:
- an atom containing alone electron pair
- the central cacrbon of an isopropyl group attached to pyridine (lipophilic)

non-PPPs:
- carbon atom in fused heterocycle
- a double bond between carbon and nitrogen

potential acceptors of hydrogen bonds:
- N in arometic heterocycle
- hydroxyls
- ketones
- a little also a rich-electron rings



# --- 7 ---
# SAR
= Structure Activity Relationship (SPR - property)
- qualitative method
- hit to lead phase (or lead optimization phase)
- e.g., benzodiazepines (BZDs): diazepam --> flurazepam

QSAR
= Quantitative Structure Activity Relationship (QSPR)
- quantitative method
- ML - fitting a model to the data
- OECD organisation - guidance documents for QSPR (whether new coloring is not similar to already known harmful substance)
- REACH regulations and legislative for protection of human health


ML
- bias = training error
  - too high = underfitting (too simple model)
- variance = test error
  - too high = overfitting
- ideally low both bias and variance
- regularization - we want to minimize the loss (therefore to penalize the higher weights)
-> bias-variance trade-off - lambda parameter:
    - lasso -> to zeros
    - ridge -> only to small values
    - elastic net - combines both above

decision trees
- pruning
- how to tune for bias-variance trade-off
- XGBoost


# --- 8 ---
# Early ADME profiling methods and tools

- "fail-early" strategy
- small datasets with large ooverlaps
- lack of knowledge in the mode of action --> noisy
- approved drugs = small molecules + biologicals (e.g., vaccines, antibodies)
- 90% -> 75% failuires
- failures because of ADME and poor efficacy
- toxicities tested separately (different organs and so)


- read-across - expert driven technical approach - grouping of substances
- QSAR - data-driven


- lead - over 5000 tools to check whether the lead can be used
- OECD - QED - toxic compounds database? - there is an official structural alert list
  - or ochem.eu/alerts platform
  - wim thinks qed is bad metric


### QSAR:
- naive QSAR: Single task learning (STL)
- Multi-task learning (MTL) - parallel - usually NN, XGBoost, ...
  - starts to be efficient at low data size
  - each task is independent
  - compute and unified embeddings for every task
  - regression -  l1 (lasso) or l21 regularization - we kill some descriptors
    - l1 - tasks should share the same value of l1 parameter
  - evaluation = evaraged RMSD
- Feature Net (FT) - sequential


- applicability domain = uncertainty of the model (if the predicted value is relevant or not)


- private datasets contain PROTAGs: QED = zero
- private - tend to be structurally related to known molecules (with modifications)
- potency is easier to model than ADMET


# --- 9 ---

# ...

### Note:
- Q2 is getting higher only when we add useful descriptors
- best performance metrics: MSE, R2, Q2

# Novelty detection

# Reinforcement learning - Tsetlin Machine



# --- 10 ---

- OECD principles
  - one of the principles: defined AD:

## Applicability Domain
= (cca) a set of datapoints that we can predict for
### Novelty detection - types of novelty detection:
- range-based methods - novel object = exceeds the range of at least one descriptor
- geometrical methods
- distance-based methods
- probability density function methods - na okraji v PCA grafu - outliers?
  - very good (e.g. KNN)
  - local density:
  - we look on neighbour of the neighbour (so we take into account the density - which distance-metrics do not take into account)
  - kernel density estimation - somehow probability of the input data
    - for every data point we select a gaussian and then we sum all the gaussians to get the density estimation (viz grafík v ppt)
    - KDE plot in seaborn - we can see the dense and less dense regions (from above)
- similarity-based methods (viz ADM na FITu)
  - classic similarity using fingeprints
  - e.g., tanimoto distance on morgan fingerprints

### Confidence estimation
- calibration curve (for binary classification)
  - axis x - mean predicted proba (division into bins)
  - axis y - fraction of positives
  - viz ppt --> the model will have more FN
  - the model should be better calibrated
- how to get posterior probabilities?
  - KNN
  - RN - %
  - NN - sigmoid? or softmax
- next level - conformal prediction
  - = we add confidence interval (for regression)



- TCP and ICP
  - calibration set (20) + proper training set (60) = training set (80%)
  - test set (20%)
  - we set confidence level (1 - alpha)
    --> "if we get 0.56 confidence for the prediction, we can say it is not reliable prediction"
    --> "we get 0.94 confidence for active class, we say it is reliable"
  - p-value is different from the statistical one
- CP variants
  - e.g., Mondrian conformal prediction (MCP), CCP
  - suitable for imbalanced dataset
- "empty and both" interpretation
  - both = it is on decision boundary
  - empty = it is outside applicability domain

- TCP = transductive - online (the model is retrained every time we get a new data point)
- ICP = inductive - offline (suitable for QSAR modeling)

### Validity and Efficiency
Validity = the true value will be within the confidence region (1-alpha)

Efficiency = how specific the predictions are
  - more efficient = low number of "both" and "empty" predictions
  - we can tune on efficiency - e.g., use it in cross-validation (CCP - cross-conformal prediction)

Conformal prediction - model agnostic :) - thats good


# Last Lecture

## Explainable AI (EX-AI)
* != how NN or RF works in general
* = how the trianed model works in a comprehensible way
  * why was the sample classified as active?
  * which features are important for the prediction?
  * how this feature affects the prediction - in which way (does the residue take part in binding or sth else?)

* explainability != interpretability
* interpretability = we can interpret what the model does (for a small model, e.g., small decision tree)
* explainability - we deal with model that we can not understand this way
* explainability = someone has to explain the model's inner working for us (e.g., NN)


* tasks:
  * understanding the model (diagnostics and trust for the model)
  * understanding the modelled process (scientific discovery, molecule optimization)


* examples:
  * AND operation - A=1, B=0 --> 0 (depends only on B value)
  * linear reggression - y = a*X + b (depends on a*X)

### SHAP
- popular method
- Shapley values - from game theory ---> X SHAP values

Shapley values:

* Fairness properties: efficiency, symmetry, dummy (no impact --> gets no reward), additivity
   * additivity = no matter if you divide the reward after first game or aafter all games, the result whould be same

game:
player ---> reward

model:
features ---> output

* Shapley value = sum of coalitions[ (S with J) - (S without j) ]

  * ...coalition = S

* Shapley value = average contribution (of a feature value to the prediction) in different coalitions
* Shapley value != the difference in prediction if we remove the feature from the model


SHAP values:

= SHAPley values of a conditional expectation function of the original model

SHAP = feature attribution method
* also SHAP = additive = if we sum all SHAP values, we get the final prediction
* both for classification and regression
* model-agnostic (works for any model)
* local explanation method (for a single prediction)
* interpretability - "global" explanation method (for the whole model) - we can average SHAP values for all - still not a precise global but used that way
* suitable for tabular data
* mathematical guarantee (from Shapley values)
* slow (especially for svm) <--- you need to sample from the training data distribution for those features that are not in current coalition
* limitation: prone to unrealistic data instances
* features should be independent
