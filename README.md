# CODE and DATA for Reconciling the potentially irreconcilable? Genotypic and phenotypic amoxicillin-clavulanate resistance in Escherichia coli

<pr> The paper relating to this code can currently be accessed at https://aac.asm.org/content/early/2020/03/18/AAC.02026-19.long </pr>
<pr> Note all genetic data available at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA540750 </pr>


### Contents
<ol>
  <li> raw_data: MIC data, files of identified genes and the original ARIBA reports (</li>
  <li> rule_based_resistance_prediction: Code and files required to recreate the rules based resistance prediction analysis </li>
  <li> model_based_resistance_prediction: Code and files required to recreate the model based resistance prediction analysis </li>
  </ol>

#### raw_data

This folder has contains the basline metadata linked to sequence files in https://www.ncbi.nlm.nih.gov/bioproject/PRJNA540750 . This includes, 

- Ampicillin and Amoxicillin-clavulanate MICs obatined using the BD Phoenix for all 976 samples (bd_phoenix_mics.csv)
- Agar dilution MICs (both EUCAST and CLSI) (261 samples) (agar_dilution_raw_mics.json)
- Agar dilution interpreted phenotypes (agar_diltion_phenotypes.csv)
- ARIBA report files (ariba_reports.tar.gz)
- ARIBA assembly files (ariba_assembly.tar.gz)
- ARIBA MLST files (ariba_mlsts.tar.gz)

#### rule_based_resistance_prediction

This folder runs the rule based resistance prediction section of the paper, which as well as making isolate specific predictions, also reconstructs images 2 and 4 from the paper.To run this code, you will need to run the **rule_based_resistance_prediction_and_analysis.ipynb** notebook which takes you through the analysis, with comments as to what each section is doing. Data required for this analysis are also present in this folder so it can be downloaded standalone

Specifically this code , written using Python 3, also requires the following python modules
- Biopython
- Numpy
- Pandas
- Matplotlib
- Seaborn
- Pylab

#### model_based_resistance_prediction

This folder runs the mixed model based prediction section of the paper. Note this section uses both python and STATA code.
Python is used (**model_interpretations_of_genetic_features.ipynb**) to first categorize elements into model components, and then create a dataset with these model components and linked MICs (model_base.csv).
This file is then imported to STATA and used as the basis for mixed modelling. Code for this is process is in (**model_fitting_and_prediction.do**) which, as well as going through the fitting proceedure also then performs prediction at the end of the model. Output of this is then passed back to a ipython notebook to produce Figure 5 in the paper (**make_model_graph.ipynb**)

To run this code you will need 
- STATA 14.2 or later
- Numpy
- Pandas
- Matplotlib

To Summarize, to run this you 
1. run model_interpretations_of_genetic_features.ipynb to create a model base
2. run model_fitting_and_prediction.do to fit the model and predict MICs
3. run make_model_graph.ipynb to reproduce figure 5


Additionally , post this code we have included the code we use for our cross validation of the model, as detailed in the paper supplement.


###### Note

The mixed model used was developed to reflect associations between resistance elements and agar-dilution MIC as accurately as possible. However, having developed the model, an obvious question is whether the fixed effect estimates could provide accurate MIC predictions. We therefore stress that the techniques could be used as basis of MIC predcition, the exact model fitted was specific to this paper. Should one wish to use this approach on their own data they will need to re-fit their own mixed model.

###### For any questions, please submit an issue on this github page.

The paper related to this data/code is currently under review, full citation will be updated once published. In the interim, 
if using this code/data, please cite

DOI: 10.1128/AAC.02026-19

<pr> Davies TJ, Stoesser N, Sheppard AE, Abuoun M, Fowler P, Swann J, Quan TP, Griffiths D, Vaughan A, Morgan M, Phan HTT, Jeffery KJ, Andersson M, Ellington MJ, Ekelund O, Woodford N, Mathers AJ, Bonomo RA, Crook DW, Peto TEA, Anjum MF, Walker AS. 2020. Reconciling the potentially irreconcilable? Genotypic and phenotypic amoxicillin-clavulanate resistance in Escherichia coli. Antimicrob Agents Chemother. <pr>
