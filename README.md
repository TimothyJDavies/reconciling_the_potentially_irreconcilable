# CODE and DATA for Reconciling the potentially irreconcilable? Genotypic and phenotypic amoxicillin-clavulanate resistance in Escherichia coli

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

This folder runs the rule based resistance prediction section of the paper, which as well as making isolate specific predictions, also reconstructs images 2 and 4 from the paper.To run this code, you will need to run the rule_based_resistance_prediction_and_analysis.ipynb notebook which takes you through the analysis, with comments as to what each section is doing. Data required for this analysis are also present in this folder so it can be downloaded standalone

Specifically this code , written using Python 3, also requires the following python modules
- Biopython
- Numpy
- Pandas
- Matplotlib
- Seaborn
- Pylab

#### model_based_resistance_prediction
T

###### For any questions, please submit an issue on this github page.

The paper related to this data/code is currently under review, full citation will be updated once published. In the interim, 
if using this code/data, please cite

DOI: https://doi.org/10.1101/511402

<pr> Davies TJ, Stoesser N, Sheppard AE, Abuoun M, Fowler P, Swann J, Quan TP, Griffiths D, Vaughan A, Morgan M, Phan HT, Jeffery KJ, Andersson M, Ellington MJ, Ekelund O, Mathers AJ, Bonomo RA, Woodford N, Crook DW, Peto TEA, Anjum MF, Walker AS. 2019. Reconciling the potentially irreconcilable? Genotypic and phenotypic amoxicillin-clavulanate resistance in Escherichia coli. bioRxiv 511402. <pr>
