# Single_patient_DMR
Set of functions to perform DMR identification in single patients.

Basically, we use the z-score to test for differential methylation at individual CpGs. Then, we aggregate individual score into DMR score using the Empirical Brown's method described in Poole, William et al. “Combining dependent P-values with an empirical adaptation of Brown's method.” Bioinformatics (Oxford, England) vol. 32,17 (2016): i430-i436. doi:10.1093/bioinformatics/btw438.

The code presented here was used in the following publication: 

Grolaux, R., Hardy, A., Olsen, C. et al. Identification of differentially methylated regions in rare diseases from a single-patient perspective. Clin Epigenet 14, 174 (2022). https://doi.org/10.1186/s13148-022-01403-7


