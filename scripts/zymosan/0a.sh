#Requirements: Local instance of HOMER (http://homer.ucsd.edu/homer/index.html), installation of mouse promotors and mm10 genome.

findMotifsGenome.pl GSE161765_DA_peaks_MEM_VS_BL_increase.bed mm10 MEM_VS_BL_increase/ -size 200 -mask -p 8
findMotifsGenome.pl GSE161765_DA_peaks_MEM_VS_BL_decrease.bed mm10 MEM_VS_BL_decrease/ -size 200 -mask -p 8
findMotifsGenome.pl GSE161765_DA_peaks_MEM_VS_BM_increase.bed mm10 MEM_VS_BM_increase/ -size 200 -mask -p 8
findMotifsGenome.pl GSE161765_DA_peaks_MEM_VS_BM_decrease.bed mm10 MEM_VS_BM_decrease/ -size 200 -mask -p 8
findMotifsGenome.pl GSE161765_DA_peaks_AP_VS_BL_increase.bed mm10 AP_VS_BL_increase/ -size 200 -mask -p 8
findMotifsGenome.pl GSE161765_DA_peaks_AP_VS_BL_decrease.bed mm10 AP_VS_BL_decrease/ -size 200 -mask -p 8
findMotifsGenome.pl GSE161765_DA_peaks_AP_VS_BM_increase.bed mm10 AP_VS_BM_increase/ -size 200 -mask -p 8
findMotifsGenome.pl GSE161765_DA_peaks_AP_VS_BM_decrease.bed mm10 AP_VS_BM_decrease/ -size 200 -mask -p 8
findMotifsGenome.pl GSE161765_DA_peaks_BL_VS_BM_increase.bed mm10 BL_VS_BM_increase/ -size 200 -mask -p 8
findMotifsGenome.pl GSE161765_DA_peaks_BL_VS_BM_decrease.bed mm10 BL_VS_BM_decrease/ -size 200 -mask -p 8