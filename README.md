# PNAS---Spatial-organization-of-lysosomal-exocytosis-relies-on-membrane-tension-gradients
Datasets and codes from the publication "Spatial organization of lysosomal exocytosis relies on membrane tension gradients", Lachuer et al. PNAS (2023).

The dataset are under .RData extension and the analysis is under R. If you notice any problem or have any question, please contact the authors. If you use our codes/datasets, please cite the original paper "Spatial organization of lysosomal exocytosis relies on membrane tension gradients", Lachuer et al. PNAS (2023).

1)	Dataset Exo
This dataset contains the exocytosis coordinates from 183 VAMP7-transfected WT RPE1 cells seeded on fibronectin-coated surface. Different files correspond to this condition:
Ctrl - Exocytosis coordinates: Each element of this R list corresponds to one cell. One element of the list is a matrix with the x, y and t coordinates of each exocytosis events. The x and y coordinates are in pixel units (1 pixel = 0.160µm), and the t corresponds to the frame number.
Ctrl - Acquisition time: This R vector contains the total acquisition times for all the cells in second. Each cell is imaged for a variable number of frames (specified in another file) giving the frame rate.
Ctrl - Cell masks: Each element of this R list corresponds to the mask of one cell in a matrix format.
Ctrl - Number of frames: This R vector contains the number of frames for each movie giving the frame rate.
Ctrl - Exocytosis rate: This R vector contains the exocytosis rate for each cell in event.s-1.µm-2.

The dataset also contains the exocytosis coordinates from 36 VAMP7-transfected WT RPE1 cells seeded on PLL-coated surface. Different files correspond to this condition:
PLL - Exocytosis coordinates: Each element of this R list corresponds to one cell. One element of the list is a matrix with the x, y and t coordinates of each exocytosis events. The x and y coordinates are in pixel units (1 pixel = 0.160µm), and the t corresponds to the frame number. Frame number can be converted in seconds with the other files.
PLL - Acquisition time: This R vector contains the total acquisition times for all the cells in seconds. Each cell is imaged for 1001 frames giving the frame rate.
PLL - Cell masks: Each element of this R list corresponds to the mask of one cell in a matrix format.

The dataset also contains the exocytosis coordinates from all the drugs used: Bafilomycin A1, Para-nitroblebbistatin, Methyl-β-cyclodextrin, Cytochalasin D, Golgicide A, Histamine, Hypo-osmotic shock, Nocodazole and PF-573228. The label XXX is used to designate commonly these drugs. Different files correspond to these conditions:
XXX - Exocytosis coordinates: Each element of this R list corresponds to one cell. One element of the list is a matrix with the x, y and t coordinates of each exocytosis events. The x and y coordinates are in pixel units (1 pixel = 0.160µm), and the t corresponds to the frame number. Frame number can be converted in seconds with the other files.
XXX - Acquisition time: This R vector contains the total acquisition times for all the cells in seconds. Each cell is imaged for 1001 frames giving the frame rate.
XXX - Cell masks: Each element of this R list corresponds to the mask of one cell in a matrix format.
XXX – Condition: This R vector specifies the condition for each cell: control or drug. Because all the drugs are tested using a paired design, the condition can be interpreted as “before” or “after” the treatment. Therefore, the first cell labeled as “Ctrl” corresponds to the same cell than the first one labeled “drug” but before the treatment.

The dataset also contains the exocytosis coordinates from VAMP7-transfected WT RPE1 cells seeded on ring-shaped micropatterns of different diameters. Different files correspond to this condition:
Ring micropattern - Exocytosis coordinates: Each element of this R list corresponds to one cell. One element of the list is a matrix with the x, y and t coordinates of each exocytosis events. The x and y coordinates are in pixel units (see the corresponding file for the pixel size), and the t corresponds to the frame number. Frame number can be converted in seconds with the other files.
Ring micropattern - Acquisition time: This R vector contains the total acquisition times for all the cells in seconds. Each cell is imaged for 1001 frames giving the frame rate.
Ring micropattern - Pattern parameters: This R matrix contains geometrical parameters of the micropattern: center coordinates, diameter, etc. These values are expressed in pixels unit.
Ring micropattern - Pixel size: This R vector specifies the pixel size (in µm) for each cell. This is the only dataset that does not use a constant pixel size of 0.160µm.

The dataset also contains the exocytosis coordinates from VAMP7-transfected mCh-Rab6A RPE1 cells seeded on rectangular micropatterns. Different files correspond to this condition:
Rectangular micropattern - Exocytosis coordinates: Each element of this R list corresponds to one cell. One element of the list is a matrix with the x, y and t coordinates of each exocytosis events. The x and y coordinates are in pixel units (see the corresponding file for the pixel size), and the t corresponds to the frame number. Frame number can be converted in seconds with the other files.
Rectangular micropattern - Acquisition time: This R vector contains the total acquisition times for all the cells in seconds. Each cell is imaged for 1001 frames giving the frame rate.
Rectangular micropattern - Pattern parameters: Each element of this R list contains the x and y coordinates of the 4 borders of the micropattern in pixel units of one cell.
Rectangular micropattern – PolarizationAxis: Each element of this R list contains the x and y coordinates of the nucleus and Golgi apparatus center (in this order and with pixel units) for one cell.


2)	Dataset Colocalization
This dataset contains the quantification of the colocalization between VAMP7 and different markers (Rab11, LAMP1 and LysoTracker) in fixed RPE1 cells. The dataset contains the file:
MandersCoefficients: This .txt table contains the Manders coefficients for colocalization between VAMP7 and different markers (Rab11, LAMP1 and LysoTracker).

3)	Datasets FA colocalization
The dataset contains the data for colocalization between Focal Adhesions (FAs) (paxillin staining) and exocytosis coordinates in VAMP7-transfected WT RPE1 cells seeded on fibronectin-coated surface. Different files correspond to this condition:
Ctrl - Exocytosis coordinates: Each element of this R list corresponds to one cell. One element of the list is a matrix with the x, y and t coordinates of each exocytosis events. The x and y coordinates are in pixel units (1 pixel = 0.160µm), and the t corresponds to the frame number. Frame number can be converted in seconds with the other files.
Ctrl - Cell masks: Each element of this R list corresponds to the mask of one cell in a matrix format.
Ctrl – Paxillin staning: Each element of this R list corresponds to the paxillin staining of one cell in a matrix format (16-bits image, pixel size = 0.160µm).

The dataset also contains the data for colocalization between FAs and exocytosis coordinates for the following treatments: Methyl-β-cyclodextrin and Hypo-osmotic shock. The label XXX is used to designate commonly these treatments. Different files correspond to these conditions:
XXX - Exocytosis coordinates: Each element of this R list corresponds to one cell. One element of the list is a matrix with the x, y and t coordinates of each exocytosis events. The x and y coordinates are in pixel units (1 pixel = 0.160µm), and the t corresponds to the frame number. Frame number can be converted in seconds with the other files.
XXX - Cell masks: Each element of this R list corresponds to the mask of one cell in a matrix format.
XXX – Paxillin staning: Each element of this R list corresponds to the paxillin staining of one cell in a matrix format (16-bits image, pixel size = 0.160µm).
XXX – Condition: This R vector specifies the condition for each cell: control or drug. Because all the drugs are tested using a paired design, the condition can be interpreted as “before” or “after” the treatment. Therefore, the first cell labeled as “Ctrl” corresponds to the same cell than the first one labeled “drug” but before the treatment.

4)	Dataset Syntaxin
This dataset contains the coordinates of segmented STX3/4 spots in WT RPE1 cells seeded on ring-shaped micropattern (diameter 37µm). In this dataset, the coordinates has been normalized in order to place the micropattern center at the position (x=1,y=1) and with a radius r=1. The XXX label is used to designate commonly STX3 and STX4 staining. The dataset contains the files:
XXX – coordinates: Each element of this R list contains the normalized coordinates of STX3 or STX4 spots.

5)	Dataset Flipper Drugs
This dataset contains the Flipper-TR fluorescence lifetime maps for WT RPE1 cells from Control, Hypo-osmotic shock and Methyl-β-cyclodextrin conditions. The label XXX is used to designate commonly these conditions. The dataset contains the files:
XXX – Lifetimes: Each element of this R list contains the Flipper-TR fluorescence lifetime (in ns) of one segmented cell (pixel size = 0.087µm).

6)	Dataset Flipper Fibronectin and PLL
This dataset contains the Flipper-TR fluorescence lifetime maps for WT RPE1 cells seeded either on fibronectin- or PLL-coated surface. The label XXX is used to designate commonly these coatings. The dataset contains the files:
XXX – Lifetimes: Each element of this R list contains the Flipper-TR fluorescence lifetime (in ns) of one segmented cell (pixel size = 0.087µm).

7)	Dataset Flipper Micropattern
This dataset contains the Flipper-TR fluorescence lifetime maps for WT RPE1 cells seeded on ring-shaped micropatterns with different diameters. The dataset contains the files:
Ring micropattern – Lifetimes: Each element of this R list contains the Flipper-TR fluorescence lifetime (in ns) of one segmented cell (pixel size = 0.087µm).
Ring micropattern - Pattern parameters: This R matrix contains geometrical parameters of the micropattern: center coordinates and diameter. These values are expressed in pixels unit.

8)	Codes:
Datasets can be analyzed thanks to different R codes:
Functions: This R code contains all the functions used in the other codes presented here or for alternative analysis used in our article.
CodeSpatioTemporalAnalysis: This R code uses files from the “Dataset Exo” to do various spatiotemporal analysis including the correlation between the Area Under the Curve of the Ripley’s K function and exocytosis frequency, NND analysis, spatial/temporal/spatiotemporal Ripley’s K functions, distances to cell borders, Fourier analysis and anisotropy analysis.
SpatialAnalysisTwoSamples: This R code uses files from the “Dataset Exo” to compare the spatial structures of exocytosis events between two conditions in a paired design. The code requires different files directly present in the dataset or easily obtained from the dataset. The code compares the exocytosis rate and the Ripley’s K function before and after a treatment.
SpatioTemporalAnalysisRingShapedMicropattern: This R code uses files from the “Dataset Exo” to analyze the spatiotemporal organization of exocytosis events for cell seeded on ring-shaped micropatterns. The code requires different files directly present in the dataset or easily obtained from the dataset. The analysis includes NND analysis, map of exocytosis intensity, density of events in the radial and angular directions using KDE, polarization analysis and Ripley’s K function.
CodeExoFACoupling: This R code uses files from the “Datasets FA colocalization” to compute a colocalization index between exocytosis events and FA staining. The code requires different files directly present in the dataset.
FLIManalysisTwoSamples: This R codes uses files from the “Dataset Flipper Drugs” or “Dataset Flipper Fibronectin and PLL” to compare the fluorescence lifetimes between two conditions. The code requires different files directly present in the dataset. The code computes average lifetimes, standard-deviation, coefficient of variation, lifetime distributions and Moran’s I index. 
