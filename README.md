# Intracellular-chips-quantification

ICC: Intracellular nanochip
PMT: Photomultiplier (device that converts photons from sample fluorecence into electrical signal, and ultimatelly into pixel value)
ROI: Region of Interest

This file is a user friendly MATLAB microscopy image analysis algortithm. It guides the user through all the steps needed and displays the results in 5 variables.
The user at the end of each session should get the data of these and copy them on an spread sheet.

Workflow:
It uses as input the bright field and the corresponding fluorescence images of cell-nanochips cultures from fluorescence/confocal microscopy acquisitions.
It allows to obtain, by automatic segmentation, the background mean level, and the cells -with no chips- mean autofluorescence level. 
These two values serve to determine the threshold for normalization and noise removal. Which are caused by individual acquisition PMT offset and gain parameters setting during acquisition as well as daya to day variability.
Through this algorithm one is able to obtain true ICCs signals and compare them coming from different images - specially important in experiments were acquisitions are taken at different time points- and different treatments conditions (normalization).
Manual bright field ROIs selection of ICCs enables the final output of this program. From each chip, we get the average, maximum and total sum intensity. But also, the differences between the average cytoplasmic fluorescence of each cell containing each ICC and the average cytoplasmic fluorescence of cells without a chip.
