# Molecule-Specific-Normalization
A molecule specific normalization algorithm which adopts a robust surface fitting strategy to minimize the molecular profile difference of a group of house-keeping molecules across samples.


## A schematic work flow of the MSN method
![pipeline](https://github.com/AmeniTrabelsi/Molecule-Specific-Normalization/blob/master/Work_Flow.jpg)

## How to run

1. You can run Test_SurfaceFitting for an example of comparison of MSN with other normalization methods using our mouse liver extraction LC-MS data.
2. To apply MSN. You can use the function Surf_Fit:
Example:

finaldata=Surf_Fit(Modified_data(:,group1), Modified_data(:,group2), data_MZ, data_RT, thr)

Input:
Modified_data: data to normalize.
data_MZ: M/Z data corresponding to the data.
data_RT: Retention Time data corresponding to the data.
thr: Boxplot Threshold (ex: 1.5, 2...)

Output:

finaldata: data after normalization using MSN.

## Dependencies
libPLS_1.95 is a library to compute the PLS-DA. It can be downloaded from: (http://libpls.net/download.php)
It is also uploaded with the code.
