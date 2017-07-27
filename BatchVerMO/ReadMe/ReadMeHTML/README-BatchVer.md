# CNMF-E (BatchVer)

CNMF-E is powerful in extracting neurons from calcium imaging data captured by endoscope, especially in longer data. However, due to the concern of photo-bleaching and other experimental limits and designs, usually our data sets consist of large amount of small(short) data. One can run CNMF-E on all the small data individually yet if you would like to track the same nuerons over days, some extra analysis would be needed. This version helps deal with this situation. Most inspirations of this version derive from ideas implemented in the normal CNMF-E. It has good scalability and manageable memory use. Please read CNMF-E principle before using this version as the documentation assumes the terminology used in normal CNMF-E. See section "Description" for method descirption, and see "Getting Started" to use this version of cnmf-e.

## Description
The package uses normal CNMF-E on two different data sets with two different "depth":

 (1) full CNMF-E on sample data---"initiation mode"---to get a representative A of this data set; 
 
 (2) simplified CNMF-E on all data individually, regressing background subtracted data onto A (representative A found in step 1) and C---"massive mode". 
 
 The following figure shows the general idea of the batched version of CNMF-E.
 
 <center><b>Figure1, General Description of Batch Version CNMF-E</center></b>
![](/Users/gushijie/Documents/MATLAB/CaImaging/cnmf-e/BatchVer/ReadMe/Overview.png)

 More detailed decription of this version is shown in Figure2 below. You may need to consult the variable summary table below alongside.
 
 <center><b>Figure2, Detailed Description of Batch Version CNMF-E</b></center>
![](/Users/gushijie/Documents/MATLAB/CaImaging/cnmf-e/BatchVer/ReadMe/July10th.png)
 
 
 <center><b>Table1, Variable Summary Table</b></center>
![](/Users/gushijie/Documents/MATLAB/CaImaging/cnmf-e/BatchVer/ReadMe/VariablesJuly10th.png)


## Prerequisites
__*Very Important*__
0. Motion corrected data

Data must not have motion. Users can tell if there is motion by inspecting the spatial footprint for each sample in the "result files". Open them as unit8 images in imageJ and paste one image on top of another (select white-zero in "paste control").

1. Basic cnmf-e

This package builds on basic cnmf-e. (A summary between the differences between the basic version and this batch version is below.) This not only means you need code from the basic version but also means that certain requirement for data of basic cnmf-e apply here as well, such as motion-corrected data.

2. Cluster-dependancy

Though this package has fairly adjustable RAM requirement provided by the flexible choice of samples, its normal requirement for 50 GB of RAM makes it cluster-dependent. Except from using small data testing your code, in which case you can set "running_on_cluster=false", you should normally set "running_on_cluster=true". Meanwhile, specify workersnum you want to have given your cluster setting in the input section. More will be discussed in section "Getting Started".

3. Data storage site accessible both on your cluster and your local machine.
You will see why in just a minute.

4. Matlab version

Matlab version may come into play in terms of whether the package will run smooth or not. The package is tested during development on MATLAB2016b for the parts run on cluster. Although, the beginning part of the package is simple input/output specification, which can be run on 2014b or onwards.

## Getting Started

This package is of two parts: (1) cnmfeBatchVer_LocalPart  (2) cnmfeBatchVer_ClusterPart
The first part is input/output specification. It runs on your local computer. The purpose of this part is to make sure sample data and experimental data are properly chosen. The variables will be saved and will be used in cnmfeBatchVer_ClusterPart, which is all the actual work of neuron extraction and output parsing.

### (1) cnmfeBatchVer_LocalPart
#### A. Important directories and parameters
1. Open the template/demp "cnmfeBatchVer_LocalPart"  __in the "BatchVer" folder__. If not, specify code directory and somthing like addpath(genpath(codeDir));

2. Fill in where the code will be __for the cluster__.

3. In InputOutput(), fill in arguments where all the directory should be with respect to __your local machine/PC__, not the cluster. This function involves an interactive process and possibly a pop-out window to choose sample files if you set 'SamplingMethod'='manual'. Alternatively, choose 'SamplingMethod'='auto'. You will only need to input a number when prompted. Type "help InputOutput()" in MATLAB command window for more infomation. Note, data types that are supported are the same as that in the original CNMF-E.

4. After this section, section1.1 helps you to change directory pre-fix from your local machine to that on cluster. Users are encouraged to look into InputOutput() and section 1.1 to change the syntax/method that easily specify directories given the user's data storage structure.

5. Section2 specifies some parameters used in normal cnmf-e. There are a few more parameters that can tuned but are not listed here. User can modify this according to their needs, but most important ones that are able to make a significant difference are listed here. One different "parameter" is the last variable, "namepattern". In cnmfe(BatchVer), after extracting each sample data, a plot of A will be saved in the working folder of the cluster. In the current release, the plot will be named with filename(namepattern). In addition, in cnmfe(BatchVer), initiation process of each sample data has a plot showing the min_PNR*min_corr curve and the seeds of neuron. This plot is intended to help user to better find a suitable min_PNR. This plot uses the same naming strategy.

6. Section3 is some parameters specifically important to BatchVer. Type help Over_Days_findAnn for explaination and help MergeAC. There is no general rule; it depends on your data. Please try some to figure out a suitable value for you. Note that the "merge_thr_2" here will be used in MergeAC, the BatchVer specific merging, while "merge_thr" in section2 is used for for normal CNMF-E neuron merging in initiation and iteration process.

7. Section4 specifies workers of matlab parallel pool. This number will depend on the cluster you use and the CPU number you request. Consult your cluster documentation on this. With this workersnum, cnmfeBatchVer_ClusterPart has functions written to create a cluster and a parpool within it. The design caters to cluster use. The functions for these are in a seperate folder "HandlingParforbyGalen" which is written by Galen Lynch.

8. Save all these variables to the outputdir with respect to your local machine, which will be read on cluster from "outputdir".

#### B. Script for cluster

Fill in the template in Part B to make a script for cluster. The example here is for clusters using SLURM for job management.
To make the script look clearer, we reproduce results here.
```
#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200000
#SBATCH -t 0-10:00
#SBATCH --time-min=0-01:00
#SBATCH -o /net/feevault/data0/shared/EmilyShijieShared/BatchResult/6938FirstFewDaysForBatch/job_%A.out
#SBATCH -e /net/feevault/data0/shared/EmilyShijieShared/BatchResult/6938FirstFewDaysForBatch/job_%A.err
cd /net/feevault/data0/shared/EmilyShijieShared/BatchResult/6938FirstFewDaysForBatch/
module add mit/matlab/2016b
```
The matlab command is essentially:
```
load(fullfile('/net/feevault/data0/shared/EmilyShijieShared/BatchResultTest/','LogisticscnmfeBatchVer.mat'))
addpath(genpath(codeDir));
cnmfeBatchVer_ClusterPart
```
First load variables made by cnmfeBatchVer_LocalPart. Then add code directory. Run cnmfeBatchVer_ClusterPart.

### (2) cnmfeBatchVer_ClusterPart

From now on, users do not need to look into the code. Just let the clsuter run. Basic ideas in this black box are shown in Figure2 above. Users who would like to improve the method can look into key functions in the figure. They all have documentations inside the functions and can be accessed with "help".

#### Output description
The output is in CNMFE_BatchVer.mat, in which neuron_batch is the structure of the result. See Key Variable Table. Other outputs are meant for parameter choice or for parameter monitor in case somthing goes wrong. See Parameter choice (PNR) for more of this.

## Summary of differences between CNMF-E (basic) and CNMF-E (BatchVer)

Compared to normal CNMF-E with many rounds of iteration of background, spatial and temporal update and merging neurons for each file, BatchVer splits the process on different scales. For sampling, normal rounds of iteration of spatial and temporal update apply. From them, we get a comprehensive and relatively decent A. Then, for all files, use this A to find C on roughly background subtracted data. After one round of temporal update (de-noise and deconvolution), output the result. The output that will be primarily used is A*C. The design of this version sacrifices some precision on the extracted result but enable users to compare activity of neurons over days.

## Parameter choice

There are many options in CNMF-E that can be used in your data's favour. Here we provide some guidelines on gSig, gSiz, PNR.

### gSig and gSiz
* **Amount**
gSiz/gSig ratio is critical. Only some range is good. For example, some range may be 1.2-2. Below 1.2, too few neurons selected, above 2, too messy. There is also one ratio, around which no matter what absolute gSiz and gSig is, results look similar in every aspect.

* **Energy/spatial spread**
(1)Sometimes neurons extracted are with "ripples of tails". Some amounts of tails give us dendrite/axon information, yet too much of this is suspicious in that noise contamination is more likely. 

(2) Sometimes neurons have "energies" "spread out" rather than focused. This can happen when neurons are out of focus or are on the edge of the FOV. Those neurons often are bigger than most of the central neurons. 
The common point of (1) (2) is how much spatial footprint spreadout.

Two factors can affect what level of spread out spatial footprint show up in the final extraction: absolute value of gSiz and gSiz/gSig ratio.
Bigger gSiz => more tails.
Bigger gSiz/gSig ratio => more tails.

In general: gSiz < some level is good.
gSiz/gSig ratio <= some level is good.

If you have a lot of data to operate on, you may end up choosing one conservative parameter for all of them.

### PNR
Compared to min_corr, min_PNR seems to have a much bigger effect. For each sample with initiation involved, a plot is drawn to show the min_corr, min_PNR threshold and the initiated seeds. It is carried out by the newly intoduced method drawPNRCn(min_pnr,min_corr). You can change PNR and see how many seeds get involved. The number on the plot next to each seed correspond to the final neuron number(column/row number) in the neuron.A and neuron.C. You will need to couple this figure with the output of viewNeuron for each sample, each in a folder (spatial and temporal trace for each neuron).

### About RAM and time choice
RAM need is roughly in linear relationship with the amount of sample data you choose. Time choice depends on the amount of both sample data and all data for extraction.

## Notes

This version keeps some features of basic version intact, one of which is the option in initiating in patches. PNR/Correlation diagnose plot is not supported in patch option though interested users are encouraged to modify it. Since we are already operating on small data sets (in sampling stage) and that patching field of view is used to reduce memory use, patching is not necessary here. In addition, there might be some artifacts with the patch option. In short, patch view is not recommended.

## Authors
Pengcheng Zhou initiated the idea. Emily Mackevicius and Shijie Gu developed the idea. Shijie Gu implemented this version.

* **Shijie Gu** [email](techel@live.cn)
* **Emily Mackevicius** [email](elm@mit.edu)
* **Pengcheng Zhou**

## Acknowledgments
Galen Lynch provides two functions that are necessary for smooth parfor function on cluster.
* **Galen Lynch** [email](glynch@mit.edu)

