# RNAdetector: a user-friendly and stand-alone pipeline for RNA-Seq data analysis
<img src="logo/unict-logo-png-3.png" width="200" hight="100">

***RNAdetector*** is a user-friendly pipeline for the analysis of protein-coding genes and ncRNAs from RNA-Seq data
## Performed analyses
***RNAdetector*** allows to perform several types of analysis such as:
- Quantification and normalization
- Differential expression analysis
- miRNA-sensitive topological pathway analysis.

***RNAdetector*** supports several species such as Human, Mouse and Rat that are available for download in our remote repository.
However, It can be also easily used with any other organism by following the step-by-step procedures detailed in the user interface.
## Non-coding RNAs analyzed
In addition to mRNAs, ***RNAdetector*** can also analyze several classes of small and long ncRNAs.
##### Small non-coding RNAs
- micro RNAs (miRNAs)
- PIWI-associated RNAs (piRNAs)
- tRNA fragments (tRFs)
- tRNA derived small ncRNAs (tsRNAs)
##### Long non-coding RNAs
- long non-coding RNAs (lncRNAs)
- transcribed UltraConserved Regions (tUCRs)
- circular RNAs (circRNAs)

## Features
***RNAdetector*** has several important features such as:
- Easy installation and dependencies management
- Cross-platform (Windows Professional, macOS, Ubuntu) 
- Remotely controllable
- Completely off line
- Graphical User Interface (GUI)
- Graphical reports with the results of the analysis
# RNAdetector installation 
## Installation on Windows Professional
- Install **Docker Desktop** in your computer by downloading the installer at this link [Docker Desktop](https://hub.docker.com/editions/community/docker-ce-desktop-windows/)
- Once the installer is downloaded, then double-click **Docker Desktop Installer.exe** to run the installer.
- When prompted, ensure the Enable Hyper-V Windows Features option is selected on the Configuration page.
- Follow the instructions on the installation wizard to authorize the installer and proceed with the install.
- When the installation is successful, click *Close* to complete **Docker** installation process.
- Now that **Docker Desktop** is installed, it is possible to proceed with **RNAdetector** installation by downloading its installer [RNAdetector installer]() and following the instructions on the installation wizard to authorize the installer and proceed with the installation.
- **RNAdetector** is installed! Now, you can find **RNAdetector** in your application list.

## Installation on macOS
- Before proceed with **Docker Desktop** installation, check the system requirements at this [link](https://docs.docker.com/docker-for-mac/install/#system-requirements)
- Install **Docker Desktop** in your computer by downloading the installer at this link [Docker Desktop](https://hub.docker.com/editions/community/docker-ce-desktop-mac/)
- Once the installer is downloaded, then double-click **Docker.dmg** to open the installer, and drag the Docker icon to the Applications folder.
- Now that **Docker Desktop** is installed, it is possible to proceed with **RNAdetector** installation by downloading its installer [RNAdetector installer]() and following the instructions on the installation wizard to authorize the installer and proceed with the installation.
- Drag the **RNAdetector icon** to the Applications folder to complete the installation.
- **RNAdetector** is installed! Now, you can find **RNAdetector** in your application list.

## Installation on Ubuntu
- Open terminal.

- Install **Docker** in your computer by writing the following commands in your terminal one by one.
```
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
```

- Set your user name in Docker. If you do not now your user name write in your terminal the following command `whoami`.
```
sudo usermod -aG docker YOUR_USER_NAME
```

- Log out and log in from Ubuntu.

- Download the *RNAdetector.deb* file from [RNAdetector web site]().

- Open again the terminal and install **RNAdetector** by writing the following command.
```
sudo dpkg -i RNAdetector.deb
```

- **RNAdetector** is installed! Now, you can find **RNAdetector** in your application list.

# RNAdetector user guide
RNAdetector allows users to perform several RNA-Seq data analysis by using our intuitive GUI. Users can select from the dashboard:

- **Small RNA-Seq analysis** for the identification and quantification of ***miRNAs***, ***piRNAs***, ***tRFs***, and ***tsRNAs***.
- **Long RNA-Seq analysis** for the identification and quatification of ***mRNAs***, ***lncRNAs***, and ***tUCR***.
- **Circular RNA analysis** for identification and quantification of ***circRNAs***.

Once one of the abovementioned analysis is performed, it is possibile to execute downstream analysis such as:

- **Differential expression analysis** to compare difference in the expression profile of such RNA molecules between case vs control samples
- **Pathway analysis** to execute a  miRNA-sensitive topological pathway analysis on the results obtained by the differential expression analysis. However, only the differentially expressed ***mRNAs*** and ***miRNAs*** can be analyzed.

Before proceeding with the analysis, users can download reference indexed genomes available in our remote repository or upload additional indexed genomes from users' computers by using the **Reference Sequences** section of our dashboard and following step-by-step the indication showed on the user interface. Additionally, *GTF* or *BED* genome annotation files can also be uploaded from the **Annotations** section of our dashboard. 

Here follow a description about how to perform each of the above mentioned analysis. 

## Small RNA-Seq analysis
This analysis allows the identification and quantification of ***miRNAs***, ***piRNAs***, ***tRFs***, and ***tsRNAs***. To start the analysis, click **Run Analysis** on the dashboard and then **SmallRNA Analysis**. After that, follow the indications step-by-step described on the user interface to set and start the analysis. Four steps are required
1. **Choose type**. Here users have to indicate the *sample Code*, the *Analysis Name*, the *Input Type* (FASTQ, BAM,SAM), *sequencing strategy* (paired-end or single-end), and the *Number of threads* to be used for the analysis. *Sample code* is used for further analyses such as the Differential Expression Analysis. Specifically, it identifies the sample during the analysis, therefore, it should be a string without any spaces (only letters, numbers, and dashes). Click **Next** to proceed with the next step.
2. **Set pipeline preferences**. Accordingly with the type of the input files, here users can choose which steps will be included in the analysis such as *trimming*, *BAM/SAM to FASTQ conversion*, *alignment*, and *quantification*. If trimming is enabled, it is possible to indicate the minimum PHRED quality and the minimum reads length for the reads quality filtering. For the alignment step, users can choose among ***Salmon***, ***TopHat***, and ***HISAT 2***, while for the quantification ***Feature Counts***, ***HT-seq***, and ***Salmon*** are available. Once everything has been selected, click **Next** to proceed with the next step.
3. **Select references**. If ***TopHat*** or ***HISAT 2*** have been selected for the alignment, at this step users have to select the reference indexed genome and genome annotation file while if ***Salmon*** has been selected, the reference indexed transcriptome must be selected. In addition to the reference genome and transcriptome already available in RNAdetector, it is also possible to add additional indexed genomes, transcriptome and genome annotation files by downloading them from our repository or uploading them from users' computer (read the above section). Once everything has been selected, click **Next** to proceed with the next step.
4. **Upload file**. Here users can add samples to the analysis and upload their files. For each sample, users will be also able to input a custom *Sample Code*. If users are uploading multiple samples for a batch analysis, they can also upload a sample description file (in TSV format) that can be used for the differential expression analysis. Once the samples have been uploaded, click on the **Start Analysis** button to start the analysis.
## Long RNA-Seq analysis
This analysis allows the identification and quantification of ***mRNAs***, ***lncRNAs***, and ***tUCR***. To start the analysis, click **Run Analysis** on the dashboard and then **LongRNA Analysis**. After that, follow the indications step-by-step described on the user interface to set and start the analysis. Four steps are required
1. **Choose type**. Here users have to indicate the *sample Code*, the *Analysis Name*, the *Input Type* (FASTQ, BAM,SAM), *sequencing strategy* (paired-end or single-end), and the *Number of threads* to be used for the analysis. *Sample code* is used for further analyses such as the Differential Expression Analysis. Specifically, It identifies the sample during the analysis, therefore, it should be a string without any spaces (only letters, numbers, and dashes). Click **Next** to proceed with the next step.
2. **Set pipeline preferences**. Accordingly with the type of the input files, here users can choose which steps will be included in the analysis such as *trimming*, *BAM/SAM to FASTQ conversion*, *alignment*, and *quantification*. If trimming is enabled, it is possible to indicate the minimum PHRED quality and the minimum reads length for the reads quality filtering. For the alignment step users can choose among ***Salmon***, ***TopHat***, and ***HISAT 2***, while for the quantification ***Feature Counts***, ***HT-seq***, and ***Salmon*** are available. Once everything has been selected, click **Next** to proceed with the next step.
3. **Select references**. If ***TopHat*** or ***HISAT 2*** have been selected for the alignment, at this step users have to select the reference indexed genome and genome annotation file while if ***Salmon*** has been selected, the reference indexed transcriptome must be selected. In addition to the reference genome and transcriptome already available in RNAdetector, it is also possible to add additional indexed genomes, transcriptome and genome annotation files by downloading them from our repository or uploading them from users' computer (read the above section). Once everything has been selected, click **Next** to proceed with the next step.
4. **Upload file**. Here users can add samples to the analysis and upload their files. For each sample, users will be also able to input a custom *Sample Code*. If users are uploading multiple samples for a batch analysis, they can also upload a sample description file (in TSV format) that can be used for the differential expression analysis. Once the samples have been uploaded, click on the **Start Analysis** button to start the analysis.
## Circular RNA analysis
This analysis allows the identification and quantification of ***circRNAs***. To start the analysis, **Run Analysis** on the dashboard and then **CircRNA Analysis**. After that, follow the indications step-by-step described on the user interface to set and start the analysis. Four steps are required
1. **Choose type**. Here users have to indicate the *sample Code*, the *Analysis Name*, the *Input Type* (FASTQ, BAM,SAM),and *sequencing strategy* (paired-end or single-end). *Sample code* is used for further analyses such as the Differential Expression Analysis. It identifies the sample during the analysis, therefore, it should be a string without any spaces (only letters, numbers, and dashes). Click **Next** to proceed with the next step.
2. **Set pipeline preferences**. Accordingly with the type of the input files, here users can choose which steps will be included in the analysis such as *trimming*, *BAM/SAM to FASTQ conversion*, which version of ***Ciri*** should be used for the identification and quantification of circRNAs, and the *Number of threads* to be used for the analysis. If trimming is enabled, it is also possible to indicate the minimum PHRED quality and the minimum reads length for the reads quality filtering. Once everything has been selected, click **Next** to proceed with the next step.
3. **Select references**. Here users have to select the reference indexed genome, the genome annotation file (for circRNAs analysis the genomic annotation of coding-protein genes must be selected), and the BED file with the Back-Spliced Junction Annotation to be used for the analysis. Once everything has been selected, click **Next** to proceed with the next step.
4. **Upload file**. Here users can add samples to the analysis and upload their files. For each sample, users will be also able to input a custom *Sample Code*. If users are uploading multiple samples for a batch analysis, they can also upload a sample description file (in TSV format) that can be used for the differential expression analysis. Once the samples have been uploaded, click on the **Start Analysis** button to start the analysis.
## Creation of a new sample group
## Differential expression analysis
To start the analysis, click **Run Analysis** on the dashboard and then **New DEGs Analysis**. After that, follow the indications step-by-step described on the user interface to set and start the analysis. For the differential expression analysis six steps are required.
1. **Choose name and type**. Here users can write the identification code of the samples, the name of the analysis, and select the previously generated sample group (read the previous section). Click **Next** to proceed with the next step.
2. **Select Variables and Contrasts**. Here users can choose which variables will be used to define the contrasts. After selecting the variables, you will be able to add one or more contrasts and defining case and control samples. Click **Next** to proceed with the next step.
3. **Common Parameters**. Here users can select several parameters which are essential for the analysis such as the p-value threshold, the Log offset, when the filtering should be applied, which method should be used for the p-values adjustment, which method should be used for the calculation of the meta-analysis p-value, which format should be used to export the final figures and the number of threads that have to be used for the analysis. Once everything has been selected, click **Next** to proceed with the next step.
4. **Normalization Parameters**. Here users can select the normalization method (DESeq2 or edgeR) and its parameters accordingly with user preferences. Once everything has been selected, click **Next** to proceed with the next step.
5. **Filtering Parameters**. Here users can choose if use several filtering criteria for the selection of transcripts to be considered for the differential expression analysis. Users can choose to filter transcripts according with their lenght, average number of mapped reads, overal expression, and the prevalence of that transcript across the samples by indicating the percentage of samples which have to express that transcripts and its threshold number of mapped reads in order to be considered expressed. After that all the filtering criteria for the transcript selection have been selected, click **Next** to proceed with the next step.
6. **Statistics Parameters**. Here users can choose among several statistical methods (LIMMA, DESeq2, and edgeR) and their parameters to be used for the differential expression analysis. After entering all the parameters, you can start the analysis by clicking on the **Start Analysis** button.
## Pathway analysis
After the differential expression analysis (read the previous section), it is also possible to perform a pathway analysis. It is important to clarify that pathway analysis can only be perfomed on *mRNAs* and\or *miRNAs* which have been found differentially expressed by the previous analysis. To start the pathway analysis, click **Run Analysis** on the dashboard and then **New Pathway Analysis**. After that, follow the indications step-by-step described on the user interface to set and start the analysis. For pathway analysis three steps are required.
1. **Choose name and type**. Here users can write the sample code, the name of the analysis, and selected the differential expression analysis results to be used for the pathway analysis. Click **Next** to proceed with the next step.
2. **DEGs Analysis Parameters**. Here users can choose the criteria for the selection of the differentially expressed transcripts (mRNAs and\or miRNAs) to be used for the pathway anaysis by writing the p-value cutoff and the Log-Fold change threshold. Click **Next** to proceed with the next step.
3. **Pathway Analysis Parameters**. Here users can choose the organism to be used for the pathway analysis and the p-values threshold for the selection of the statistically significant perturbated pathways. After entering all the parameters, you can start the analysis by clicking on the **Start Analysis** button.
## Final results
# Cite us
if you use ***RNAdetector*** cite:
# Contact us 
[Alessandro La Ferlita](https://www.researchgate.net/profile/Alessandro_La_Ferlita2) (alessandro.laferlita@unict.it)

[Salvatore Alaimo](https://www.researchgate.net/profile/Salvatore_Alaimo) (alaimos@dmi.unict.it)

[Alfredo Pulvirenti](https://www.researchgate.net/profile/Alfredo_Pulvirenti) (alfredo.pulvirenti@unict.it)
