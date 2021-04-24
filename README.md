# RNAdetector: a free user-friendly stand-alone and cloud-based system for RNA-Seq data analysis
<img src="logo/logoRNAdetectorBlack.png" width="300" hight="200">

Developed by

<img src="logo/new_unict_logo.png" width="200" hight="100">

in collaboration with

<img src="logo/osu_logo.jpg" width="200" hight="100"><img src="logo/nms_logo.png" width="200" hight="100"><img src="logo/ior_logo.jpg" width="200" hight="100">

***RNAdetector*** is a user-friendly software for the analysis of protein-coding genes and ncRNAs from RNA-Seq data. It can be run as stand-alone application or as a cloud based system.
## Performed analyses
***RNAdetector*** allows to perform several types of analysis such as:
- Quantification and normalization
- Differential expression analysis
- miRNA-sensitive topological pathway analysis.

***RNAdetector*** supports any organisms whose genomes\transcriptoms have been sequenced. However, some commonly studied species such as ***Human***, ***Mouse***, and ***C.elegans*** are available for download in our remote repository. The list of organisms will keeped updated.

***Any other organisms*** can be easily added by uploading their genomes and\or transcriptomes and the genomic coordinates of the RNA molecules intended to be analyzed following the step-by-step procedures detailed in the user interface.
## Non-coding RNAs analyzed
In addition to mRNAs, ***RNAdetector*** can also analyze several classes of small and long ncRNAs. Specifically, for Human, Mouse and C.elegans RNA-Seq data the following ncRNA classes can be analyzed:
##### Small non-coding RNAs
- micro RNAs (miRNAs)
- PIWI-associated RNAs (piRNAs)
- Small nucleolar RNAs (snoRNAs)
- tRNA derived small ncRNAs (tRFs and tsRNAs)
##### Long non-coding RNAs
- long non-coding RNAs (lncRNAs)
- transcribed UltraConserved Regions (tUCRs) (only for human)
- circular RNAs (circRNAs)

The genomic annotations of all these above-mentioned classes of ncRNAs are already available for downloading from our remote repository. However, additional ncRNAs classes for human and other organism can be analyzed by uploading their genomic annotation or indexed transcriptome on ***RNAdetector*** following the step-by-step procedures detailed in the user interface. Therefore, any class of ncRNAs of any biological species can be analyzed by ***RNAdetector***.

## Features
***RNAdetector*** has several features such as:
- Easy installation and dependencies management
- Cross-platform (Windows Professional, macOS, Ubuntu) and deployable as ***cloud-based system***
- Remotely controllable <!-- Completely off line -->
- Graphical User Interface (GUI)
- Interactive graphical reports with the results of the analysis
# RNAdetector installation
**Detailed information concerning the [system requirements and installation](https://github.com/alessandrolaferlita/RNAdetector/wiki/Requirements-and-Setup) can be found in the wiki section**
## Installation on Windows Professional
- Install **Docker Desktop** in your computer by downloading the installer at this link [Docker Desktop](https://hub.docker.com/editions/community/docker-ce-desktop-windows/)
- Once the installer has been downloaded, then double-click **Docker Desktop Installer.exe** to run it.
- When prompted, select "Enable Hyper-V Windows Features" option in the Configuration page.
- Follow the instructions on the installation wizard to authorize the installer and proceed with the installation.
- When the installation has been done successfully, click *Close* to complete **Docker** installation process.
- Now that **Docker Desktop** is installed, it is possible to proceed with **RNAdetector** installation by downloading its installer [RNAdetector installer](https://rnadetector.atlas.dmi.unict.it/download.html) and following the instructions through the installation wizard.
- **RNAdetector** is installed! Now, you will find **RNAdetector** in your application list.

## Installation on macOS
- Before proceed with **Docker Desktop** installation, check the system requirements at this [link](https://docs.docker.com/docker-for-mac/install/#system-requirements)
- Install **Docker Desktop** in your computer by downloading the installer at this link [Docker Desktop](https://hub.docker.com/editions/community/docker-ce-desktop-mac/)
- Once the installer has been downloaded, then run the **Docker.dmg** installer, and drag and drop the Docker icon the Applications folder.
- Now that **Docker Desktop** is installed, it is possible to proceed with **RNAdetector** installation by downloading its installer [RNAdetector installer](https://rnadetector.atlas.dmi.unict.it/download.html) and following the instructions through the installation wizard.
- Drag and drop the **RNAdetector icon** to the Applications folder to complete the installation.
- **RNAdetector** is installed! Now, you will find **RNAdetector** in your application list.

## Installation on Ubuntu
- Open terminal.

- Install **Docker** in your computer by writing the following commands in your terminal (one by one).
```
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
```

- Set your user name in Docker. If you do not know your user name, write the following command `whoami` in your terminal.
```
sudo usermod -aG docker YOUR_USER_NAME
```

- Log out and log in from Ubuntu.

- Download the *RNAdetector.deb* file from [RNAdetector web site](https://rnadetector.atlas.dmi.unict.it/download.html).

- Open again the terminal and install **RNAdetector** by writing the following command.
```
sudo dpkg -i RNAdetector.deb
```

- **RNAdetector** is installed! Now, you can find **RNAdetector** in your application list.

# RNAdetector user guide
RNAdetector allows users to perform several RNA-Seq data analysis by using our intuitive GUI. Users can select from the dashboard:

- **Small RNA-Seq analysis** for the analysis of small RNA-Seq data. It can be used for identification and quantification of small ncRNAs such as ***miRNAs***, ***piRNAs***, ***tRFs***, and ***tsRNAs***.
- **RNA-Seq analysis** for the analysis of standard RNA-Seq data. It can be used for the identification and quatification of RNA species longer than 200 nt such as ***mRNAs***, ***lncRNAs***, and ***tUCR***.
- **Circular RNA analysis** for identification and quantification of ***circRNAs***.

Once one of the abovementioned analysis is performed, it is possibile to execute downstream analysis such as:

- **Differential expression analysis** to compare difference in the expression profile of such RNA molecules between case vs control samples;
- **Pathway analysis** to execute a  miRNA-sensitive topological pathway analysis on the results obtained by the differential expression analysis. However, currently only the differentially expressed ***mRNAs*** and ***miRNAs*** can be analyzed.

**A more detailed [user guide](https://github.com/alessandrolaferlita/RNAdetector/wiki) is available in the wiki section**

# Cite us
if you use ***RNAdetector*** cite:
# Contact us 
[Alessandro La Ferlita](https://www.researchgate.net/profile/Alessandro_La_Ferlita2) (alessandro.laferlita@unict.it)

[Salvatore Alaimo](https://www.researchgate.net/profile/Salvatore_Alaimo) (alaimos@dmi.unict.it)

[Alfredo Pulvirenti](https://www.researchgate.net/profile/Alfredo_Pulvirenti) (alfredo.pulvirenti@unict.it)
