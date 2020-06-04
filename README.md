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
- Graphical reports with the results of the analyses
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
RNAdetector allows users to perform several RNA-Seq data analysis by using our intuitive GUI. Users can select:

- **Small RNA-Seq analysis** for the identification and quantification of ***miRNAs***, ***piRNAs***, ***tRFs***, and ***tsRNAs***.
- **Long RNA-Seq analysis** for the identification and quatification of ***mRNAs***, ***lncRNAs***, and ***tUCR***.
- **Circular RNA analysis** for identification and quantification of ***circRNAs***.

Once one of the abovementioned analysis is performed, it is possibile to execute downstream analysis such as:

- **Differential expression analysis** to compare difference in the expression profile of such RNA molecules between case vs control samples
- **Pathway analysis** to execute a  miRNA-sensitive topological pathway analysis on the results obtained by the differential expression analysis. However, only the differentially expressed ***mRNAs*** and ***miRNAs*** can be analyzed.

Here follow a description of how to perform each of the mentioned analysis. 

## Small RNA-Seq analysis
## Long RNA-Seq analysis
## Circular RNA analysis
## Differential expression analysis
## Pathway analysis
# Cite us
if you use ***RNAdetector*** cite:
# Contact us 
[Alessandro La Ferlita](https://www.researchgate.net/profile/Alessandro_La_Ferlita2) (alessandro.laferlita@unict.it)

[Salvatore Alaimo](https://www.researchgate.net/profile/Salvatore_Alaimo) (alaimos@dmi.unict.it)

[Alfredo Pulvirenti](https://www.researchgate.net/profile/Alfredo_Pulvirenti) (alfredo.pulvirenti@unict.it)
