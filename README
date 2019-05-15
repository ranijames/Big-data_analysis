General dependencies
General needs for the tools to work:
python3 installed. Best option is to have conda version installed.
Most important dependencies:
pandas: conda install pandas
numpy: conda install numpy
Installation of this repository
Install this repository (in develop mode) as follows:
pip install -e .
After the installation you will have a series of comands available, all with a prefix according to their role:
rna-trimandalign
exon-bcbio
exon-copywriter
exon-sciclone
For more commands just try with rna- and exon- and hit tab in the terminal.
Configuration file
A configuration is needed at the exact location: ~/.agbaldus/agbaldus-seq.conf. 
You may copy the template file found in this repository at allpy/config/agbaldus-seq.conf.template and adapt the configuration.
Running the tools
Each tool has a help which can be printed to the output with -h.
RNA-Seq
Make sure the STAR genome is available
Exon-Seq
Installation of bcbio-nextgen
If you're sitting behind a proxy, make sure you connect to GitHub via the https protocol:
git config --global url.https://github.com/.insteadOf git://github.com/
Download the installer script of bcbio_nextgen
wget https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
Define a directory in your home dir where you want to install bcbio_nextgen and install:
BCB_INSTALL_DIR=~/bin/bcbionextgen
python2 bcbio_nextgen_install.py $BCB_INSTALL_DIR/data/ --tooldir $BCB_INSTALL_DIR/tools/ --genomes GRCh37 --aligners star --aligners bwa --distribution ubuntu
Add data rich supplemental software (not installed by default)
bcbio_nextgen.py upgrade --tools --toolplus data
Download and extract GATK and then link it to bcbio-nextgen
GATK_DIR=~/bin/GATK
mkdir $GATK_DIR
bunzip2 -dc GenomeAnalysisTK-3.4-0.tar.bz2 | tar xvf - -C $GATK_DIR
bcbio_nextgen.py upgrade --tools --toolplus gatk=$GATK_DIR/GenomeAnalysisTK.jar

