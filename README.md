# FastQC Simulator

FastQC is a program designed to do the quality control on raw sequence data coming from high throughput sequencers. The program evaluates the quality of reads in different analysis modules, providing various statistics plots about the quality of reads. This program replicates the original [FastQC](https://github.com/s-andrews/FastQC/tree/fdc27fe360af44dc675e7d67bc3ec89c21f273e6) program written in Java and released in 2010, making it available in Python. 

## Table of Contents
1. [Team: Fasta and Curious](#Team)
2. [Installation and Usage](#instaus)
3. [Test data](#test_data)
4. [Software Requirements](#Software)

<a name="Team"></a>
## Team: Fasta and Curious

Team:

1. Anna Toidze:  
 	* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/github.svg height=20>: [AnnaToi01](https://github.com/AnnaToi01)
	* Tasks:
		* Basic Statistics (part of FastQC)
		* Per tile sequence quality (part of FastQC, together with Ivan Semenov)
		* Adapter content (part of FastQC)
		* Sequence length distribution (part of FastQC)
		* README.md (together with Anton Muromtsev)
2. Anton Muromtsev:
	* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/github.svg height=20>: [AntonMuromtsev](https://github.com/AntonMuromtsev)
	* Taks:
      * Overrepresented sequences (Part of FastQC)
      * Sequence Duplication Levels (Part of FastQC)
      * Binding images into .pdf (service function)
      * README.me (together with Anna Toidze)
3. Ivan Semenov 
	* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/github.svg height=20>: [ipsemenov](https://github.com/ipsemenov)
	* Tasks:
      * Per tile sequence quality (part of FastQC, together with Anna Toidze)
      * Quality scores across all bases (Part of FastQC)
      * Quality score distribution over all sequences (Part of FastQC)
4. Mikhail Fofanov
	* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/github.svg height=20>: [MVFofanov](https://github.com/MVFofanov)
	* Tasks:
        * Sequence content across all bases (Part of FastQC)
        * N content across all bases (Part of FastQC)
        * GC distribution over all sequences (Part of FastQC)

<a name="instaus"></a>
## Installation and Usage

### Pipeline structure
Pipeline is separated into several files located in `scripts` folder:

* `caclulations.py` contains functions required for calculations

* `plotting.py` contains functions for plotting graphs

* `main.py` contains main piece pf code for running all computations




### Preliminary settings
Clone repository
```
$ git clone git@github.com:ipsemenov/FastQC_simulator.git
```
Move to project directory 
```
$ cd FastQC
```
Set up virtual environment in working directory:
1. Create virtual environment
    * Via `virtualenv`

       * Install virtualenv if it is not installed.
         ```
         $ pip install virtualenv
         ```
       * Create virtual environment
         ```
         $ virtualenv venv --python=3.8
         ```
       * Activate it
         ```
         $ source ./venv/bin/activate
         ```
    * Via `conda`
        * [Install Anaconda](https://docs.anaconda.com/anaconda/install/index.html)
        * Create virtual environment
           ```
           $ conda create --name <env_name> python=3.8
           ```
        * Activate it
           ```
           $ conda activate <env_name>
2. Install necessary libraries
 ```
$ pip install -r requirements.txt
 ```
3. Install package wkhtmltopdf
```
$ sudo apt-get install wkhtmltopdf
```

### Console interface
This instrument is a console utility maintaining following parameters:
```
  -i ,  --input     path to fastq file`
  -o , --output     path to output folder for storing results
  -a , --adapters   path to file with adapters. Default: ./adapters.txt
```


### Running utility
Example workflow:
```
$ python main.py -i <path_to_fastq> -o <path_to_ouptut_dir>  -a <path_to_adapters> 
```
To show brief information about parameters, execute following command:
```	
$ python main.py -h
```

<a name="test_data"></a>
## Test data

Test data can be found in `test_data/` folder. The `amp_res_1.fastq.gz` has to be unzipped:
```
$ gunzip amp_res_1.fastq.gz
```
Run the program on test data (from the `scripts` folder):
```
$ python main.py -i ../test_data/amp_res_1.fastq -o ../test_data/results -a ./adapters.txt
```
The test result, `results.pdf`, can be found in 'test_data/results', along with single images of the statistics modules.

<a name="Software"></a>
## Software Requirements
<img src=https://img.shields.io/badge/FastQC%20Simulator-FASTQ%20Quality%20Check-informational height = 20>

* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/python.svg height=20> Python 3.8
* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/ubuntu.svg height = 20> Ubuntu 20.04 and 21.04
* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/git.svg height = 20> Git 2.30.2
* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/markdown.svg height=20> Markdown
* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/html5.svg height=20> HTML
* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/github.svg height=20> GitHub
* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/gnubash.svg height=20> Bash
* Rest of the requirements are in requirements.txt.

