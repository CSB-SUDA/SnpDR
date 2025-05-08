# SnpDR

## Framework
![Framework](Framework.png)

## How to use the script?
### Requirements
* R
* Anaconda3

These codes run in R, but some of their functions depend on a Python environment. Therefore, it is need to create a virtual environment named `pyUser` using anaconda3, with the Python version no lower than 3.7. Please open the Anaconda Prompt and enter the following command:
```
conda create -n pyUser python=3.9
conda activate pyUser
pip install rdkit
pip install descriptastorus 
pip install DeepPurpose
pip install seaborn
pip install goatools
pip install prody
conda deactivate
```
The Python environment is created to use the DeepPurpose tool and PRS analysis. For more details about DeepPurpose, please visit https://github.com/kexinhuang12345/DeepPurpose.

> **Note that the path to the pyUser environment should be obtained using `conda env list`. In subsequent analyses, this path will be used as a parameter in some functions.**

### Usage

#### Step1: Prepare
Download this repository to your local machine and extract it. Then, in R, set the working directory to the extracted folder and load the required packages:
```
setwd("this_file_path/SnpDR")
source("Rscript/0.Load_packages.R")
```

#### Step2: Signatrue Module
Preprocess the expression profile data.
```
source("Rscript/1.Signatrue.R")
expr_process(df_file="path/expression.txt")
```
* _expr_file_: The path of TXT file storing expression profile. 
  
  |protein| sample1 | sample2 |... |
  | --- | --- | --- | --- |
  | protein1 | 2.345 | 6.480 | ... |
  | protein2 | 7.985 | 4.621  | ... |
  | ... | ... | ...  |...|
