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
