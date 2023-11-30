# Interaction-analysis
## Lasso
### Input file:
1. Predictor (x) file: a csv file, samples are rows and variables are columns.
2. Independent (y) file: a csv file, samples are rows and variables are columns.

### Script usage:
``` R
./lasso.R \
-s . \ # Path to the script
--input_predictor_file /INPUT/X/FILE/NAME \ # Name of  the input predictor(x) file (csv file,column:taxa name, row:sample)
--input_independent_file /INPUT/Y/FILE/NAME \ # Name of the input independent(y) file (csv file,column:gene name, row:sample)
-n 6 \ # Number of threads
-o  /OUTPUT/PATH \ # Path to the output
-x taxa \ # Name of predictor variable
-y gene # Name of independent variable
```
```
./lasso.R --help # check for usage information of the script.
```

### Output files:
1. output_path/Rawfile

each file represents for a single independent variable with all the predictors.

| column name | explaination|
|-------------|-------------|
| independent variable name | the parameter specified as **-x** |
| predictor variable name | the parameter specified as **-y** |
| r.sqr | R-squared of the *smallest mean squared error* model giving by lambda.min during the **cross validation process**|
| r.sqr.adj | Adjusted R-squared of the *smallest mean squared error* model giving by lambda.min during the cross validation process |
| pval | p value for each predictor varible of the model after desparse lasso estimation analysis |
| padj | Adjusted p value for each predictor varibale |
| ci.lower | The lower confidence interval of coeffient of each predictor variable |
| ci.upper | The upper confidence interval of coeffient of each predictor variable |
| sigma | The standard error of the model on the **whole dataset** instead of cross validation |
| sigma.flag | 1 if the number of non-zero coeffient predictor variables is too big (more than n-2 , n is the sample number) |
| selected | yes if the predictor variable is stably selected by stablity selection |

2. Name_lasso_results_FDR_0.05_stable.csv

combine all the independent variable files and adjusted for the p value, only those with FDR < 0.05 and stably choosed rows are retained.

## SparseCCA
### Input files:
1. Predictor (x) file: a csv file, samples are rows and variables are columns.
2. Independent (y) file: a csv file, samples are rows and variables are columns.
   
### Script usage:
``` R
./sparseCCA.R \
--script_path /SCRIPTS/PATH \
--input_predictor_file /INPUT/X/FILE/NAME \
--input_independent_file /INPUT/Y/FILE/NAME \
-n 8 \
-o /OUTPUT/PATH \
--name Control 

```
### Output files:
1. Name/xx_xx.txt
2. Name/components/independent_predictor_component_x.txt



