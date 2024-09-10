## ccAFv2 code repository
This github repository contains all the codes used to construct the ccAFv2 classifier. The scripts used to generate the code can be found in the 'scripts' directory.

# Building the ccAFv2 classifier
The construction of the ccAFv2 classifier takes place in the file:  [04_building_ccAFv2.py](https://github.com/plaisier-lab/ccAFv2/blob/main/scripts/04_building_ccAFv2.py)
This script contains the ANN model for ccAFv2, which can be repurposed to build classifiers to predict cell types in other classifiers by exchanging the U5-hNSC scRNA-seq data and labels with new data and labels. The classifier_ccAFv2 class can take in any dataset, and provides parameters to adjust the construction of the classifier:

```python
def __init__(self,
             training_sets,
             label_sets,
             variable_genes = 1000,
             epochs = 10,
             training_iterations = 10,
             validation_split = 0.2,
             activation = 'relu',
             dropout_rate = 0.4,
             model_specs = [200,100],
             model_name = calendar.timegm(time.gmtime()))
```
