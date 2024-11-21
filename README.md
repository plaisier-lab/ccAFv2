## ccAFv2 code repository
This github repository contains all the codes used to construct the ccAFv2 classifier. The scripts used to generate the code can be found in the 'scripts' directory.

# Building the ccAFv2 classifier
The construction of the ccAFv2 classifier takes place in the file:  [04_building_ccAFv2.py](https://github.com/plaisier-lab/ccAFv2/blob/main/code/04_building_ccAFv2.py)
This script contains the ANN model for ccAFv2, which can be repurposed to build classifiers to predict cell types in other classifiers by exchanging the U5-hNSC scRNA-seq data and labels with new data and labels. The classifier_ccAFv2 class can take in any dataset, and provides parameters to adjust the construction of the classifier. It does this by taking those parameters in through it's init funciton:

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

- **training_sets**: a dictionary of training scRNA-seq Scanpy objects, no default
- **label_sets**: a dictionary of training scRNA-seq cell label lists, no default
- **variable_genes**: number of genes to include in classifier, default is 1000 (after filtering may be fewer)
- **epochs**: number of training epochs, default is 10
- **validation_split**: what fraction of cells should be set aside for validation, defaults to 0.2
- **activation**: activation function type for artificial neurons, defaults to 'relu'
- **dropout_rate**: rate to use for dropout nodes between layers to avoid overfitting, defaults to 0.4
- **model_specs**: size of the layers of the network, defaults to [200,100] (for ccAFv2 we used [600,200])
- **model_name**: adds the calendar time to the outputs for the classifier, default date and time

The script details how to import the data, initialize, and construct a classifier. We also offer all the other scripts needed to test the classifier and demonstrate its utility.
