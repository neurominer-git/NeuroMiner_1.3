# function that simulates data with the same structure as input data
import pandas as pd
import numpy as np
from sdv.metadata import SingleTableMetadata
from sdv.single_table import GaussianCopulaSynthesizer
from sdv.constraints import OneHotEncoding, FixedCombinations
from sdv.sampling import Condition

data = pd.read_csv(data_file)
#data['label'] = labels 

metadata = SingleTableMetadata()
metadata.detect_from_dataframe(data)

# preallocate space for simulated dataset
sim_sample = pd.DataFrame(columns = data.columns)

# constraints         
constraints = []

# keep relationship between outer and inner leave-one-group-out categorization the same
if cv1lco > 0 and cv2lco > 0:
    fixed_cvlco_constraint = FixedConstrained(column_names = ['CV1LCO', 'CV2LCO'])
    constraints.append(fixed_cvlco_constraint)

# consider dummy-coded sites    
if type(sitesCols) is not int:
    sitesColsPy = [x-1 for x in sitesCols]
    auxSitesColsNames = data.columns[sitesColsPy]
 
    sites_constraint =  {
        'constraint_class': 'OneHotEncoding',
        'constraint_parameters': {
            'column_names': auxSitesColsNames
        }
    }

    print(sites_constraint)
    constraints.append(sites_constraint)
        
model = GaussianCopulaSynthesizer(metadata)
model.add_constraints(constraints = constraints)    
model.fit(data)

# conditions 
print(len(condVals))
if len(condVals) > 0:
    condition = pd.DataFrame({condColName : condVals})#, num_rows = condN)
    sim_sample = model.sample_remaining_columns(condition)
else: 
    sim_sample = model.sample(n_obs)


out_path = f'{rootdir}/simData.csv'
sim_sample.to_csv(out_path,index=False)
