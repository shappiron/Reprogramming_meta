import pandas as pd
import json
from copy import deepcopy

# save structure
def save_structure(ultradf, path='test.json'):
    ultrasave = deepcopy(ultradf)
    for e in ultrasave.keys():
        ultrasave[e]['data'] = ultrasave[e]['data'].reset_index().to_csv()
        ultrasave[e]['pheno'] = ultrasave[e]['pheno'].to_csv()
    with open(path, 'w+') as f:
        json.dump(ultrasave, f)
        
# load structure
def load_structure(path='ultradf.json'):
    import sys
    if sys.version_info[0] < 3: 
        from StringIO import StringIO #convert string to IO monade
    else:
        from io import StringIO
    #load after all
    with open(path, 'r') as f:
        ultradf = json.load(f)
        
    for e in ultradf.keys():
        ultradf[e]['data'] = pd.read_csv(StringIO(ultradf[e]['data']), index_col=0).set_index('index')
        ultradf[e]['pheno'] = pd.read_csv(StringIO(ultradf[e]['pheno']), index_col=0)
    
    return ultradf