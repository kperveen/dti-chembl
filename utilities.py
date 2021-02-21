import pandas as pd
import numpy as np

class Util:
    def __init__(self, data, target_ids):
        self.mol_data = data
        self.target_abbr = {
            'CHEMBL1910': 'ADA',
            'CHEMBL203': 'EGFR-erbB1',
            'CHEMBL2189121': 'GTPase KRas',
            'CHEMBL1957': 'IGF-1',
            'CHEMBL614725': 'MIA PaCa-2',
            'CHEMBL3580494': 'MUC-1',
            'CHEMBL3474': 'sPLA2-IIA',
            'CHEMBL2842': 'mTOR',
            'CHEMBL4026': 'STAT-3',
            'CHEMBL209': 'Trypsin I',
            'CHEMBL1955': 'VEGFR3'
        }
        self.targets_dict = {}
        self.target_ids = target_ids

    def get_target_abbr(self, target_id):
        return self.target_abbr[target_id]

    def set_target_dict(self):
        for target in self.target_ids:
            self.targets_dict[target] = {'percentiles': [], 25: 0, 50: 0, 75: 0, 100: 0}
        return self.targets_dict

    def get_standard_values(self, target_id):
        svalues = self.mol_data[self.mol_data["Target ChEMBL ID"] == target_id]["Standard Value"].tolist()
        return svalues

    def target_name_from_id(self, target_id):
        name = self.mol_data[self.mol_data["Target ChEMBL ID"] == target_id]["Target Name"].tolist()[0]
        return name

    def target_abbr_from_id(self, target_id):
        #name = self.mol_data[self.mol_data["Target ChEMBL ID"] == target_id]["Target Name"].tolist()[0]
        abbr = self.target_abbr[target_id]
        return abbr

    def set_percentile_count(self):
        for i in self.target_ids:
            self.percentile_count(i)
        return self.targets_dict

    def percentile_count(self, target_id):
        svalues = self.get_standard_values(target_id)
        self.targets_dict[target_id]["percentiles"].append(min(svalues))
        val_25 = np.percentile(svalues, 25)
        self.targets_dict[target_id]["percentiles"].append(val_25)
        self.targets_dict[target_id][25] = len([i for i in svalues if i <= val_25])
        val_50 = np.percentile(svalues, 50)
        self.targets_dict[target_id]["percentiles"].append(val_50)
        self.targets_dict[target_id][50] = len([i for i in svalues if (i > val_25 and i <= val_50)])
        val_75 = np.percentile(svalues, 75)
        self.targets_dict[target_id]["percentiles"].append(val_75)
        self.targets_dict[target_id][75] = len([i for i in svalues if (i > val_50 and i <= val_75)])
        val_100 = max(svalues)
        self.targets_dict[target_id]["percentiles"].append(val_100)
        self.targets_dict[target_id][100] = len([i for i in svalues if (i > val_75 and i <= val_100)])

