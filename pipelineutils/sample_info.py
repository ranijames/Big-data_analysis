import os
import pandas as pd

dtypes = {'purity': int, 'cna_log2_center_correction': float, 'gender': str}
__author__ = 'ARJ'


def get_sample_info(patient_dir):
    info_file = os.path.join(patient_dir, 'sample_info.txt')
    if os.path.exists(info_file):
        df = pd.read_table(info_file, index_col=0).T
        df = df.astype({k:v for k,v in dtypes.items() if k in df.index})
        return df
    else:
        return pd.DataFrame(columns=['purity','gender'])
