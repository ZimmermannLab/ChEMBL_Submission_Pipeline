from rdkit import Chem
import pandas as pd
from rdkit.Chem import PandasTools
import os
from rdkit import RDConfig
import itertools

def check():
    pass

def generate_compound_files(output_dir, input_file, RIDX, prefix):
    try:
        if "sdf" in input_file:
            sdfFile = os.path.join(RDConfig.RDDataDir, input_file)
            frame = PandasTools.LoadSDF(sdfFile, includeFingerprints=True)
            RIDX = ", ".join(itertools.repeat(RIDX, len(frame))).split(', ')
            data_list = list(range(1,len(frame)+1))
            num_to_use = len(str(len(frame))) + 1
            CIDX = [prefix + str(item).zfill(num_to_use) for item in data_list]
            frame["CIDX"] = CIDX
            frame["RIDX"] = RIDX
            if "COMPOUND_NAME" not in frame.columns:
                frame["COMPOUND_NAME"] = frame["ID"]
            if "COMPOUND_KEY" not in frame.columns:
                frame["COMPOUND_KEY"] = frame["ID"]
            PandasTools.WriteSDF(frame, output_dir + "/COMPOUND_CTAB.sdf", properties=['CIDX'])
            if "COMPOUND_SOURCE" in frame:
                frame_record = frame[["CIDX", "RIDX", "COMPOUND_NAME", "COMPOUND_KEY", "COMPOUND_SOURCE"]]
                frame_record.to_csv(output_dir + "/COMPOUND_RECORD.tsv", sep ="\t")
                return(frame)
            else: 
                frame_record = frame[["CIDX", "RIDX", "COMPOUND_NAME", "COMPOUND_KEY"]]
                frame_record.to_csv(output_dir + "/COMPOUND_RECORD.tsv", sep ="\t")
                return(frame)

        elif "csv" in input_file:
            frame = pd.read_csv(input_file)
            frame['ROMol'] = frame['SMILES'].apply(lambda x: Chem.MolFromSmiles(x))
            RIDX = ", ".join(itertools.repeat(RIDX, len(frame))).split(', ')
            data_list = list(range(1,len(frame)+1))
            num_to_use = len(str(len(frame))) + 1
            CIDX = [prefix + str(item).zfill(num_to_use) for item in data_list]
            frame["CIDX"] = CIDX
            frame["RIDX"] = RIDX
            if "COMPOUND_NAME" not in frame.columns:
                frame["COMPOUND_NAME"] = frame["ID"]
            if "COMPOUND_KEY" not in frame.columns:
                frame["COMPOUND_KEY"] = frame["ID"]
            PandasTools.WriteSDF(frame, output_dir + "/COMPOUND_CTAB.sdf", properties=['CIDX'])
            if "COMPOUND_SOURCE" in frame:
                frame_record = frame[["CIDX", "RIDX", "COMPOUND_NAME", "COMPOUND_KEY", "COMPOUND_SOURCE"]]
                frame_record.to_csv(output_dir + "/COMPOUND_RECORD.tsv", sep ="\t")
                return(frame)
            else: 
                frame_record = frame[["CIDX", "RIDX", "COMPOUND_NAME", "COMPOUND_KEY"]]
                frame_record.to_csv(output_dir + "/COMPOUND_RECORD.tsv", sep ="\t")
                return(frame)
    except:
        print("input_file should either be a csv or sdf as file, hence name should contain '.csv' or '.sdf'")

