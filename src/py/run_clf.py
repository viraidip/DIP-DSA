import os
import RNA
import json
import joblib
import numpy as np
import pandas as pd

from Bio.Seq import Seq
from sklearn.preprocessing import OneHotEncoder

CHARS = 'ACGU'
CHARS_COUNT = len(CHARS)
MAX_LEN = 2361

def run_classification(s, e, strain, segment, sequence, clf)-> str:
    features = ["Segment", "DI_Length", "Direct_repeat","Junction", "3_5_ratio", "length_proportion", "delta_G", "Peptide_Length"]

    strain_dict = dict({"A/California/07/2009": "Cal07",
                        "A/NewCaledonia/1999": "NC",
                        "A/Perth/16/2009": "Perth",
                        "A/PuertoRico/8/1934": "PR8",
                        "B/Lee/1940": "BLEE",
                        "A/WSN/33": "WSN"
                        })
    strain = strain_dict[strain]

    # create df
    data = dict({"Start": s, "End": e, "Strain": strain, "Segment": segment, "Sequence": sequence})
    df = pd.DataFrame(data, index=[1])

    # create artificial samples to let segment OHE run correctly
    sgms = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
    n_df = [{"Start": 1,"End":2,"Strain":strain,"Segment":seg,"Sequence": sequence} for seg in sgms]
    n_df = pd.DataFrame(n_df)
    
    df = pd.concat([n_df, df], ignore_index=True)

    # create features
    df, feature_cols = generate_features(df, features)

    # remove artificial samples
    df = df[df["Start"] != 1]
    X = df[feature_cols]

    # load pickle file
    p = os.path.join("..", "data", "classifiers", f"{clf}.pkl")
    clf = joblib.load(p)
    X = X.reindex(columns=clf.feature_names_in_)

    # make prediction
    label = clf.predict(X)
    return label[0]

def generate_features(df: pd.DataFrame, features: list)-> (pd.DataFrame, list):
    feature_cols = ["Start", "End"]
    if "Segment" in features:
        df, segment_cols = segment_ohe(df)
        feature_cols = feature_cols + segment_cols
    if "DI_Length" in features:
        df["DI_Length"] = df.apply(get_dirna_length, axis=1)
        feature_cols.append("DI_Length")
    if "Direct_repeat" in features:
        df["Direct_repeat"] = df.apply(get_direct_repeat_length, axis=1)
        feature_cols.append("Direct_repeat")
    if "Junction" in features:
        df, junction_start_cols = junction_site_ohe(df, "Start")
        df, junction_end_cols = junction_site_ohe(df, "End")
        feature_cols = feature_cols + junction_start_cols + junction_end_cols
    if "3_5_ratio" in features:
        df["3_5_ratio"] = df.apply(get_3_to_5_ratio, axis=1)
        feature_cols.append("3_5_ratio")
    if "length_proportion" in features:
        df["length_proportion"] = df.apply(get_length_proportion, axis=1)
        feature_cols.append("length_proportion")
    if "delta_G" in features:
        df["delta_G"] = df.apply(get_delta_G, axis=1)
        feature_cols.append("delta_G")
    if "Peptide_Length" in features:
        df["Peptide_Length"] = df.apply(get_peptide_length, axis=1)
        feature_cols.append("Peptide_Length")

    return df, feature_cols

def segment_ohe(df: pd.DataFrame)-> (pd.DataFrame, list):
    ohe = OneHotEncoder()
    segment_df = pd.DataFrame(ohe.fit_transform(df[["Segment"]]).toarray())
    ohe_cols = ohe.get_feature_names_out().tolist()
    segment_df.columns = ohe_cols
    df = df.join(segment_df)
    return df, ohe_cols

def get_dirna_length(row: pd.Series)-> int:
    seq_len = len(row["Sequence"])
    return row["Start"] + (seq_len - row["End"] + 1)

def get_direct_repeat_length(row: pd.Series)-> int:
    n = calculate_direct_repeat(row["Sequence"], row["Start"], row["End"])
    return n

def calculate_direct_repeat(seq: str, s: int, e: int)-> (int, str):
    w_len = 15
    start_window = seq[s-w_len: s]
    end_window = seq[e-1-w_len: e-1]
    # if they are the same return directly to avoid off-by-one error
    if start_window == end_window:
        return len(start_window)

    counter = 0
    for i in range(w_len - 1, -1, -1):
        if start_window[i] == end_window[i]:
            counter += 1
        else:
            break

    return counter

def junction_site_ohe(df: pd.DataFrame, pos: str)-> (pd.DataFrame, list):
    res = np.zeros((df.shape[0], CHARS_COUNT * 10), dtype=np.uint8)

    for i, r in df.iterrows():
        seq = r["Sequence"][r[pos]-5:r[pos]+5]
        for j, char in enumerate(seq):
            res[i][j*CHARS_COUNT+CHARS.rfind(char)] = 1

    # format as data frame and create columns names of OHE
    encoded_df = pd.DataFrame(res)
    col_names = [f"{pos}_{i}_{ch}" for i in range(1, 11) for ch in CHARS]
    encoded_df.columns = col_names
    return df.join(encoded_df), col_names

def get_3_to_5_ratio(row: pd.Series)-> float:
    seq_len = len(row["Sequence"])
    len3 = row["Start"]
    len5 = seq_len - row["End"] + 1
    return len3/len5

def get_length_proportion(row: pd.Series)-> float:
    seq_len = len(row["Sequence"])
    dirna_len = row["Start"] + (seq_len - row["End"] + 1)
    return dirna_len/seq_len

def full_sequence_ohe(df: pd.DataFrame)-> (pd.DataFrame, list):
    res = np.zeros((df.shape[0], CHARS_COUNT * MAX_LEN), dtype=np.uint8)

    for i, r in df.iterrows():
        seq = r["Sequence"]
        seq = seq + "*" * (MAX_LEN - len(seq))
        for j, char in enumerate(seq):
            res[i][j*CHARS_COUNT+CHARS.rfind(char)] = 1

    # format as data frame and create columns names of OHE
    encoded_df = pd.DataFrame(res)
    col_names = [f"{i}_{ch}" for i in range(1, MAX_LEN+1) for ch in CHARS]
    encoded_df.columns = col_names
    return df.join(encoded_df), col_names

def get_delta_G(row: pd.Series)-> float:
    seq = row["Sequence"]
    del_seq = seq[:row["Start"]] + seq[row["End"]-1:]
    mfe = RNA.fold_compound(del_seq).mfe()[1]
    return mfe/len(del_seq)

def get_peptide_length(row: pd.Series)-> float:
    strain = row["Strain"]
    segment = row["Segment"]
    seq = row["Sequence"]
    f = open(os.path.join("py", "translation_indices.json"))
    indices = json.load(f)
    f.close()
    s = indices[strain][segment]["start"]-1
    e = indices[strain][segment]["end"]
    seq_obj = Seq(seq[s:row["Start"]] + seq[row["End"]-1:e])
    pep_seq = seq_obj.translate(to_stop=True)
    return len(pep_seq)

