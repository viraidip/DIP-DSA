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
    '''

    '''
    features = ["Segment", "DI_Length", "Direct_repeat","Junction", "3_5_ratio", "length_proportion", "delta_G", "Peptide_Length"]

    strain_dict = dict({"A/California/07/2009": "Cal07"
                        })
    strain = strain_dict[strain]

    # create df
    data = dict({"Start": s, "End": e, "Strain": strain, "Segment": segment, "Sequence": sequence})
    df = pd.DataFrame(data, index=[1])


    # create artificial samples to let segment OHE run correctly
    sgms = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
    n_df = [{"Start": 0,"End":1,"Segment":seg,"Strain":strain, "Sequence": sequence} for seg in sgms]
    df = df.append(n_df, ignore_index=True)

    # create features
    df, feature_cols = generate_features(df, features)

    # remove artificial samples
    df = df[df["Start"] != 0]
    X = df[feature_cols]

    # load pickle file
    p = os.path.join("..", "data", "classifiers", f"{clf}.pkl")
    clf = joblib.load(p)

    # make prediction
    label = clf.predict(X)
    return label[0]


def generate_features(df: pd.DataFrame,
                      features: list
                      )-> (pd.DataFrame, list):
    '''
        Main function to generate/load the features.
        :param df: data frame including the deletion site data
        :param features: list indicating which features to include
        :param load_precalc: if True features are loaded from precalculated
                             file
        :return: data frame including the features,
                 list with the names of the columns for the features
    '''
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
    '''
        Calculates the length of the DI RNA sequence given a row of a data
        frame with the necessary data.
        :param row: data frame row including Strain, Segment, Start, and End

        :return: length of DI RNA sequence
    '''
    seq_len = len(row["Sequence"])
    return row["Start"] + (seq_len - row["End"] + 1)

def get_direct_repeat_length(row: pd.Series)-> int:
    '''
        Calculates the length of the direct repeat given a row of a data frame
        with the necessary data.
        :param row: data frame row including Strain, Segment, Start, and End

        :return: length of direct repeat
    '''
    n, _ = calculate_direct_repeat(row["Sequence"], row["Start"], row["End"], 15, 1)
    return n

def calculate_direct_repeat(seq: str,
                            s: int,
                            e: int,
                            w_len: int,
                            m: int
                            )-> (int, str):
    '''
        counts the number of overlapping nucleotides directly before start and
        end of junction site --> direct repeats
        :param seq: nucleotide sequence
        :param w_len: length of window to be searched
        :param s: start point
        :param e: end point
        :param m: mode how the overlap is created
                    full junction site is 1234J6789
                    1: S 1234J | E 1234J (overlap beginning of both sites)
                    2: S 1234J | E J6789 (overlap what stays in deletion seq)
        :return: Tuple with two entries:
                    Integer giving the number of overlapping nucleotides
                    String of the overlapping nucleotides
    '''
    counter = 0

    if m == 1:
        start_window = seq[s-w_len: s]
        end_window = seq[e-1-w_len: e-1]
        
        #if they are the same return directly to avoid off-by-one error
        if start_window == end_window:
            return len(start_window), start_window

        for i in range(w_len - 1, -1, -1):
            if start_window[i] == end_window[i]:
                counter += 1
            else:
                break
        overlap_seq = str(start_window[i+1:w_len])

    elif m == 2:
        for i in range(w_len):
            if seq[s-i:s] == seq[e-1:e-1+i]:
                counter = i
                overlap_seq = str(seq[s-i:s])

    assert counter == len(overlap_seq), f"{counter=}, {len(overlap_seq)}"
    if len(overlap_seq) == 0:
        overlap_seq = "_"

    return counter, overlap_seq

def junction_site_ohe(df: pd.DataFrame,
                      position: str
                      )-> (pd.DataFrame, list):
    '''
        Gets the sequence around the start or end of a given junction site and
        converts the sequence into an one hot encoding.
        :param df: data frame including Start, End, Strain, and Segment
        :param position: is either 'Start' or 'End' to indicate which site

        :return: Tuple with two entries:
                    data frame including original data and OHE data
                    list with the column names of the OHE
    '''
    # initializing matrix
    n = df.shape[0]
    res = np.zeros((n, CHARS_COUNT * 10), dtype=np.uint8)


    # getting sequence window for each row and do OHE
    for i, r in df.iterrows():
        s = r[position]
        sequence = r["Sequence"]
        seq = sequence[s-5:s+5]
        # Write down OHE in numpy array
        for j, char in enumerate(seq):
            pos = CHARS.rfind(char)
            res[i][j*CHARS_COUNT+pos] = 1

    # format as data frame and create columns names of OHE
    encoded_df = pd.DataFrame(res)
    col_names = [f"{position}_{i}_{ch}" for i in range(1, 11) for ch in CHARS]
    encoded_df.columns = col_names
    df = df.join(encoded_df)

    return df, col_names

def get_3_to_5_ratio(row: pd.Series)-> float:
    '''
        Calculates the proportion of the 3' sequence to the 5' sequence given
        a row of a data frame.
        :param row: data frame row including Strain, Segment, Start, and End

        :return: ratio of 3' to 5' sequence length
    '''
    seq_len = len(row["Sequence"])
    len3 = row["Start"]
    len5 = seq_len - row["End"] + 1
    return len3/len5

def get_length_proportion(row: pd.Series)-> float:
    '''
        Calculates the proportion of the length of the DI RNA sequence to the
        full length sequence given a row of a data frame.
        :param row: data frame row including Strain, Segment, Start, and End

        :return: ratio of DI RNA lenght to full length sequence
    '''
    seq_len = len(row["Sequence"])
    dirna_len = row["Start"] + (seq_len - row["End"] + 1)
    return dirna_len/seq_len

def full_sequence_ohe(df: pd.DataFrame)-> (pd.DataFrame, list):
    '''
        Gets the whole sequence as an one hot encoding. Sequences get
        normalized to the longest sequence length by adding * at the end
        :param df: data frame including Start, End, Strain, and Segment

        :return: Tuple with two entries:
                    data frame including original data and OHE data
                    list with the column names of the OHE
    '''
    # defining initializing matrix
    n = df.shape[0]
    res = np.zeros((n, CHARS_COUNT * MAX_LEN), dtype=np.uint8)

    # getting sequence window for each row and do OHE
    for i, r in df.iterrows():
        seq = r["Sequence"]
        seq = seq + "*" * (MAX_LEN - len(seq))
        # Write down OHE in numpy array
        for j, char in enumerate(seq):
            pos = CHARS.rfind(char)
            res[i][j*CHARS_COUNT+pos] = 1

    # format as data frame and create columns names of OHE
    encoded_df = pd.DataFrame(res)
    col_names = [f"{i}_{ch}" for i in range(1, MAX_LEN+1) for ch in CHARS]
    encoded_df.columns = col_names
    df = df.join(encoded_df)

    return df, col_names

def get_delta_G(row: pd.Series)-> float:
    '''
        :param row: data frame row including Strain, Segment, Start, and End

        :return: ratio of DI RNA lenght to full length sequence
    '''
    seq = row["Sequence"]
    del_seq = seq[:row["Start"]] + seq[row["End"]-1:]
    mfe = RNA.fold_compound(del_seq).mfe()[1]
    return mfe/len(del_seq)

def get_peptide_length(row: pd.Series)-> float:
    '''
        :param row: data frame row including Strain, Segment, Start, and End

        :return: length of resulting peptid
    '''
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

