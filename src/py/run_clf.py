import os
import joblib

def run_classification(s, e, strain, segment, clf)-> str:
    '''

    '''
    # create df

    # create features

    # load pickle file
    p = os.path.join("..", "data", "classifiers", f"{clf}.pkl")
    clf = joblib.load(p)

    # make prediction

    label = "low"

    return label
