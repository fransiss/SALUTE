import pickle
import math
import unicodedata
import os
import pandas as pd
import numpy as np
from statistics import mode
from collections import Counter
import sklearn
import math

def get_first_mode(a):
    c = Counter(a)  
    mode_count = max(c.values())
    mode = {key for key, count in c.items() if count == mode_count}
    first_mode = next(x for x in a if x in mode)
    return first_mode

def BetaScore(x,n):
    a1=2+x
    b1=2+n-x
    scoreBayesiano=(a1-1)/(a1+b1-2)
    return scoreBayesiano      

def initialize_results():
    return {
        'drugs2direction_all': {},
        'drugs2direction_clinical': {},
        'drugs2points_RS': {},
        'drugs2points_RS_clinical': {},
        'drugs2negative_clinical': {},
        'drugs2direction_CS_clinical': {},
        'drugs2negative': {},
        'drugs2direction_CS': {},
        'drugs2positive_clinical': {},
        'drugs2tot_clinical': {},
        'drugs2tot': {},
        'drugs2positive': {},
        'drugs2num_positive': {},
        'drugs2num_negative': {},
        'drugs2num_positive_clinical': {},
        'drugs2num_negative_clinical': {},
        'drugs2direction_clinistr': {}
    }

def parse_pointstr(pointstr):
    if pointstr != 'NA':
        id2pointslist = pointstr.split(';')
        id2points = dict([(sp.split(':')[0], int(sp.split(':')[1])) for sp in id2pointslist])
    else:
        id2points = {}
    return id2points


def calculate_relevance_score(all_points):
    scors = 1
    all_points_sc = [s for s in all_points if s != 0]
    for sco in all_points_sc:
        scorp = sco / 100
        scors *= (1 - scorp)
    return 1 - scors

#fucntion calculate confidence score
def calculate_scores(pos, neg, tot):
    x = len(pos)
    xneg = len(neg)
    n = x + xneg
    scoreBayes_all = BetaScore(x, n)
    direction = '+'
    if scoreBayes_all < 0.5:
        scoreBayes_all = 1 - scoreBayes_all
        direction = '0'
    elif scoreBayes_all == 0.5:
        direction = '0'
    if len(pos) > len(neg):
        direction = '+'
    elif len(pos) < len(neg):
        direction = '-'
    else:
        direction = '0'
    return x, xneg, n, direction, scoreBayes_all

# if the drug has no clincal ids - reset the dictionaries to 0 or empty
def reset_clinical_results(drug, results):
    results['drugs2direction_clinical'][drug] = '0'
    results['drugs2direction_CS_clinical'][drug] = 0
    results['drugs2num_positive_clinical'][drug] = 0
    results['drugs2num_negative_clinical'][drug] = 0
    results['drugs2tot_clinical'][drug] = 0
    results['drugs2negative_clinical'][drug] = ''
    results['drugs2positive_clinical'][drug] = ''
    
# Function to update dataframe with results
def update_dataframe(df, results, clinical):
    if clinical:
        df["Direction"] = df["Drug"].map(results['drugs2direction_clinical'])
        df["RS"] = df["Drug"].map(results['drugs2points_RS_clinical'])
        df["CS"] = df["Drug"].map(results['drugs2direction_CS_clinical'])
        df["NegativeTot"] = df["Drug"].map(results['drugs2num_negative_clinical'])
        df["PositiveTot"] = df["Drug"].map(results['drugs2num_positive_clinical'])
        df["clinicalTot"] = df["Drug"].map(results['drugs2tot_clinical'])
    else:
        df["Direction"] = df["Drug"].map(results['drugs2direction_all'])
        df["RS"] = df["Drug"].map(results['drugs2points_RS'])
        df["CS"] = df["Drug"].map(results['drugs2direction_CS'])
        df["NegativeTot"] = df["Drug"].map(results['drugs2num_negative'])
        df["PositiveTot"] = df["Drug"].map(results['drugs2num_positive'])
        df["num_paperTot"] = df["Drug"].map(results['drugs2tot'])

def save_to_excel(df, filepath, results, clinical):
     with pd.ExcelWriter(filepath) as writer:
        directions = df["Direction"].unique()
        for direction in directions:
            df_direction = df[df["Direction"] == direction]
            df_direction.to_excel(writer, sheet_name=f"{direction}", index=False)


def load_data_if_none(drugs2posneg, drug2pubid_clinical, df_clinical, df_allDrug, fp_dict, fpath_txt):
    if drugs2posneg is None:
        drugs2posneg = pickle.load(open(os.path.join(fp_dict, 'drugs2posneg'), 'rb'))

    if drug2pubid_clinical is None:
        drug2pubid_clinical = pickle.load(open(os.path.join(fp_dict, 'drug2article_clinical'), 'rb'))

    if df_clinical is None:
        df_clinical = pd.read_excel(os.path.join(fpath_txt, "Drugs_Clinical_papers_pubmedID.xlsx"))

    if df_allDrug is None:
        df_allDrug = pd.read_excel(os.path.join(fpath_txt, "Drugs_papers_pubmedID.xlsx"))

    return drugs2posneg, drug2pubid_clinical, df_clinical, df_allDrug


def main(drugs2posneg=None, drug2pubid_clinical=None, df_clinical=None, df_allDrug=None):
    path_base = os.getcwd()
    fpath_txt = os.path.join(path_base, 'txt')
    fp_dict = os.path.join(path_base, 'dictionary')
    drugs2posneg, drug2pubid_clinical, df_clinical, df_allDrug= load_data_if_none(drugs2posneg, drug2pubid_clinical, df_clinical, df_allDrug, fp_dict, fpath_txt)
    results = initialize_results()
    
    for drug,lists in drugs2posneg.items():
        pos, neg, pointstr, tot = lists[0], lists[1], lists[2], lists[3]
        if tot==0:
            continue
        id2points = parse_pointstr(pointstr)

        if drug in drug2pubid_clinical.keys():
            clinical_ids=drug2pubid_clinical[drug]
        
            if isinstance(clinical_ids, str):
                clinical_ids_split=clinical_ids.split(',')
                all_points_clin = [id2points[clin] for clin in clinical_ids_split if clin in id2points]
                
                relevanceScore_clinical = calculate_relevance_score(all_points_clin)

                results['drugs2points_RS_clinical'][drug]=relevanceScore_clinical
                
                clinical_ids_w_abst=set(clinical_ids.split(',')).intersection(id2points.keys())
                if clinical_ids_w_abst:
                    
                    pos_clin = clinical_ids_w_abst.intersection(pos)
                    neg_clin = clinical_ids_w_abst.intersection(neg)
    
                    x_clinical = len(pos_clin)
                    x_clinical_neg = len(neg_clin)
                    n_clinical = x_clinical + x_clinical_neg
                    
            
                    scoreBayes_clinical=BetaScore(x_clinical,n_clinical)
                    direction_clin='+'
                    direction_clin = '+' if scoreBayes_clinical >= 0.5 else '-'
                    if scoreBayes_clinical==0.5:
                        direction_clin='0'
    
                    results['drugs2direction_clinical'][drug] = direction_clin
                    results['drugs2direction_CS_clinical'][drug] = scoreBayes_clinical
                    results['drugs2num_positive_clinical'][drug] = x_clinical
                    results['drugs2num_negative_clinical'][drug] = x_clinical_neg
                    results['drugs2tot_clinical'][drug] = n_clinical
                    results['drugs2negative_clinical'][drug] = ','.join([f for f in neg_clin])
                    results['drugs2positive_clinical'][drug] = ','.join([f for f in pos_clin])
                else:
                   reset_clinical_results(drug, results)

        all_points = [id2points[idsl] for idsl in id2points.keys()]
        relevanceScore = calculate_relevance_score(all_points)
        results['drugs2points_RS'][drug]=relevanceScore
    
        x, xneg, n, direction, scoreBayes_all = calculate_scores(pos, neg, tot)
        results['drugs2direction_all'][drug] = direction
        results['drugs2direction_CS'][drug] = scoreBayes_all
        results['drugs2num_positive'][drug] = x
        results['drugs2num_negative'][drug] = xneg
        results['drugs2tot'][drug] = n
        results['drugs2negative'][drug] = ','.join([f for f in neg])
        results['drugs2positive'][drug] = ','.join([f for f in pos])
    
    update_dataframe(df_clinical, results, clinical=True)
    update_dataframe(df_allDrug, results, clinical=False)
    
    df_clinical = df_clinical[df_clinical['#paper'] > 0]
    save_to_excel(df_clinical, os.path.join(fpath_txt, "Clinical_Drugs_scores.xlsx"), results, clinical=True)
    save_to_excel(df_allDrug, os.path.join(fpath_txt, "all_Drugs_Scores.xlsx"), results, clinical=False)
    
    print( "Clinical_Drugs_scores.xlsx and all_Drugs_Scores.xlsx have been saved to the path \n"+ fpath_txt)
    
