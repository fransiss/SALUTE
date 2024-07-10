import matplotlib.pyplot as plt
import matplotlib
import os
import orangecontrib.text.pubmed as pb
import pickle
import metapub
from collections import Counter
import xlsxwriter
from Bio import Entrez
from Bio import Medline
import pandas as pd
from wordcloud import WordCloud
import re
from nltk.tokenize import word_tokenize 
from nltk.corpus import stopwords
import nltk
import analyze_and_compute_scores

import unicodedata
import numpy as np
import sklearn
import math


def setup_matplotlib():
    font = {'family': 'sans-serif', 'sans-serif': 'Arial', 'weight': 'bold'}
    matplotlib.rc('font', **font)
    plt.rcParams['svg.fonttype'] = 'none'
    
def load_DoI_terms(file_path):
   df = pd.read_csv(file_path, header=None)
   return list(df[0])

def load_Drug2Indications(file_path):
   df = pd.read_csv(file_path,  encoding='ISO-8859-1',sep='\t')
   # txt file with  drug name int the first column and the indications separated by comma
   drug2indication_list={}
   for dr, inidstr in zip(df['Drug'],df['Indications']):
       drug2indication_list[dr]=inidstr.split(',')
   return drug2indication_list


def build_query_string(terms, term_factor):
    query_parts = []
    for term in terms:
        query_part = f'({term}[MeSH Terms] AND {term_factor}[MeSH Terms]) OR ({term}[Title/Abstract] AND {term_factor}[Title/Abstract])'
        query_parts.append(query_part)
    return ' OR '.join(query_parts)

def process_pubmed_records(fetch,ids_all, qualifier2mesh, qualifier2keep):
    all_meshs_ret = []
    for nni, ii in enumerate(ids_all):
        review = 'no'
        handle = Entrez.efetch(db="pubmed", id=ii, rettype="medline", retmode="text")
        records = list(Medline.parse(handle))
        if any('Review' in record.get('PT', []) for record in records):
            review = 'yes'
            #do not consider reviews
        if review == 'no':
            try:
                article = fetch.article_by_pmid(ii)
            except:
                article = fetch.article_by_pmid(ii)

            if article.mesh:
                for mesh in article.mesh.values():
                    desc = mesh['descriptor_name']
                    qualifier = mesh.get('qualifier_name', 'None')
                    if qualifier in qualifier2keep:
                        qualifier2mesh.setdefault(qualifier, []).append(desc)
                        all_meshs_ret.append(desc)

    return all_meshs_ret

def save_mesh_counts(all_meshs, output_path):
    count_meshes = Counter(all_meshs)
    with open(output_path, 'w') as file:
        file.write('MeshTerm\tcounting\n')
        for mesh, count in count_meshes.items():
            file.write(f'{mesh}\t{count}\n')
            
def save_mesh_to_excel(qualifier2mesh, output_path):
    workbook = xlsxwriter.Workbook(output_path)
    for qualifier, mesh_list in qualifier2mesh.items():
        worksheet = workbook.add_worksheet(qualifier if qualifier != 'history' else 'history_1')
        worksheet.write(0, 0, 'MeshTerm')
        worksheet.write(0, 1, 'Count')
        for row, (mesh, count) in enumerate(Counter(mesh_list).items(), start=1):
            worksheet.write(row, 0, mesh)
            worksheet.write(row, 1, count)
    workbook.close()
            
def save_data(fp_dict, ii_done, no_mesh, ids_notREV, qualifier2mesh, all_meshs):
    def save_set(data_set, file_name):
        with open(os.path.join(fp_dict, file_name), 'w') as file:
            file.write('\n'.join(data_set))

    def save_pickle(data, file_name):
        with open(os.path.join(fp_dict, file_name), 'wb') as file:
            pickle.dump(data, file)

def contains_keywords_or_patterns(sentence, keywords, regex_patterns):
    return any(word in sentence for word in keywords) or any(re.search(pattern, sentence) for pattern in regex_patterns)
    


if __name__ == "__main__":
    setup_matplotlib()
    
    #initialize pubmed with thr 
    youremail='francescavitali@email.arizona.edu'
    pubmed=pb.Pubmed(youremail)
    
    path_base = os.getcwd()
    fpath_txt = os.path.join(path_base, 'txt')
    fp_dict = os.path.join(path_base, 'dictionary')
    fpath_images=os.path.join(path_base, 'figures')
    fp_cache = os.path.join(path_base, 'cache')
    
    term_factor = 'risk factor*'
    DoI_terms = load_DoI_terms(os.path.join(fpath_txt, "disease_terms.txt"))
    
    doi_key = ["stroke", "cerebral vascular accident"]
    doi_regex_patterns = [r'\bcva\b']
    
    
    #keep only the mesh terms beloning to these categories
    mesh_qualifier2keep=['prevention & control','drug therapy','therapy','diagnosis']
    
    # query PubMed with DoI terms and the risk factors to identify 
    query_string = build_query_string(DoI_terms, term_factor)

    PUBMED_FIELD_HEADINGS = 'MESH headings'
    
    fetch = metapub.PubMedFetcher()
    _, pubmed_ids_all = pubmed._search_for_records(advanced_query=query_string)

    qualifier2mesh={}
    all_meshs=process_pubmed_records(fetch, pubmed_ids_all, qualifier2mesh, mesh_qualifier2keep)
    #save txt file
    outp=os.path.join(fpath_txt,'mesh_counting.xlsx')
    save_mesh_to_excel(qualifier2mesh,outp)
    
    #generate world cloud 
    all_mesh_risk_factors=[mesh for qual,meshv in qualifier2mesh.items() for  mesh in meshv]  
    
    word_could_dict=Counter(all_mesh_risk_factors)
    wordcloud = WordCloud(font_path='arial',width = 2000, height = 1000, background_color="white", colormap="Blues").generate_from_frequencies(word_could_dict)

    plt.figure(figsize=(15,8))
    plt.imshow(wordcloud)
    plt.axis("off")
    plt.savefig(os.path.join(fpath_images,'wordcould.png'), dpi=300)
    plt.close()

    
    
    #load Drug Indication database - here txt file provided  - function return a dictionary
    drug2indications= load_Drug2Indications(os.path.join(fpath_txt, "drug_indication_database.txt"))
    #dictionary from mesh risk factor if it is an indication of a drug
    RF_mesh2drugs={}
    for drug, indicv in drug2indications.items():
        for mm in all_mesh_risk_factors:
            for indic in indicv:
               if  mm.lower() in indic.lower():
                   drugc=RF_mesh2drugs.get(mm, [])
                   drugc.append(drug)
                   RF_mesh2drugs[mm]=list(set(drugc))
                   

    #all drugs found associated to at least on Risk Factors
    all_drugs=list(set([dr for drv in RF_mesh2drugs.values() for dr in drv]))
    # all_drugs=['Certoparin']
    
    positive_words=['delay','lower', 'decrease', 'protective','reduce','improve','reduced','reduction', 'preventing','prevent','preventative']
    negative_words=['increase', 'elevated','elevate', 'increased','worse']
    strtoremove=r"""\.|,|:|;|!|\?|\(|\)|\||\+|'|"|‘|’|“|”|'|\’|…|\-|–|—|\$|&|\*|>|<|\/|\[|\]|\n|\|%|=|"""
    regex = re.compile(strtoremove)
    nltk.download('punkt')
    nltk.download('stopwords')
    stopWords = set(stopwords.words('english'))
    
    # save excel file with drugs - number of paper relating the drug to the disease of interest and the related Pubmed IDs
    workbook = xlsxwriter.Workbook(os.path.join(fpath_txt,"Drugs_papers_pubmedID.xlsx")) 
    worksheet = workbook.add_worksheet()
    bold = workbook.add_format({'bold': True})
 
    first_row =['Drug', '#paper', 'paperIDs']
    for col_num, data in enumerate(first_row):
        worksheet.write(0, col_num, data, bold)
    
    # save excel file with drugs - number of CLINICAL paper relating the drug to the disease of interest and the related Pubmed IDs
    workbook1 = xlsxwriter.Workbook(os.path.join(fpath_txt,"Drugs_Clinical_papers_pubmedID.xlsx")) 
    worksheet1 = workbook1.add_worksheet()
    bold = workbook1.add_format({'bold': True})
    
    first_row =['Drug', '#paper', 'paperIDs']
    for col_num, data in enumerate(first_row):
        worksheet1.write(0, col_num, data, bold)
    i=1 
    
    terms='Stroke*'
    
    ii=1
    drug2articles={} 
    drug2article_clinical={}
    drugs2posneg={}
    for nnnd,drug in enumerate(all_drugs):
        drug_v=drug.split('+')
        string_drugMESH = ' AND '.join([f'"{drugc.strip(" ")}"[MeSH Terms]' for drugc in drug_v])
        string_drugTIT = ' AND '.join([f'"{drugc.strip(" ")}"[Title/Abstract]' for drugc in drug_v])
        query = f'(({terms}[MeSH Terms] AND {string_drugMESH}) OR ({terms}[Title/Abstract] AND {string_drugTIT}))'
        numart, ids = pubmed._search_for_records(advanced_query=query)
        ids_positive, ids_negative, ids_incognito = set(), set(), set()
        numart_clinical = []
        points = 0
        if ids:
            pbIds2abstract={}
            ids2points={}
            for idPub in ids:
                try:
                    handle = Entrez.efetch(db="pubmed", id=idPub, rettype="medline", retmode="text")
                    records = Medline.parse(handle)
                    for record in records:
                        abstr = record.get('AB', '')
                        titl = record.get('TI', '')
                        type_pub = record.get('PT', [])
                        if 'Review' in type_pub:
                            titl, abstr = '', ''
                except Exception as e:
                    print(f'Error fetching record: {e}')
            
                abstr = abstr.replace('"', '')
                if re.search(rf' {drug.lower()} ', abstr.lower()):
                    no_drug_include = 'include'
                else:
                    no_drug_include = 'exclude'
                print(idPub)
                # print(no_drug_include)
                if abstr and titl and no_drug_include == 'include':
                   points, drugPoints, adPoints = 0, 0, 0
                   # split by sentence divided by dot - but not if the dot is related to a numner 
                   periods = [p for p in re.sub(r"(\d)\.(\d)", r"\1[DOT]\2", abstr).split('.') if p]
                   if len(periods)==1:
                       print(abstr)
                       periods=abstr.lower().split('discussion')
                       first_sentence = periods[0]
                       last_sentence = periods[-1]
                       last3periods_str = ' '.join(periods[-3:])
                       last3periods = periods[-3:]
                   else:
                       first_sentence = periods[0]
                       last_sentence = periods[-1]
                       last3periods_str = ' '.join(periods[-3:])
                       last3periods =periods[-3:]
                   
                   if drug.lower() in titl or drug in titl:
                       drugPoints += 50
                   elif drug.lower() in last_sentence.lower() or drug in last_sentence or drug.lower() in first_sentence.lower() or drug in first_sentence:
                       drugPoints += 30
                   
                   last_sentence_lower=last_sentence.lower()
                   first_sentence_lower=first_sentence.lower()
                   if any(word in titl.lower() for word in doi_key) or any(re.search(pattern, titl.lower()) for pattern in doi_regex_patterns):
                       adPoints += 25
                   elif (contains_keywords_or_patterns(last_sentence_lower, doi_key, doi_regex_patterns) or 
                         contains_keywords_or_patterns(first_sentence_lower, doi_key, doi_regex_patterns)):
                       adPoints += 15
                   points_to_add=0
                   for period in periods:
                       if any(word in period.lower() for word in doi_key) or any(re.search(pattern, period.lower()) for pattern in doi_regex_patterns):
                        words_ok = word_tokenize(period)
                        period = ' '.join([w for w in words_ok if w.lower() not in stopWords])
                        words = period.partition("risk")
                        if words[1]:
                            before5 = ' '.join(words[0].split()[-5:])
                            after5 = ' '.join(words[2].split()[:5])
                            if (contains_keywords_or_patterns(before5.lower(), doi_key, doi_regex_patterns) or
                                 contains_keywords_or_patterns(after5.lower(), doi_key, doi_regex_patterns)):
                                points_to_add = 25
                        if not points_to_add:
                            words = period.partition("diagnosis")
                            if words[1]:
                                before5 = ' '.join(words[0].split()[-5:])
                                after5 = ' '.join(words[2].split()[:5])
                                if (contains_keywords_or_patterns(before5.lower(), doi_key, doi_regex_patterns) or
                                     contains_keywords_or_patterns(after5.lower(), doi_key, doi_regex_patterns)):
                                    points_to_add = 25
                
                   adPoints=adPoints+points_to_add
                   points=adPoints+drugPoints
                   ids2points[idPub] = points
                  
                    #check drug in first last sentence\
                   print(drug_v[0])
                   print(idPub)
                   print(str(points) + 'points')

                   ids2points[idPub]=points
                   points_str = ';'.join([f'{pp}:{p}' for pp, p in ids2points.items()])
                    
                   titl = regex.sub('', titl)
                   positive, negative = 'no', 'no'
                   for lastPer in last3periods:
                        words = word_tokenize(lastPer)
                        lastPer = ' '.join([w for w in words if w.lower() not in stopWords])
                        for p in positive_words:
                            if p in lastPer:
                                for spli in re.split(p, lastPer):
                                    before5 = ' '.join(spli.split()[-5:])
                                    after5 = ' '.join(re.split(p, lastPer)[1].split()[:5])
                                    if (contains_keywords_or_patterns(before5.lower(), doi_key, doi_regex_patterns) or
                                         contains_keywords_or_patterns(after5.lower(), doi_key, doi_regex_patterns)):
                                        positive = 'yes'
                        for p in negative_words:
                            if p in lastPer:
                                for spli in re.split(p, lastPer):
                                    before5 = ' '.join(spli.split()[-5:])
                                    after5 = ' '.join(re.split(p, lastPer)[1].split()[:5])
                                    if (contains_keywords_or_patterns(before5.lower(), doi_key, doi_regex_patterns) or
                                         contains_keywords_or_patterns(after5.lower(), doi_key, doi_regex_patterns)):
                                        negative = 'yes'
                         
                   if positive=='yes':
                        ids_positive.add(idPub)
                   if negative=='yes':
                        ids_negative.add(idPub)
                        
                   if abstr!='':
                        abstr=regex.sub('',abstr)
                        ab_strings_all=abstr.split()
                   else:
                        ab_strings_all=''
                        
                   if titl:
                        titl_strings=titl.split()
                   else:
                        titl_strings=''
                        
                   upper_title=[word for word in titl_strings if word.isupper()]
                   upper_abstr=[word for word in ab_strings_all if word.isupper()]
                    
                   lower_title=[word for word in titl_strings if word not in upper_title]
                   lower_abstr=[word for word in ab_strings_all if word not in upper_abstr]
                    
                   lower_title=' '.join([s.lower() for s in lower_title])
                   lower_abstr=' '.join([s.lower() for s in lower_abstr])
                   jj=0
                   cond2='HR' in upper_title or 'HR' in upper_abstr
                   cond3='RR' in upper_title or 'RR' in upper_abstr
                   cond3='relative risk' in lower_title or 'relative risk' in lower_abstr
                   cond4='hazard ratio' in lower_title or 'hazard ratio' in lower_abstr
                   cond5='odd ratio' in lower_title or 'odd ratio' in lower_abstr
                   if (cond2) or (cond3) or (cond4) or (cond5) :
                        numart_clinical.append(idPub)
                   numart_cl=len(list(numart_clinical))
                   print(numart_cl)               
                else:
                    numart_cl=0
                    points_str='NA'
    
            jj=0 
            worksheet1.write(ii, jj, drug)
            jj=jj+1
            worksheet1.write(ii, jj, numart_cl)
            jj=jj+1
            worksheet1.write(ii, jj, ','.join([t for t in numart_clinical]))
            ii=ii+1
            drug2article_clinical[drug]=numart_clinical
            
            
            numart=len(ids)    
            drug2articles[drug]=(numart, ids)
            j=0
            
            worksheet.write(i, j, drug)
            j=j+1
            worksheet.write(i, j, numart)
            j=j+1
            worksheet.write(i, j, ','.join([t for t in ids]))
            i=i+1
            lista=[list(ids_positive),list(ids_negative),points_str, len(ids)]
            drugs2posneg[drug]=lista
            
    workbook.close()
    workbook1.close()
    pickle.dump(drugs2posneg, open(os.path.join(fp_dict,'drugs2posneg'),'wb'))
    pickle.dump(drug2article_clinical, open(os.path.join(fp_dict,'drug2article_clinical'),'wb'))
    pickle.dump(drug2articles, open(os.path.join(fp_dict,'drug2articles'),'wb'))
    
    analyze_and_compute_scores.main(drugs2posneg, drug2article_clinical )