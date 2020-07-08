

#### import libraries
from Bio import motifs
from Bio.Seq import Seq
import os
import pandas as pd
import numpy as np
import re

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import linear_kernel
from nltk.corpus import PlaintextCorpusReader



#### import JASPAR matrices ####
################################
def JASPAR_matrix_import(directory):
    BCRANKfile = os.listdir(directory)
    MOTIF_files = []
    for x in BCRANKfile:
        if len(re.findall("\\.matrix$", x)) == 1:
            MOTIF_files.append(x)
    return MOTIF_files

#### identify transcription factor motifs that match identified motifs ####
###########################################################################
def motifannot(jspmotifdb, jspmotifdbopen, MOTIF_files, motif_dir):
    #### load file names obtained from JASPAR 
    file = os.listdir(jspmotifdb)#
    ##############################
    temp_df = pd.DataFrame() 
    for b in MOTIF_files:
        #### file from BCRANK in R
        with open(motif_dir + b) as handle:
           BCRANK = motifs.read(handle, "JASPAR")
        #### loop through files and make motif comparisons
        off = []
        dist = []
        PCC = []
        BC_name = []
        name_comp = []
        for f in file:
            #### open file ####
            with open(jspmotifdbopen + f) as handle:
               m = motifs.read(handle, "JASPAR")
            #### normalive psuedocounts and background between motifs ####
            m.background = {"A":0.3,"C":0.2,"G":0.2,"T":0.3} #m.background = BCRANK.background
            m.pseudocounts = {"A":0.6, "C": 0.4, "G": 0.4, "T": 0.6} #m.pseudocounts = BCRANK.pseudocounts
            BCRANK.background = {"A":0.3,"C":0.2,"G":0.2,"T":0.3}
            BCRANK.pseudocounts = {"A":0.6, "C": 0.4, "G": 0.4, "T": 0.6}
            #### Get Position Specific Scoring Matrices ###################
            BCRANK_PSSM = BCRANK.pssm
            m_PSSM = m.pssm
            #### make comparisons #########################################
            distance, offset = BCRANK_PSSM.dist_pearson(m_PSSM)
            #### Return offset, distance, and PCC #########################
            off.append(offset)# = print(offset)
            dist.append(distance)# = print(distance)
            PCC.append(1-distance)# = print("PCC =", 1-distance)
            #### return names #############################################
            BC_name.append(BCRANK.name + "_" + BCRANK.matrix_id)
            name_comp.append(m.name + "_" + m.matrix_id)
        #### transform into a data frame ####
        dic = {"query name":BC_name, "offset":off, "distance":dist, "PCC":PCC, "comparison_name":name_comp}
        df = pd.DataFrame(dic)
        #### find the motif that has the greatest match ####
        df = df.sort_values(by='PCC', ascending=False)
        df = df.reset_index()
        df = df[0:10]
        temp_df = temp_df.append(df, ignore_index=True)
    return temp_df

#### identify transcription factor motifs that match the reverse complement of the identified motifs ####
#########################################################################################################
def revcompmotifannot(jspmotifdb, jspmotifdbopen, MOTIF_files, motif_dir):
    #### load file names obtained from JASPAR
    file = os.listdir(jspmotifdb)#[0]
    temp_df = pd.DataFrame() 
    for b in MOTIF_files:
        #### file from BCRANK in R ###################################################################
        with open(motif_dir + b) as handle:
           BCRANK = motifs.read(handle, "JASPAR")
        ###############################################################################################

        #### loop through files and make motif comparisons #################
        off = []
        dist = []
        PCC = []
        BC_name = []
        name_comp = []
        for f in file:
            #### open file ####
            with open(jspmotifdbopen + f) as handle:
               m = motifs.read(handle, "JASPAR")
        #    #### take the reverse complement of JASPAR database motif
        #    m = m.reverse_complement()
            #### normalive psuedocounts and background between motifs ####
            m.background = {"A":0.3,"C":0.2,"G":0.2,"T":0.3} #m.background = BCRANK.background
            m.pseudocounts = {"A":0.6, "C": 0.4, "G": 0.4, "T": 0.6} #m.pseudocounts = BCRANK.pseudocounts
            BCRANK.background = {"A":0.3,"C":0.2,"G":0.2,"T":0.3}
            BCRANK.pseudocounts = {"A":0.6, "C": 0.4, "G": 0.4, "T": 0.6}
           #### Get Position Specific Scoring Matrices ###################
            BCRANK_PSSM = BCRANK.pssm
            m_PSSM = m.pssm
            #### take the reverse complement of JASPAR database motif
            m_PSSM= m_PSSM.reverse_complement()
            #### make comparisons #########################################
            distance, offset = BCRANK_PSSM.dist_pearson(m_PSSM)
            #### Return offset, distance, and PCC #########################
            off.append(offset)# = print(offset)
            dist.append(distance)# = print(distance)
            PCC.append(1-distance)# = print("PCC =", 1-distance)
            #### return names #############################################
            BC_name.append(BCRANK.name + "_" + BCRANK.matrix_id)
            name_comp.append(m.name + "_" + m.matrix_id)
        #### transform into a data frame ####
        dic = {"query name":BC_name, "offset":off, "distance":dist, "PCC":PCC, "comparison_name":name_comp}
        df = pd.DataFrame(dic)
        #df.columns = ['index', 'f']
        #### find the motif that has the greatest match ####
        df = df.sort_values(by='PCC', ascending=False)
        df = df.reset_index()
        df = df[0:10]
        temp_df = temp_df.append(df, ignore_index=True)
    return temp_df




##########################################################################################
#### Function that identifies what transcription factors correspond to selected terms ####
##########################################################################################
def TFIDFgenescore( directory, terms):
    #### open text files containing compiled GO term information ####
    #################################################################
    corpusdir = directory # Directory of corpus.
    newcorpus = PlaintextCorpusReader(corpusdir, '.*')
    #newcorpus

    #### create a corpus list ####
    ##############################
    fincorpus = [terms]
    for infile in sorted(newcorpus.fileids()):
        #print(infile) # The fileids of each file.
    #### Select method of opening individual corpus ####
        #newcorpus.open("1_AHR.txt").read().strip() # Opens the file. # Prints the content of the file
        corp = newcorpus.raw(infile).strip() # Access the plaintext; outputs pure string/basestring.
        # newcorpus.paras() # Each element in the outermost list is a paragraph, and Each paragraph contains sentence(s), and Each sentence contains token(s)
        # newcorpus.paras(newcorpus.fileids()[0]) # To access pargraphs of a specific fileid.
        # newcorpus.sents() # Access sentences in the corpus. (list of list of strings) # NOTE: That the texts are flattened into sentences that contains tokens.
        # newcorpus.sents(newcorpus.fileids()[0]) # To access sentences of a specific fileid.
        # newcorpus.words()# Access just tokens/words in the corpus. (list of strings)
        # newcorpus.words(newcorpus.fileids()[0]) # To access tokens of a specific fileid.
        fincorpus.append(corp)

    #### fit and transform newly created corpus ####
    ################################################
    tfidf = TfidfVectorizer().fit_transform(fincorpus)

    #### the cosine distances of one document (e.g. the first in the dataset) and all of the others you just need to 
    #### compute the dot products of the first vector with all of the others as the tfidf vectors are already row-normalized.
    #### Extract count features and apply TF-IDF normalization and row-wise euclidean normalization
    #from sklearn.metrics.pairwise import linear_kernel
    cosine_similarities = linear_kernel(tfidf[0:1], tfidf).flatten()
    #cosine_similarities

    #### Construct document similarities ####
    #########################################
    related_docs_indices = cosine_similarities.argsort()[:-(21):-1]
    #print(related_docs_indices)
    #print(cosine_similarities[related_docs_indices])
    # fincorpus[0] #### Print out the contents of the corpus match

    gene = []
    for m in related_docs_indices[1:20]:
        gene.append(re.findall("\n.+ ", fincorpus[m])[0].replace('\n', '').replace(' ', ''))
    return(pd.DataFrame({'Gene_Name' : gene, 'cosine_similarity' : cosine_similarities[related_docs_indices][1:20]}))



