#!/usr/bin/env python3
#coding:utf-8

#default libs
from glob import glob
from collections import defaultdict
from math import sqrt
import time 
from statistics import mean

class GI:
    def __init__(self, contig, ID, start, stop):
        self.contig = contig
        self.id = ID
        self.start= int(start)
        self.stop = int(stop)

def get_references(filename):
    ref = defaultdict(list)
    with open(filename,"r") as f:
        for line in f:
            start, stop = line.split(':')[1].split('..')
            link = line.split('.')[0]
            ref[link].append((int(start), int(stop)))
    myref = {}
    for key, val in ref.items():
        myref[key] = sorted(val, key = lambda x : x[0])
    return myref

def sort_predicted(predicted):
    new = {}
    for key, val in predicted.items():
        new[key] = sorted(val, key = lambda x : x.start)
    return new

def get_intersect_pairs(pairs1, pairs2):
    i = 0
    j = 0
    inter = 0
    sorted_pairs1 = sorted(pairs1,key = lambda x : x[0])
    sorted_pairs2 = sorted(pairs2,key = lambda x : x[0])
    while i < len(sorted_pairs1) and j < len(sorted_pairs2):
        curr1 = sorted_pairs1[i]
        curr2 = sorted_pairs2[j]
        if curr1[0] > curr2[1]:
            j+=1
        elif curr2[0] > curr1[1]:
            i+=1
        else:#there is some intersect
            inter += min(curr1[1], curr2[1]) - max(curr1[0], curr2[0])
            if curr1[1] < curr2[1]:
                i+=1
            elif curr1[1] > curr2[1]:
                j+=1
            else:
                i+=1
                j+=1
    return inter

def overlap(pairs):
    over = 0
    for pair in pairs:
        over += pair[1] - pair[0]
    return over

def compute_score(pred, pos ,neg, name):
    # print(f"Computing the scores of {name} compared to positive and negative datasets.")
    TP = 0
    FP = 0
    TN = 0
    FN = 0
    nbkeys= 0
    recalls = []
    precisions = []
    accuracies = []
    F1scores = []
    MCCs = []
    for key in pos.keys():#for each reference genome
        GIsPos = []
        GIs = pos.get(key)
        if GIs is not None:
            for gi in GIs:
                GIsPos.append((gi[0], gi[1]))
        GIsneg = []
        GIs = neg.get(key)
        if GIs is not None:
            for gi in GIs:
                GIsneg.append((gi[0],gi[1]))
        GIs = pred.get(key)
        GIspred = []
        if GIs is not None:
            nbkeys+=1
            for gi in GIs:
                GIspred.append((gi.start, gi.stop))
        interpos = get_intersect_pairs(GIspred, GIsPos)
        TP = interpos
        FN = overlap(GIsPos) - interpos
        interneg = get_intersect_pairs(GIspred, GIsneg)
        FP = interneg
        TN = overlap(GIsneg) - interneg

        recalls.append(TP / (TP + FN))
        if TP + FP == 0:
            precisions.append(1)
        else:
            precisions.append( TP / (TP + FP))
        accuracies.append((TP + TN) / (TP + FN + FP + TN))
        F1scores.append( 2*TP / (2*TP + FN + FP))
        if (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) == 0:
            MCCs.append(0)
        else:
            MCCs.append((( TP * TN ) - (FP * FN )) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
    results = {"Recall" : round(mean(recalls),3), "Precision": round( mean(precisions),3), "Accuracy": round( mean(accuracies),3), "F1 score":round( mean(F1scores),3) , "MCC": round(mean(MCCs),3), "name":name}
    return results

def link_GI_to_predicted(filename, predicted):
    curr_GI = defaultdict(list)
    c = 0
    with open(filename,"r") as f:
        f.readline()
        prec = ""
        for i in f:
            line = i.split()
            if not line[0] == prec:
                link = None
                if any(gbk in line[0] for gbk in predicted.keys()):
                    link = line[0].split('.')[0]
            prec = line[0]
            if link is not None:
                curr_GI[link].append(GI(line[0],str(c), line[1], line[2]))
    return sort_predicted(curr_GI)

def print_scores(*args):
    args = sorted(args, key = lambda x : x["MCC"], reverse=True)
    print("Tool\tMCC\tF1 score\tAccuracy\tPrecision\tRecall")
    for arg in args:
        print(f"{arg['name']}\t{arg['MCC']}\t{arg['F1 score']}\t{arg['Accuracy']}\t{arg['Precision']}\t{arg['Recall']}")


def run_dataset(dataset, neg, comment):
    set_to_use = dataset
    print(comment)
    refseqpanRGP = link_GI_to_predicted('dataset_panRGP.txt',set_to_use)

    keys_to_use = { key: val for key, val in refseqpanRGP.items() if key in set_to_use }

    new_dataset = {}
    for key in dataset:
        if key in keys_to_use:
            new_dataset[key] = dataset[key]
    set_to_use = new_dataset

    new_neg = {}
    for key in neg:
        if key in keys_to_use:
            new_neg[key] = neg[key]
    negative = new_neg


    test_use = {}
    for key, val in set_to_use.items():
        test_use[key] = []
        for x in val:
            test_use[key].append(GI(None, None, x[0], x[1]) )

    islandviewer = link_GI_to_predicted("all_gis_islandviewer_iv4.txt", keys_to_use)
    alienHunter = link_GI_to_predicted("all_gis_alienHunter.txt", keys_to_use)
    SigiHMM = link_GI_to_predicted("all_gis_SigiHMM.txt", keys_to_use)
    SigiCFR = link_GI_to_predicted("all_gis_SigiCFR.txt", keys_to_use)
    islandpathdimob = link_GI_to_predicted("all_gis_dimob.txt", keys_to_use)
    predictBias = link_GI_to_predicted("all_gis_PredictBias.txt", keys_to_use)
    zislandExplorer = link_GI_to_predicted("all_gis_ZislandExplorer.txt", keys_to_use)
    IslandCafe = link_GI_to_predicted("all_gis_IslandCafe.txt", keys_to_use)
    xenoGI = link_GI_to_predicted("all_gis_xenoGI.txt",keys_to_use)
    GI_cluster = link_GI_to_predicted("all_gis_GI_cluster.txt", keys_to_use)

    GI_clusterres = compute_score(GI_cluster, set_to_use, negative, "GI-Cluster")
    islandviewerres = compute_score(islandviewer, set_to_use, negative, "islandviewer4")
    alienhunterres = compute_score(alienHunter, set_to_use, negative, "AlienHunter")
    sigiHMMres = compute_score(SigiHMM, set_to_use, negative, "SigiHMM")
    SigiCFRres = compute_score(SigiCFR, set_to_use, negative, "SigiCRF")
    islandpathdimobres = compute_score(islandpathdimob, set_to_use, negative, "IslandPath-DIMOB")
    predictBiasres = compute_score(predictBias, set_to_use, negative, "PredictBias")
    zislandExplorerres = compute_score(zislandExplorer, set_to_use, negative, "ZislandExplorer")
    IslandCaferes = compute_score(IslandCafe, set_to_use, negative, "IslandCafe")
    xenoGIres = compute_score(xenoGI, set_to_use, negative, "xenoGI")
    refpanres = compute_score(refseqpanRGP, set_to_use, negative, 'panRGP-refseq')

    print_scores(refpanres, islandviewerres, alienhunterres, GI_clusterres, sigiHMMres, SigiCFRres, islandpathdimobres, predictBiasres, zislandExplorerres, IslandCaferes, xenoGIres)


def main():

    #Expects many files to be existing in the working directory.
    #if one of them is missing, errors will be raised.

    positive = get_references("positive_dataset.txt")
    negative = get_references("negative_dataset.txt")
    literature = get_references("literature_dataset.txt")

    for dset, comment in [(literature, "L-dataset"),(positive,"C-dataset")]:
        run_dataset(dset, negative, comment)

if __name__=="__main__":
    main()
