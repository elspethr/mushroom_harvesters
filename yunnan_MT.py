### set up ###

import sys, os

#still don't know why I need this :/
sys.path.append('/Users/elspethready/repos/multitensor/build')

import multitensor
import pandas as pd
import networkx as nx
import os
import numpy as np
import tools as tl
import cv_functions as cvfun
import importlib


### some functions for analysis ###

def prep_test_train(data0, ly_idx, nodeId2Name, cv_type='kfold', NFold=5, rseed=10, fold=0, out_mask=False):
    data = tl.filter_layers_from_sktensor(data0, ly_idx)
    nnzB_subs = np.where(data.subs[0] == 0)[0] # account only for the non-zero entries of the 0-th layer (test layer)
    nnzB = len(nnzB_subs)
    if nnzB > 7500: # arbitrary threshold
        samples_test = 7500 // NFold
        samples_train = 7500 - samples_test
    else:
        samples_test = nnzB // NFold
        samples_train = nnzB - samples_test
    if cv_type == 'kfold':
        idxG = cvfun.shuffle_indicesG(nnzB_subs, rseed=rseed)
    else:
        idxG = nnzB_subs.copy()
    data_train, testG = cvfun.extract_masks_layer_interdependence(data, idxG=idxG, cv_type=cv_type, NFold=NFold,
                                                                  fold=fold, rseed=rseed, out_mask=out_mask)
    df, data1, adjacency_file = tl.write_data_to_file(data_train, nodeId2Name = nodeId2Name)
    return adjacency_file, data, idxG, testG, samples_test, samples_train

def run_mt_with_cv(adjacency_file, data, idxG, testG, K, samples_train, samples_test, directed, rseed=10, assortative=True):        
   
   u0, v0, w0, report = multitensor.run(adjacency_file, K, seed=rseed, nof_realizations=10,
                                        assortative=assortative, max_nof_iterations=500,
                                        directed=directed, nof_convergences=2)
   u = tl.sort_membership_by_nodeID(u0, method='asint')
   if directed == True:
       v = tl.sort_membership_by_nodeID(v0, method='asint')
   else:
       v = u.copy()
   w = np.array(w0)
   auc_train = cvfun.calculate_AUC(data, u, v, w, mask=idxG[np.logical_not(testG)], verbose=0,
                                   samples=samples_train, n_comparisons=1e5, rseed=rseed)
   auc_test = cvfun.calculate_AUC(data, u, v, w, mask=idxG[testG], samples=samples_test,
                                  verbose=0, n_comparisons=1e5, rseed=rseed)
   
   return u, w, auc_train, auc_test

def normalize_membership(u):
  row_sum = u.sum(axis=1)
  u_norm = np.copy(u)
  u_norm[row_sum>0] = u_norm[row_sum>0]/row_sum[row_sum>0,np.newaxis]
  return u_norm

def binary_assignment_communities(u,threshold=0.1):
  u_norm = normalize_membership(u)
  return [np.where(u_norm[i]>threshold)[0] for i in range(u.shape[0])]


### main ###

def main():
  ### set up variables ###

  #network info
  undirected = True
  in_folder = '../yunnan_data/input/'
  out_folder = '../yunnan_data/output/'
  years = ["2014", "2015"]
  assortative = True

  #mt and cv params
  rseed = 78953
  directed = np.logical_not(undirected)

  #0 sup
  #1 social
  #2 info
  #3 labor
  #4 kin
  #5 clan
  #6 dist
  #7 harv

  #what layer/group combinations to run
  ly_combs = [[0], [7], [3],
              [5, 6], [4, 6],
              # original sup models
              [0, 4], [0, 5], [0, 6], [0, 4, 6], [0, 5, 6], [0, 4, 5], [0, 4, 5, 6],
              #[0, 5], [0, 6], [0, 7],
              # new sup models
              [3, 4], [3, 5], [3, 6], [3, 4, 6], [3, 5, 6],
              [3, 4, 5], [3, 4, 5, 6],
              [3, 1], [3, 2], [3, 1, 2],
              # original harv models
              [7, 4], [7, 5], [7, 6], [7, 4, 5, 6], [7, 0], [7, 0, 4, 5, 6],
              [7, 4, 6], [7, 4, 5], [7, 4, 0], [7, 0, 4, 6], [7, 4, 5, 0],
              [7, 0, 5], [7, 0, 6], [7, 5, 6], [7, 0, 5, 6],
              # new harv models
              [7, 1], [7, 2], [7,3],
              [7, 1, 2], [7, 1, 3], [7, 2, 3], [7, 1, 2, 3],
              # new harv models with dist
              [7, 1, 6], [7, 2, 6], [7, 3, 6],
              [7, 1, 2, 6], [7, 1, 3, 6], [7, 2, 3, 6],
              [7, 1, 2, 3, 6],
              # new harv models with kin
              [7, 1, 4], [7, 2, 4], [7, 3, 4],
              [7, 1, 2, 4], [7, 1, 3, 4], [7, 2, 3, 4],
              [7, 1, 2, 3, 4],
              # new harv models with kin and dist
              [7, 1, 4, 6], [7, 2, 4, 6], [7, 3, 4, 6],
              [7, 1, 2, 4, 6], [7, 1, 3, 4, 6], [7, 2, 3, 4, 6], [7, 1, 2, 3, 4, 6]
              ]
  comb_name = ["sup", "harv", "labor",
               "dist_clan", "dist_kin",
               # original sup models
               "sup_kin", "sup_clan", "sup_dist", "sup_kindist", "sup_clandist", "sup_kinclan", "sup_kinclandist",
               #"sup_kin125", "sup_kin25", "sup_kin5",
               # new sup models
               "labor_kin", "labor_clan", "labor_dist", "labor_kindist", "labor_clandist",
               "labor_kinclan", "labor_kinclandist",
               "labor_social", "labor_info", "labor_socialinfo",
               # original harv models
               "harv_kin", "harv_clan", "harv_dist", "harv_kinclandist", "harv_sup", "harv_all",
               #"harv_kin125", "harv_kin25", "harv_kin5",
               "harv_kindist", "harv_kinclan", "harv_kinsup", "harv_kinsupdist", "harv_kinclansup",
               "harv_supclan",  "harv_supdist", "harv_clandist", "harv_supclandist",
               # new harv models
               "harv_social", "harv_info", "harv_labor",
               "harv_socialinfo", "harv_sociallabor", "harv_infolabor", "harv_socialinfolabor",
               # new harv models with dist
               "harv_socialdist", "harv_infodist", "harv_labordist",
               "harv_socialinfodist", "harv_sociallabordist", "harv_infolabordist",
               "harv_socialinfolabordist",
               # new harv models with kin
               "harv_socialkin", "harv_infokin", "harv_laborkin",
               "harv_socialinfokin", "harv_sociallaborkin", "harv_infolaborkin",
               "harv_socialinfolaborkin",
               # new harv models with kin and dist
               "harv_socialkindist", "harv_infokindist", "harv_laborkindist",
               "harv_socialinfokindist", "harv_sociallaborkindist", "harv_infolaborkindist",
               "harv_socialinfolaborkindist"
                ]
  min_grp = 2
  max_grp = 10

  ### loop for layer interdependence/link prediction ###

 # for year in years:
 #  #arrange data
 #  #label = '_harv'
 #  adj_name = 'yunnan_networks_' + year + '.dat'
 #  A, data0, nodes, nodeName2Id, nodeId2Name = tl.import_data(in_folder, adj_name=adj_name, sep='\s+', header=None,
 #                                                             ego=0, alter=1, undirected=undirected, verbose=0)
 #  for nf in [2,3,4,5]:
 #    n = -1
 #    for ly_idx in ly_combs:
 #      n +=1
 #      print(comb_name[n])
 #      #prep test/train split
 #      adjacency_file, data, idxG, testG, samples_test, samples_train = prep_test_train(data0, ly_idx, nodeId2Name, NFold=nf, rseed=rseed)
 #      #run multitensor for n groups
 #      auc_array = np.empty((0,3), float)
 #      for K in range(min_grp, (max_grp+1)): 
 #        u, w, auc_train, auc_test = run_mt_with_cv(adjacency_file, data, idxG, testG, K, samples_train,
 #                                                     samples_test, rseed=rseed, directed=directed, assortative=assortative)
 #        auc_array = np.append(auc_array, np.array([[K, auc_train, auc_test]]), axis=0)
 #        u_norm = normalize_membership(u)
 #        #cm = binary_assignment_communities(u, threshold=0.1)
 #        fname1 = out_folder + year +"/cv_membership_" + comb_name[n] + "_%s_%sFCV.txt" % (K, nf)
 #        with open(fname1, 'w') as f:
 #          for node in u_norm:
 #            for prob in node:
 #              f.write("%s," % prob)
 #            f.write("\n")
 #        fname2 = out_folder + year +"/cv_auc_" + comb_name[n] + "_%sFCV.txt" % nf
 #        np.savetxt(fname2, auc_array, delimiter=",")

  ### loop for straight up group membership ###

  best_models = ["kin", "sup", "clan", 
                 "supdist", "labor", "labordist", 
                 "socialdist","socialinfolabor"]
  for year in years:
    #arrange data
    for model in best_models:
      model2 = model + '_harv'
      model3 = model + '_harv_first' #for some debugging...
      for mod in [model, model2, model3]:
        #print(mod)
        adj_name = in_folder + 'yunnan_networks_' + year + "_" + mod + '.dat'
        for K in range(min_grp, (max_grp+1)):
          u0, v0, w0, report = multitensor.run(adj_name, K, seed=rseed, nof_realizations=10,
                                               assortative=assortative, max_nof_iterations=500, directed=directed,
                                               nof_convergences=2)
          fname3 = out_folder + year +"/membership" + mod + "_%s.txt" % K
          with open(fname3, 'w') as f:
            for node in u0: #note multitensor.run gives diff output than run_mt_with_cv
              for prob in node:
                f.write("%s," % prob)
              f.write("\n")
          fname4 = out_folder + year + "/affinity" + year + "_" + mod + "_%s.txt" % K
          with open(fname4, 'w') as f:
            for group in w0:
              for aff in group:
                f.write("%s," % aff)
              f.write("\n")    

######

if __name__ == "__main__":
    main()
