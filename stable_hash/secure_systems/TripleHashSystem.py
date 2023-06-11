import numpy as np
from scipy.spatial import distance
import time



class TripleHashIdentificationSystem:

    def __init__(self):

        self.hash_table = {} 

        self.dataset = []

        self.sujects_to_stable_hash = {} # key is subject, value is cluster id

        self.cluster_to_subjects = {} # key is cluster id, value is a list of subjects id

    
    def enrol(self, g_hash, g_features, labels): # enrol to hash table 

        for h, f, l in zip(g_hash, g_features, labels): 

            h_tmp = ''.join(str(e) for e in h) 

            if h_tmp in self.hash_table: # if the index already exist 

                self.hash_table[h_tmp].append((f, l))
            
            else:

                self.hash_table[h_tmp] = [(f, l)] # initialize from start --> value is a list corresponding of two values (file values and subject identifier)

            # Initalize subjects to stable hash dictionary and cluster to subject dictionary
            if l in self.sujects_to_stable_hash:
                if not (h_tmp in self.sujects_to_stable_hash[l]):
                    self.sujects_to_stable_hash[l].append(h_tmp)
            else: 
                self.sujects_to_stable_hash[l] = [h_tmp]

            if h_tmp in self.cluster_to_subjects: 
               if not (l in self.cluster_to_subjects[h_tmp]):
                   self.cluster_to_subjects[h_tmp].append(l)
            else: 
               self.cluster_to_subjects[h_tmp] = [l]
            

    def enrol_exhaustive(self, g_features, labels):

        for f, l in zip(g_features, labels):

            self.dataset.append((f, l))


    def search(self, q_hash, q_features):

        result = []

        for h, f in zip(q_hash, q_features): 

            h_tmp = ''.join(str(e) for e in h)

            if h_tmp in self.hash_table: # if index already exists as a key in hash table

                entry =  self.hash_table[h_tmp] # extract out the values for the already defined key (index)

                dist = [] 

                for e, l in entry: # differentiate the two values for the value entry of the key

                    value = distance.sqeuclidean(f, e)

                    dist.append((value, l)) # append the distance in the new list as a tupple corresponding with the distance value and the subject

                dist.sort(key=lambda tup: tup[0]) # sort the list of the distances according to the first value (the distance)

                result.append(dist) 
                
            else:

                result.append([]) #if h_tmp not exist as a key in hash table, append empty list to result list         
    
        return result # return the result (list consisting of tuple of distance and subject identifier)
    
    
    def search_exhaustive(self, q_features):

        result = []

        for q in q_features:

            dist = []

            for element_e in self.dataset:

                feat_e = element_e[0]

                label_feat_e = element_e[1]

                value = distance.sqeuclidean(q, feat_e)

                dist.append((value, label_feat_e))
            
            # dist.sort(key=lambda tup: tup[0])

            result.append(dist)
        
        return result


    def find_label_lost(self, label):

        for g_hash in self.hash_table:

            elements_colissions = self.hash_table[g_hash]

            for value in elements_colissions:

                v = value[0]

                l = value[1]

                if label == value[1]:

                    found = True

                    print ("Found in hashtable {}".format(str(g_hash)))

                    list_elements = self.hash_table[g_hash]

                    for e in list_elements:
                        
                        print ("Colissioned candidates for {}".format(str(e[1])))


    def counting_samples(self):

        average = 0

        count_entry = 0

        count_feat = 0

        for g_hash in self.hash_table:

            count_entry +=1 

            elements_colissions = self.hash_table[g_hash] 

            count_feat += len(elements_colissions)
        
        average = count_feat / count_entry 

        return average 
    

    def getting_keys_hash_table(self):

        count_entities = 0

        count_entities = len(self.hash_table)

        return count_entities 
    
    # Function to store mated and non-mated scores (only look at one sample per subject to find non-mated samples)
    def func_score_vales_v3(self, embeddings_path, mated_file, non_mated_file, samples, samples_id):
        mated = 0
        non_mated = 0
        first = True
        with open(mated_file, "w") as mated_f, open(non_mated_file, "w") as non_mated_f:
            for i in range(len(samples)):

                if (i>0):
                    if (embeddings_path[i].stem.split('d')[0] == embeddings_path[i-1].stem.split('d')[0]): # frgc
                    # if (embeddings_path[i].stem.split('_')[0] == embeddings_path[i-1].stem.split('_')[0]): # feret
                        first = False
                    else:
                        first = True

                for j in range(len(samples)):

                    if (embeddings_path[i] != embeddings_path[j]):

                        value = distance.euclidean(samples[i], samples[j])

                        if (samples_id[i] == samples_id[j]):
                            mated_f.write('{}\n'.format(value))
                            mated += 1
                    
                        else:
                            if first: 
                                non_mated_f.write('{}\n'.format(value))
                                non_mated += 1
    
        return (mated, non_mated)
    
    # Function to store mated and non-mated scores 
    def func_score_vales_v2(self, embeddings_path, mated_file, non_mated_file, samples, samples_id):
        mated = 0
        non_mated = 0

        with open(mated_file, "w") as mated_f, open(non_mated_file, "w") as non_mated_f:

            for i in range(len(samples)):

                for j in range(len(samples)):

                    if (embeddings_path[i] != embeddings_path[j]):

                        value = distance.euclidean(samples[i], samples[j])

                        if (samples_id[i] == samples_id[j]):
                            mated_f.write('{}\n'.format(value))
                            mated += 1
                    
                        else:
                            non_mated_f.write('{}\n'.format(value))
                            non_mated += 1
        
        return (mated, non_mated)







                



