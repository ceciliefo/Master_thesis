import numpy as np
import os
import sys, random
import argparse
from quantisation.kmean_quantisation import KMeansQuantisation
from secure_systems.TripleHashSystem import TripleHashIdentificationSystem
from pathlib import Path
import time


parser = argparse.ArgumentParser(description='Kmeans-based hash quantisation',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-e', '--embeddings', type=str,
                     default='/Users/cecilie/Documents/Code/identification/Features/frgc_arcface512_float/frgc_arcface512_npy')

parser.add_argument('-cl', '--cluster_centers', type=str,
                     default='/Users/cecilie/Documents/Code/identification/Features/frgc_arcface512_float/workload_k_means/frgc_cluster_centers')

parser.add_argument('-o', '--output', type=str,
                     default='/Users/cecilie/Documents/Code/identification/Features/frgc_arcface512_float/workload_k_means')

parser.add_argument('-m', '--measure', type=str,
                     default='/Users/cecilie/Documents/Code/identification/Features/frgc_arcface512_float/workload_k_means/executions')

parser.add_argument('-n', '--name', type=str,
                     default='new_resnet-100')

parser.add_argument('-k', '--k-fold', type=int,
                     default=1) 

parser.add_argument('-c', '--centers', type=int,
                     default=64)

parser.add_argument('-s', '--sub-spaces', type=int,
                     default=1)

parser.add_argument('-t', '--soft', type=int,
                     default=1)

args = parser.parse_args()


embeddings_path = list(Path(args.embeddings).glob('*.npy'))

embeddings_path.sort()

subjects = {}


for p in embeddings_path:

    l = p.stem.split('d')[0]

    if l in subjects:

        subjects[l].append(p)

    else:

        subjects[l] = [p]

average_comparisons = 0

total_false_negative = 0

hit_rate = 0

total_average_entities = 0

ave_rank = 0

fpath_txt = os.path.join(args.output, "results_{}_{}_{}_{}.txt".format(args.name, args.centers, args.sub_spaces, args.soft))

with open(fpath_txt, 'w') as f: 

    # k_means_times = [] 

    for i in range(args.k_fold):

        f.write('Iteration {}\n'.format(i))

        # start_kmean = time.time() 

        total_entities = 0

        search_l = []

        enrol_l = []

        enrol_f = []

        search_f = []

        train_f = []

        train_l = []

        samples_f = [] 

        samples_l = [] 

        search_samples_list = [] 
 
        for s in subjects:

            samples = subjects[s] 

            # random.shuffle(samples)

            enrol_samples = samples[:1] 

            # train_samples = samples[1:] 

            search_samples = samples[1:] 

            search_samples_list.append(len(search_samples))

            # for e in train_samples:

                # train_f.append(np.loadtxt(str(e))) 

                # train_l.append(s) 

            for e in enrol_samples:

                enrol_f.append(np.loadtxt(str(e))) 

                enrol_l.append(s) 
        
            for e in search_samples:

                search_f.append(np.loadtxt(str(e))) 

                search_l.append(s) 
 
            for e in samples:  

                samples_f.append(np.loadtxt(str(e))) 

                samples_l.append(s) 
        
        
        samples_f = np.asarray(samples_f)

        enrol_f = np.asarray(enrol_f)

        count_enrol = enrol_f.shape[0]

        search_f = np.asarray(search_f)
        
        count_search = search_f.shape[0]

        f.write('Count of enrol {}\n, Count of search {}\n'.format(count_enrol, count_search))

        print('building hash for K-fold {}'.format(i))

        model = KMeansQuantisation(K = args.centers, sub_spaces = args.sub_spaces) 

        model.train_model(enrol_f) # train the model using list enrol_f (runs KMeans approach and adds the result into a codebook / cluster list)

        print('Encodes features for K-fold {}'.format(i))

        enroll_codes = model.encode(enrol_f) # return the codes matrix (each row has a 1 on the position to the lowest distance to its closest cluster --> corresponds to the value of the nearest cluster)

        hash_pos_enroll = np.asarray([np.where(enroll_codes[i, :])[0] for i in range(len(enrol_l))]) 
        # --> one index per file in the returning list (with the index of the lowest distance)

        hash_codes_enroll = np.asarray([list(map(lambda e: e % args.centers, hash_pos_enroll[i, :])) for i in range(len(enrol_l))])


        # end_kmean = time.time() 

        # elapsed_time_kmean = end_kmean - start_kmean

        # k_means_times.append(elapsed_time_kmean)

        # print("Time performing k-means on enrollment features: {}".format(elapsed_time_kmean))


        search_codes = model.encode(search_f) 

        hash_pos_search = np.asarray([np.where(search_codes[i, :])[0] for i in range(len(search_l))]) 

        hash_codes_search = np.asarray([list(map(lambda e: e % args.centers, hash_pos_search[i, :])) for i in range(len(search_l))])

        # Compute the number of false negatives: compare hash_codes_enroll and hash_codes_search
        new_hash_codes_enroll = []
        new_hash_codes_search = []
        for i in range(len(search_samples_list)): 
            for j in range(search_samples_list[i]): 
                new_hash_codes_enroll.append(hash_codes_enroll[i][0])

        for i in range(len(hash_codes_search)):
            new_hash_codes_search.append(hash_codes_search[i][0])

        true_match = 0
        false_match = 0
        for i in range(len(new_hash_codes_enroll)):
            if (new_hash_codes_enroll[i] == new_hash_codes_search[i]):
                true_match +=1
            else:
                false_match +=1
        
        print("Assigned to same cluster id: {}".format(true_match))
        print("Assigned to different cluster id: {}".format(false_match))


        model_identification = TripleHashIdentificationSystem() # initialize hash table (dictionary) and dataset (list)

        print('Enrol features for K-fold {}'.format(i))

        model_identification.enrol(hash_codes_enroll, enrol_f, enrol_l) # enrol to hash table


        # Link subject id to cluster id in a new file
        subjects_and_hashes = model_identification.sujects_to_stable_hash 

        subjects_and_hashes_file_txt = os.path.join(args.output, "frgc_results_subject_id_to_cluster.txt")

        with open(subjects_and_hashes_file_txt, 'w') as file: 

            file.write('Link stable hashes and subjects for iteration {} \n'.format(i)) 
            
            for key, value in subjects_and_hashes.items(): 
                file.write("{}:{} \n".format(key, value))


        # Link cluster id to subjects in a new file
        cluster_and_subjects = model_identification.cluster_to_subjects
        new_cluster_and_subjects = {int (key) : (value) for key, value in cluster_and_subjects.items()}

        cluster_id_to_subjects_file_txt = os.path.join(args.output, "frgc_results_cluster_id_to_subjects.txt")

        with open(cluster_id_to_subjects_file_txt, 'w') as file: 

            file.write('Link cluster key (from cluster 0 to cluster 63) to their subjects \n\n') 

            for key in sorted(new_cluster_and_subjects):
                file.write("{} : {} \n".format(key, new_cluster_and_subjects[key]))
            
        # Extract the centers 
        kmeans_centers = model.C[0].cluster_centers_

        # Make each cluster a new file, where each value is on a new line
        for i in range(len(kmeans_centers)): 

            centers_file_txt = os.path.join(args.cluster_centers, "frgc_cluster_{}.npy".format(i))

            with open(centers_file_txt, "w") as file: 

                for value in kmeans_centers[i]: 

                    file.write("{}\n".format(value))

        # Copy the cluster centers to a new file 
        kmeans_centers_file_txt = os.path.join(args.output, "frgc_results_kmeans_centers")

        with open(kmeans_centers_file_txt, 'w') as file: 

            for i in range(args.k_fold):

                file.write('Kmeans centers for iteration {} \n'.format(i)) 

                file.write(str(kmeans_centers))

        # Create a file of the cluster for each subject
        clustes_file_txt = os.path.join(args.output, "frgc_clusters.txt")

        with open(clustes_file_txt, "w") as file: 

            for i in range(len(model_identification.sujects_to_stable_hash.keys())):

                file.write('{}\n'.format(list(model_identification.sujects_to_stable_hash.values())[i][0]))

        # Compute average comparisons, candidate list, ...
        average_comparisons += model_identification.counting_samples() # return average as the count feet (number of collissions) divided by numbers of keys in the hash table

        total_entities = model_identification.getting_keys_hash_table() # returning the number of keys in the hash table (i.e. the size)
        
        f.write('Total entities: {} \n'.format(total_entities))

        total_average_entities += total_entities

        candidate_list = model_identification.search(hash_codes_search, search_f) # return the result (list consisting of tuple of distance and subject identifier)
 

        # Save mated and non-mated scores
        mated_scores_file_txt_v3 = os.path.join(args.output, "frgc_mated_scores_v3")

        non_mated_scores_file_txt_v3 = os.path.join(args.output, "frgc_non_mated_scores_v3")

        m,n= model_identification.func_score_vales_v3(embeddings_path, mated_scores_file_txt_v3, non_mated_scores_file_txt_v3, samples_f, samples_l)


        print('Check candidate list for K-fold {}'.format(i))

        false_negative = 0 

        rank_1 = 0

        for real, cand in zip(search_l, candidate_list): 
            if not real in map(lambda e: e[1], cand): # check if the original subject identifier (real) is not in the candidate liste (cand)

                false_negative += 1 
            
            elif len(cand) > 0 and cand[0][1] == real: # check if the subject in the candidate list corresponds to the original subject

                rank_1 += 1 

        total_false_negative += false_negative 

        hit_rate += (len(search_l) - false_negative)/len(search_l) # proportion of subjects where corresponding subject is in pre-selection

        ave_rank += rank_1/(len(search_l) - false_negative) 

        f.write('Total non-found: {}, accuracy:{}\n'.format(false_negative, (len(search_l) - false_negative)/len(search_l))) # compute accuracy as amount of subject minus false negative divided by amount of subjects (one computation of hit_rate)
        # compute accuracy as amount of subject minus false negative divided by amount of subjects


    average_comparisons /= args.k_fold 

    total_false_negative /= args.k_fold 

    ave_rank /= args.k_fold

    total_average_entities /= args.k_fold 

    f.write('average number of comparisons: {}\n'.format(average_comparisons))

    f.write('Total of average recognition rate for rank 1: {}\n'.format(ave_rank))

    # total_non_found = false_negative / args.k_fold    

    total_hit_rate = hit_rate / args.k_fold 

    f.write('average hit rate: {} \naverage non found: {}\n'.format(total_hit_rate, total_false_negative))

    f.write('Total average of entities: {} \n'.format(total_average_entities))

    print('...done')


    # k_means_times.sort() # Test measuring time for kmeans on enrollment
    # middle_kmeans = len(k_means_times) // 2
    # median_kmeans = (k_means_times[middle_kmeans] + k_means_times[~middle_kmeans])/2

    # print("Median time for kmeans procedure for 10 executions on enrollment references is: {}".format(median_kmeans))


# Measure the time for finding the corresponding stable hash, i.e., the cluster index of a new feature 
print("Measure one stable hash computations after the enrollment phase is finished")

number_of_executions = 10

probe_samples = []

for key in subjects:
    probe_samples += subjects[key][1:] # Remove the first sample for each subject (used for enrolment)

print(len(probe_samples))

probe_feature_measure_file = os.path.join(args.measure, "frgc_results_measure_probe_stable_hash_for_{}_executions".format(number_of_executions)) 

with open(probe_feature_measure_file, 'w') as file: 

    file.write("Measuring the time for computing the stable hash of a probe feature \n")

    file.write("Running {} executions for random probe features\n".format(number_of_executions)) 

    file.write("Running {} executions of cluster computation and finding the median time".format(number_of_executions)) 

    start_execution = time.time()

    for i in range (number_of_executions): 

        probe_number = random.randint(0, (len(probe_samples)-1))

        probe = probe_samples[probe_number] 

        file.write("\n\nIteration: {}\n".format(i)) 
        file.write("Probe feature: {}\n".format(probe.stem.split('.')[0])) 

        probe_cluster = []

        file.write("[") 

        times = []

        for j in range(number_of_executions):

            probe_comparison = {}

            start_time = time.time()

            probe_feature = np.loadtxt(str(probe)) # Load the content

            for center in range(len(kmeans_centers)): 
                sum_sq = np.sum(np.square(probe_feature - kmeans_centers[center]))
                probe_comparison[center] = sum_sq

            sorted_probe_comparison = sorted(probe_comparison.items(), key = lambda kv: (kv[1], kv[0]))

            probe_cluster.append(sorted_probe_comparison[0][0])

            end_time = time.time()

            elapsed_time = end_time - start_time

            times.append(elapsed_time)

            print("Execution time of computating a stable hash value for a new probe is:", elapsed_time, "seconds.")

            file.write("Execution time for iteration {} : {} seconds\n".format(j, elapsed_time)) 

            file.write("Cluster id: {}\n".format(sorted_probe_comparison[0][0])) 

            file.write("{}, {}, {}\n".format(probe.stem.split('.')[0], elapsed_time, sorted_probe_comparison[0][0])) 
    
        times.sort()
        middle = len(times) // 2
        median = (times[middle] + times[~middle])/2

        file.write("probe:{}, median:{}, cluster:{}]\n\n".format(probe.stem.split('.')[0], median, np.unique(probe_cluster))) # for compromised version

    end_execution = time.time()

    elapsed_execution = end_execution - start_execution

    print("Total execution time for {} executions: {}".format(number_of_executions, elapsed_execution))

