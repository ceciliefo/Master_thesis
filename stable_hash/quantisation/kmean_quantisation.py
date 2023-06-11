import numpy as np
from sklearn.cluster import KMeans
import pickle, os

class KMeansQuantisation:

    def __init__(self, K = 256, sub_spaces = 4, soft_assignment = 1):

        self.C = []

        self.K = K

        self.sub_spaces = sub_spaces 

        if soft_assignment >= K:

            raise Exception('ERROR: assignment = {} is greater than the number of centers {}'.format(soft_assignment, K))

        self.soft_assignment = soft_assignment

    
    def train_model(self, features):

        n, m = features.shape # rows and columns of the features list (e.g., enrol_f)

        if m % self.sub_spaces != 0: # feature dimension (the number of colums) must a multiple of the sub_spaces (P)

            print('feature dimension {} must be multiple of {}'.format(m, self.sub_spaces))

            return None

        num_sub_spaces = int(m/self.sub_spaces) 

        for i in range(self.sub_spaces): 

            sub_features = features[:, i*num_sub_spaces:(i + 1)*num_sub_spaces] 
    
            kmeans = KMeans(n_clusters=self.K, random_state=0).fit(sub_features) # iterative divides data points into K clusters (e.g. 64 clusters) by minimizing variance in each cluster 

            self.C.append(kmeans) # append the result to the list of clusters / codebook  


    def __predict(self, features, centers):
        
        distance = np.asarray([np.linalg.norm(centers - f, axis=1) for f in features]) # new array of the distance between the centers and the sub_features list 

        result = []

        for dist in distance:
            index = np.arange(0, self.K) 

            merge = list(zip(dist, index))  

            merge.sort(key=lambda tup: tup[0]) # sort the list from lowest distance (sort it based on the first argument (dist))
            
            merge = merge[:self.soft_assignment] #when soft_assignment = 1 --> take out the first element --> extract the one with the lowest distance
            
            result.append([l for d,l in merge]) # append the the lowest value in the result list --> only append the index value (which K value from the index list)
 
        return np.asarray(result) # return the list of indexes of the minimal distances after going through all the distances 

    
    def encode(self, features):

        n, m = features.shape 

        if m % self.sub_spaces != 0:

            print('feature dimension {} must be multiple of {}'.format(m, self.sub_spaces))

            return None

        num_sub_spaces = int(m/self.sub_spaces)

        codes = np.zeros((n, self.K*self.sub_spaces)) # create an array of zeros with shape of n rows (rows from features list / amount of files) and K*p columns (elements in each row)
       
        for i in range(self.sub_spaces):

            sub_features = features[:, i*num_sub_spaces:(i + 1)*num_sub_spaces] 

            kmeans = self.C[i] 

            top_nearest = self.__predict(sub_features, kmeans.cluster_centers_) # run predict function and return list of the index values with the lowest distance

            top_nearest += i*self.K 

            # nearest = kmeans.predict(sub_features)

            # nearest = nearest + i*self.K

            for j in range(n):

                codes[j, top_nearest[j, :]] = 1 
                # change out one zero in the zero matrix with a 1 for the value in the row j and on the colum value from index of the nearest element from top_nearest

        return codes 

        
    def save_model(self, output_path):

        with open(os.path.join(output_path, 'centers.pkl'), 'wb') as f:
            pickle.dump(self.C, f)

    def load_model(self, input_file):

        with open(os.path.join(input_file, 'centers.pkl'), 'rb') as f:
            self.C = pickle.load(f)



