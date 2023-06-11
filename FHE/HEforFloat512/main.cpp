// @file  simple-integers.cpp - Simple example for BFVrns (integer arithmetic).
// @author TPOC: contact@palisade-crypto.org
//
// @copyright Copyright (c) 2019, New Jersey Institute of Technology (NJIT))
// All rights reserved.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution. THIS SOFTWARE IS
// PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
// EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "openfhe.h" 
#include "biometric.h"
#include "comparison.h"
#include <fstream>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <stdio.h>
#include "dirent.h"


using namespace lbcrypto;
using namespace std;
using Element = DCRTPoly;

Ciphertext<Element> encDB[nrSubj] = {};
string samples[nrSamp + 4] = {};
string subjects[nrSubj] = {};
int originalIDs[nrSubj] = {};
double plainProbeIden[templateSize] = {};
double plainProbeIdenBatched[8*templateSize] = {};
double tempArr[templateSize] = {};
double tempArr2[templateSize] = {};
double cleanUpArray[templateSize] = {};

int candidateList[nrSubj] = {};

#include "Test.cpp"


int main(int argc,char *argv[]) {

    //TestComparison(argv);
    
    //argv[2] = 124, argv[3] = "940128_fb" is FERET sample 208_940128_fb.txt (where 208_940128_fa.txt is enrolled)
    //bin512_identification ../candidatelist.txt 124 940128_fb
    //argv[2] = 903, argv[3] = "d38" is FRGC sample 4763d38.txt (where 04763d02.txt is enrolled)
    //bin512_identification ../candidatelist.txt 903 d38
    //
    //int512_identification ../candidatelist3.txt 124 940128_fb
    //int512_identification ../candidatelist4.txt 0 940928_rc
     
    //
    //Aquire candidate list
    //
    
    int nrCandidates = readCandidateList(candidateList, argc, argv);
    cout << "Number of candidates = " << nrCandidates << endl;
    
    
    //
    //Key generation
    //
    auto start = std::chrono::steady_clock::now();
    
    /*
    uint32_t multDepth = 1;
    uint32_t scaleFactorBits = 50;
    uint32_t batchSize = 4096;
    SecurityLevel securityLevel = HEStd_128_classic;
    
    
    CryptoContext<DCRTPoly> cc =
      CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(
          multDepth, scaleFactorBits, batchSize, securityLevel);
    
    cc->Enable(ENCRYPTION); 
    cc->Enable(SHE);
    
    KeyPair<Element> keyPair = cc->KeyGen();
    */
    
    CCParams<CryptoContextCKKSRNS> parameters; 
        
    uint32_t multDepth = 1;
    uint32_t scaleFactorBits = 30;
    uint32_t batchSize = 4096;
    SecurityLevel securityLevel = HEStd_128_classic;
    
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(scaleFactorBits);
    parameters.SetBatchSize(batchSize);
    parameters.SetSecurityLevel(securityLevel);
    
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    
    //cc->Enable(ENCRYPTION); 
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(PKE);
    //cc->Enable(SHE);
    
    KeyPair<Element> keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    cc->EvalAtIndexKeyGen(keyPair.secretKey, {1});
    
    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    cout << "Key Generation took " << elapsed.count() << " milliseconds" << endl;
    
    //
    // Database Setup
    //
    
    start = std::chrono::steady_clock::now();
    
    setupEncDB(subjects, samples, encDB, originalIDs, cc, keyPair);
    //setupEncDBwBatching(subjects, samples, encDB, originalIDs, cc, keyPair); //used for executions

    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    cout << "Encrypted Database Setup (Combined FERET + FRGC) of " << nrSubj << " subjects took " << elapsed.count() << " milliseconds" << endl;
      
    //
    //k-batching database setup
    //
    
    start = std::chrono::steady_clock::now();
    
    int searchForSubject = atoi(argv[2]);
    int x = 8;
    vector<double>B;
    Ciphertext<Element> encProbeIdenBatched;
    
    for(int j = 0; j < 1062; j++){
    readTemplate(searchForSubject, argv[3], plainProbeIden, originalIDs);
    }
    while(x--){          
       for(int i=0;i<templateSize;i++){
           B.push_back(plainProbeIden[i]);
       }
    }
    copy(B.begin(), B.end(), plainProbeIdenBatched);

    for(int j = 0; j < 34; j++){
    encryptBatchedTemplate(encProbeIdenBatched, plainProbeIdenBatched, cc, keyPair);
    }
    
    //cout << sizeof(plainProbeIdenBatched) << endl;
    
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    cout << "Batched Database Setup for 133 ciphertexts took " << elapsed.count() << " milliseconds" << endl;
    
    int nb_executions = 1000, median, sum = 0; 
    //int probe_median;
    //long probe_executions[nb_executions];
    int success = 1; 
    long executions[nb_executions];
    bool multiple_executions = false; 
    ofstream ExecutionFile("/users/cecilie/documents/code/identification-HE/result_frgc_" + to_string(nb_executions) + "_executions_" + argv[2] + "_" + argv[3] + ".txt");
    
    if (!multiple_executions) {nb_executions = 1;}

    for (int i = 0; i < nb_executions; i++) {
        if (multiple_executions) {cout << "\nIteration " << i << endl;}
        
        //
        //Prepare Identification
        //
        
        cleanUpArray[0] = 1;
        vector<double> cleanUpVector (cleanUpArray, cleanUpArray + sizeof(cleanUpArray) / sizeof(cleanUpArray[0]));
        Plaintext plaintextCleanUp = cc->MakeCKKSPackedPlaintext(cleanUpVector);
        Ciphertext<Element> encCleanUp = cc->Encrypt(keyPair.publicKey, plaintextCleanUp);
        
        //
        //Identification (working with internal subject IDs 0..1061, where 0..528 = FERET and 529..1061 = FRGC)
        //
        
        //bool k_packing = false; //test true for execution?
        
        start = std::chrono::steady_clock::now();
        
        //int searchForSubject = atoi(argv[2]);
        readTemplate(searchForSubject, argv[3], plainProbeIden, originalIDs);
        
        Ciphertext<Element> encProbeIden;
        encryptTemplate(encProbeIden, plainProbeIden, cc, keyPair);

        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        //auto rawProbeValue = elapsed.count();
        //probe_executions[i] = rawProbeValue;

        cout << "Probe Encryption took " << elapsed.count() << " milliseconds" << endl;
        
        start = std::chrono::steady_clock::now();

        int identifiedSubject = Identify(nrCandidates, originalIDs, encProbeIden, encCleanUp, encDB, candidateList, cc, keyPair);
        
        if(searchForSubject == identifiedSubject){
            cout << "Identification successful" << endl;
            success +=1; 
        }
        else cout << "ERROR: Subject " << searchForSubject << " incorrectly identified as " << identifiedSubject << endl;
        
        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        auto rawValue = elapsed.count();
        executions[i] = rawValue;
        sum += rawValue;
        
        
        if (!multiple_executions) {
            cout << "Identification for " << nrCandidates << " subjects took " << elapsed.count() << " milliseconds" << endl;
        }
        else {
            ExecutionFile << rawValue << endl;
            cout << "Identification took: " << elapsed.count() << " milliseconds" << endl;
        }
        //cout << rawValue << endl;
    }
    
    /*for (int i = 0; i< nb_executions; i++) {
        cout << executions[i] << endl; 
    }*/
    
    // Calculating median time
    
    if (multiple_executions) {
        sort(executions, executions + nb_executions);
        //sort(probe_executions, probe_executions + nb_executions);
        
        
        if (nb_executions % 2 != 0) {
            median = executions[nb_executions/2];
            //probe_median = probe_executions[nb_executions/2];
        }
        else {
            median = (executions[(nb_executions-1)/2] + executions[nb_executions/2])/2.0;
            //probe_median = (probe_executions[(nb_executions-1)/2] + probe_executions[nb_executions/2])/2.0;
        }
        
        cout << "\nTotal execution time for " << nb_executions << " executions was: " << sum/1000.0 << " seconds" << endl;
        cout << "Total successful identifications: " << success << endl;
        cout << "\nThe median time with " << nb_executions << " executions is: " << median << " milliseconds" << endl;
        ExecutionFile << "\nThe median time with " << nb_executions << " executions is: " << median << " milliseconds" << endl;
        
        
        //cout << "\nThe median time with " << nb_executions << " executions is: " << probe_median << " milliseconds" << endl;
    }
    
    
    return 0;
    
    
}