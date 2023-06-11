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


//#include "palisade.h" 
#include "openfhe.h" // OpenFHEreplace PALISADE
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

Ciphertext<Element> EncDB[nrSubj] = {};
string Samples[nrSamp + 4] = {};
string Subjects[nrSubj] = {};
int OriginalIDs[nrSubj] = {};
double PlainProbeIden[templateSize] = {};
double plainReferenceIden[templateSize] = {};
double PlainProbeIdenBatched[8*templateSize] = {};
double PlainReferenceIdenBatched[8*templateSize] = {};

void TestComparison(char *argv[]) {
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
    
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(PKE);
    
    KeyPair<Element> keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    cc->EvalAtIndexKeyGen(keyPair.secretKey, {1});
    
    setupEncDB(Subjects, Samples, EncDB, OriginalIDs, cc, keyPair);

/// TESTING WITHOUT BATCHING
    // read template before encrypt 
    // Reference template creation without batching
    int referenceSubject = atoi(argv[2]);
    string referenceSample = "940128_fa_a";
    readTemplate(referenceSubject, referenceSample, plainReferenceIden, OriginalIDs); 
    
    Ciphertext<Element> encReferenceIden;
    encryptTemplate(encReferenceIden, plainReferenceIden, cc, keyPair); //reference template
    
    // Start measuring
    auto start = std::chrono::steady_clock::now();
    
    // Probe template creation
    int searchForSubject = atoi(argv[2]); // this is the probe template
    readTemplate(searchForSubject, argv[3], PlainProbeIden, OriginalIDs);
    
    Ciphertext<Element> encProbeIden;
    encryptTemplate(encProbeIden, PlainProbeIden, cc, keyPair); 
    
    // Compute ED between probe and reference
    cleanUpArray[0] = 1;
    vector<double> cleanUpVector (cleanUpArray, cleanUpArray + sizeof(cleanUpArray) / sizeof(cleanUpArray[0]));
    Plaintext plaintextCleanUp = cc->MakeCKKSPackedPlaintext(cleanUpVector);
    Ciphertext<Element> encCleanUp = cc->Encrypt(keyPair.publicKey, plaintextCleanUp);
    
    Ciphertext<Element> encED;
    encryptED(encED, encProbeIden, encReferenceIden, encCleanUp, cc); 
    //encryptED between reference and probe template (where rotation happens)
    
    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    //decryptED(encProbeIden, cc, keyPair);
    decryptED(encED, cc, keyPair); // outside of the timings --> get the first entry of the ciphertext
    
    cout << "Computing the Euclidean distance between probe and reference ciphertext without batching took " << elapsed.count() << " milliseconds" << endl;
    //cout << "Ciphertext decryption without batching took " << elapsed.count() << " milliseconds" << endl;
    
    
/// TESTING WITH BATCHING
    // Reference template creation with batching
    // reference template --> with encryptBatchedTemplate
    // read in the template and then encrypt --> also concatenate the template
    // int referenceSubject = atoi(argv[2]);
    // string referenceSample = "940128_fa_a";
    int x = 8;
    vector<double>B;
    Ciphertext<Element> encReferenceIdenBatched;
    
    readTemplate(referenceSubject, referenceSample, plainReferenceIden, OriginalIDs); 
    while(x--){          
       for(int i=0;i<templateSize;i++){
           B.push_back(plainReferenceIden[i]);
       }
    }
    copy(B.begin(), B.end(), PlainReferenceIdenBatched);

    encryptBatchedTemplate(encReferenceIdenBatched, PlainReferenceIdenBatched, cc, keyPair);
    
    // Start measuring
    // With batching
    start = std::chrono::steady_clock::now();
    
    //int searchForSubject = atoi(argv[2]);
    int y = 8;
    vector<double>C;
    Ciphertext<Element> encProbeIdenBatched;
    
    readTemplate(searchForSubject, argv[3], PlainProbeIden, OriginalIDs);
    while(y--){          
       for(int i=0;i<templateSize;i++){
           C.push_back(PlainProbeIden[i]);
       }
    }
    copy(C.begin(), C.end(), PlainProbeIdenBatched);

    encryptBatchedTemplate(encProbeIdenBatched, PlainProbeIdenBatched, cc, keyPair);
    
    // Compute ED between probe and reference
    Ciphertext<Element> encED_batched;
    encryptED(encED_batched, encProbeIdenBatched, encReferenceIdenBatched, encCleanUp, cc); 
    // encryptED --> only time 101 and 155 --> the encED computations
    
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    //decryptED(encProbeIdenBatched, cc, keyPair);
    decryptED(encED_batched, cc, keyPair);
    
    cout << "Computing the Euclidean distance between probe and reference ciphertext with batching took " << elapsed.count() << " milliseconds" << endl;
    //cout << "Ciphertext decryption with batching took " << elapsed.count() << " milliseconds" << endl;
    
}
