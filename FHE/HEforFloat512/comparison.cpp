#include "comparison.h"
#include <iostream>

using namespace std;

void encryptED(Ciphertext<Element> &encED, Ciphertext<Element> &encProbe, Ciphertext<Element> &encRef, Ciphertext<Element> &encCleanUp, CryptoContext<Element> &cc){
    Ciphertext<Element> ciphertextSubTemplates = cc->EvalSub(encRef, encProbe);
    Ciphertext<Element> ciphertextSquaredTemplates = cc->EvalMult(ciphertextSubTemplates, ciphertextSubTemplates);
  
    Ciphertext<Element> ciphertextRot = cc->EvalAtIndex(ciphertextSquaredTemplates, 1);
    encED = cc->EvalAdd(ciphertextSquaredTemplates, ciphertextRot);
    for(int i = 0; i < templateSize - 1; i++){
        ciphertextRot = cc->EvalAtIndex(ciphertextRot, 1);
        encED = cc->EvalAdd(encED, ciphertextRot);
    }
    //encED = cc->EvalMult(encED, encCleanUp);
}

double decryptED(Ciphertext<Element> &encED, CryptoContext<Element> &cc, KeyPair<Element> &keyPair){
    Plaintext decED;
    cc->Decrypt(keyPair.secretKey, encED, &decED);
    return real(decED->GetCKKSPackedValue()[0]);
}

int Identify(int &nrCandidates, int (&originalIDs)[nrSubj], Ciphertext<Element> &encProbe, Ciphertext<Element> &encCleanUp, Ciphertext<Element> (&encDB)[nrSubj], int (&candidateList)[nrSubj], CryptoContext<Element> &cc, KeyPair<Element> &keyPair){
    Ciphertext<Element> encED;
    double ED = 0;
    int j = 0;
    double thresholdList[20][2] = {};
    for(int i = 0; i < nrCandidates; i++){
        if(j > 19){
            cout << "ERROR: thresholdList is too small" << endl;
            exit(1);
        }
        encryptED(encED, encProbe, encDB[candidateList[i]], encCleanUp, cc);
        ED = decryptED(encED, cc, keyPair);
        if(ED < delta){
            thresholdList[j][0] = i+1;
            thresholdList[j][1] = ED;
            j++;
        }
    }
    if(j == 0){
        cout << "Subject could not be identified" << endl;
        return 0;
    }
    int pos = 0;
    int smallest = thresholdList[pos][1];
    for (int i = 1; i < j; i++)
        if (smallest > thresholdList[i][1]) {
            smallest = thresholdList[i][1];
            pos = i;
        }
    return candidateList[(int) thresholdList[pos][0] - 1]; //working with internal subject IDs 0..1061
}

int readCandidateList(int (&candidateList)[nrSubj], int &argc, char *argv[]){
    int j = 0;
    string temp;
    if(argc<2){
        cout << "Please provide a candidate list txt file, a subject ID and a sample identifier" << endl;
        cout << "Try e.g. ../candidatelist.txt 124 940128_fb for FERET sample 208_940128_fb.txt or" << endl;
        cout << "         ../candidatelist.txt 903 d38 for FRGC sample 4763d38.txt" << endl;
        exit(1);
    } else{
        FILE* f;
        char stream[nrSubj] = {};
        f = fopen(argv[1], "r");
        if (f==NULL) perror ("Error opening file that was passed as argv[1]");
        while(fgets(stream, nrSubj, f) != NULL){
            if(stream[1] == 10){
                temp = stream[0];
            } else if(stream[2] == 10){
                temp = string(1, stream[0]) + string(1, stream[1]);
            } else if(stream[3] == 10){
                temp = string(1, stream[0]) + string(1, stream[1]) + string(1, stream[2]);
            } else{
                temp = string(1, stream[0]) + string(1, stream[1]) + string(1, stream[2]) + string(1, stream[3]);
            }
            candidateList[j] = stoi(temp);
            j++;
        }
    }
    int nrCandidates = 0;
    if(candidateList[nrSubj - 1] != 0){
        nrCandidates = nrSubj;
    } else{
        for(int i = 1; i < nrSubj; i++){
            if(candidateList[i] != 0){
                nrCandidates = i+1;
            }
        }
    }
    return nrCandidates;
}

