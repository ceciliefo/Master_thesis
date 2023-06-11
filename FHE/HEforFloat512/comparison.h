#ifndef MYOPENFHE_COMPARISON_H
#define MYOPENFHE_COMPARISON_H

#include "openfhe.h" 
#include "biometric.h"


using namespace lbcrypto;
using Element = DCRTPoly;

const double delta = 0.5;

void encryptED(Ciphertext<Element> &encED, Ciphertext<Element> &encProbe, Ciphertext<Element> &encCleanUp, Ciphertext<Element> &encRef, CryptoContext<Element> &cc);
double decryptED(Ciphertext<Element> &encED, CryptoContext<Element> &cc, KeyPair<Element> &keyPair);
int Identify(int &nrCandidates, int (&originalIDs)[nrSubj], Ciphertext<Element> &encProbe, Ciphertext<Element> &encCleanUp, Ciphertext<Element> (&encDB)[nrSubj], int (&candidateList)[nrSubj], CryptoContext<Element> &cc, KeyPair<Element> &keyPair);
int readCandidateList(int (&candidateList)[nrSubj], int &argc, char *argv[]);


#endif //MYOPENFHE_COMPARISON_H