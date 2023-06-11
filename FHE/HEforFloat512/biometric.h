#ifndef MYOPENFHE_BIOMETRIC_H
#define MYOPENFHE_BIOMETRIC_H

#include "openfhe.h" 

using namespace lbcrypto;
using Element = DCRTPoly;
using std::string; 

const int templateSize = 512;
const int nrSubjFERET = 529;
const int nrSampFERET = 1413;
const int nrSubjFRGC = 533;
const int nrSampFRGC = 3165;
const int nrSubj = nrSubjFRGC + nrSubjFERET;
const int nrSamp = nrSampFRGC + nrSampFERET;

double sciToDub(const string& str);
void readTemplate(int subjectID, string sample, double (&tempArr)[templateSize], int (&originalIDs)[nrSubj]);
void encryptTemplate(Ciphertext<Element> &encTemplate, double (&plainTemplate)[templateSize], CryptoContext<Element> &cc, KeyPair<Element> &keyPair);
void encryptBatchedTemplate(Ciphertext<Element> &encTemplate, double (&plainTemplate)[8*templateSize], CryptoContext<Element> &cc, KeyPair<Element> &keyPair);
void setupEncDB(string (&subjects)[nrSubj], string (&samples)[nrSamp + 4], Ciphertext<Element> (&encDB)[nrSubj], int (&originalIDs)[nrSubj], CryptoContext<Element> &cc, KeyPair<Element> &keyPair);
void setupEncDBwBatching(string (&subjects)[nrSubj], string (&samples)[nrSamp + 4], Ciphertext<Element> (&encDB)[nrSubj], int (&originalIDs)[nrSubj], CryptoContext<Element> &cc, KeyPair<Element> &keyPair);

#endif //MYOPENFHE_BIOMETRIC_H