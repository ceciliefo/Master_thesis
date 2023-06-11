#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>
#include <iostream> //new for Biometric_Test
#include <fstream> //test for getline with stable hash procedure
#include <string> // test for read files 
#include "dirent.h" // test for read filenames in directory
#include <filesystem> // test for read directory
#include <vector> //test read subject database
#include <algorithm> //test to sort files

#include "Sampling.h"
#include "params.h"
#include "FFT.h"
#include "Random.h"
#include "Algebra.h"

using namespace std;
using namespace NTL;
namespace fs = std::__fs::filesystem; 

const ZZX phi = Cyclo();


//==============================================================================
//Generates from parameters N and q :
// - a public key : polynomial h
// - a private key : polynomials f,g,F,G
//==============================================================================
void Keygen(ZZ_pX& PublicKey, ZZX* PrivateKey)
{
    ZZ SqNorm;
    ZZX f,g,F,G;

    SqNorm = conv<ZZ>(1.36*q0/2);

    GenerateBasis(f, g, F, G, SqNorm);
    PrivateKey[0] = f;
    PrivateKey[1] = g;
    PrivateKey[2] = F;
    PrivateKey[3] = G;

    for(unsigned int i=0; i<4; i++)
    {
            PrivateKey[i].SetLength(N0);
    }

    PublicKey = Quotient(f, g);
}

//==============================================================================
//Computes the private basis B from private key PrivateKey and parameter N
//==============================================================================
void CompletePrivateKey(mat_ZZ& B, const ZZX * const PrivateKey)
{
    ZZX f,g,F,G;
    f = PrivateKey[0];
    g = PrivateKey[1];
    F = PrivateKey[2];
    G = PrivateKey[3];

    f = -f;
    F = -F;

    B = BasisFromPolynomials(g, f, G, F);
}





void GPV(RR_t * v, const RR_t * const c, const RR_t s, const MSK_Data * const MSKD) //lattice-based trapdoor
{

    int i;
    unsigned j;
    RR_t ci[2*N0], zi, cip, sip, aux;

    for(j=0; j<2*N0;j++)
    {
        ci[j] = c[j];
    }

//    for(j=0; j<2*N0; j++)
//    {
//
//    }    

    for(i=2*N0-1; i>=0; i--)
    {
        aux = (MSKD->GS_Norms)[i];
        cip = DotProduct(ci, MSKD->Bstar[i])/(aux*aux);
        sip = s/aux;
        zi = Sample4(cip, sip*PiPrime);

        for(j=0; j<2*N0; j++)
        {
            ci[j] -= zi*(MSKD->B)[i][j];
        }
    }

    for(j=0; j<2*N0; j++)
    {
        v[j] = c[j] - ci[j];
    }

}



//==============================================================================
//==============================================================================
//                            MAIN PROGRAMS
//==============================================================================
//==============================================================================


void CompleteMSK(MSK_Data * MSKD, ZZX * MSK)
{
    unsigned int i, j;
    mat_ZZ B0;

    for(i=0; i<4; i++)
    {
        MSKD->PrK[i] = MSK[i];
        ZZXToFFT(MSKD->PrK_fft[i], MSK[i]);
    }

    CompletePrivateKey(B0, MSK);

    for(i=0; i<2*N0; i++)
    {
        for(j=0; j<2*N0; j++)
        {
            MSKD->B[i][j] = ( (RR_t) conv<double>(B0[i][j]) );
        }
    }

    for(i=0; i<1; i++)
    {
        FastMGS(MSKD->Bstar, MSKD->B);
    }

    for(i=0; i<2*N0; i++)
    {
        MSKD->GS_Norms[i] = sqrt( DotProduct(MSKD->Bstar[i], MSKD->Bstar[i]) );
    }

    MSKD->sigma = 2*MSKD->GS_Norms[0];

}



void CompleteMPK(MPK_Data * MPKD, ZZ_pX MPK)
{
    MPKD->h = MPK;
    ZZXToFFT(MPKD->h_FFT, conv<ZZX>(MPK));
}



void PEKS_Trapdoor(ZZX SK_tr[2], vec_ZZ kw, const MSK_Data * const MSKD) //takes secret key and keyword vector --> produces trapdoor for keyword vector?
{
    unsigned int i;
    RR_t c[2*N0], sk[2*N0], sigma;
    ZZX f,g,aux;

    f = MSKD -> PrK[0];
    g = MSKD -> PrK[1];
    sigma = MSKD->sigma;
    SK_tr[0].SetLength(N0);
    SK_tr[1].SetLength(N0);

    for(i=0;i<N0;i++)
    {
        c[i] = ((RR_t) conv<double>(kw[i])) ;
        c[i+N0] = 0;
    }

    GPV(sk, c, sigma, MSKD); //lattice-based trapdoor
    //uses c, sigma, MSKD to write in sk

    for(i=0; i<N0; i++)
    {
        sk[i] = c[i] - sk[i];
        sk[i+N0] = - sk[i+N0];
    }

    for(i=0; i<N0; i++)
    {
        SK_tr[0][i] = sk[i]; //first half of sk
        SK_tr[1][i] = sk[i+N0]; //second half of sk
    }
    
}


unsigned long PEKS_Verify_Trapdoor(const ZZX SK_tr[2], const vec_ZZ kw, const MSK_Data * const MSKD)
{
    unsigned int i;
    ZZX f,g,t,aux;

    f = MSKD -> PrK[0];
    g = MSKD -> PrK[1];
    
    t = conv<ZZX>(kw);
    aux = ((SK_tr[0] - t)*f + g*SK_tr[1])%phi;

    for(i=0; i<N0; i++)
    {
        aux[i] %= q1;
    }

    if( IsZero(aux) != 0)
    {
        cout << "The signature (s1,s2) doesn't verify the required equality [ (s1 - t)*f + g*s2 = 0 ] !\nActually, (s1 - t)*f + g*s2 = " << aux << endl << endl;
    }
    return IsZero(aux);
}


void PEKS_Enc(long C[3][N0],  const long id0[N0],  const MPK_Data * const MPKD) 
{

    unsigned long i;
    long r[N0], e1[N0], e2[N0];
    CC_t r_FFT[N0], t_FFT[N0], aux1_FFT[N0], aux2_FFT[N0];


    for(i=0; i<N0; i++)
    {
        e1[i] = (rand()%3) - 1;
        e2[i] = (rand()%3) - 1;
        r[i] = (rand()%3) - 1;
        C[2][i] = (rand()%2);
    }


    MyIntFFT(r_FFT, r);
    MyIntFFT(t_FFT, id0);

    for(i=0; i<N0; i++)
    {
        aux1_FFT[i] = r_FFT[i]*((MPKD->h_FFT)[i]);
        aux2_FFT[i] = r_FFT[i]*t_FFT[i];
    } 

    MyIntReverseFFT(C[0], aux1_FFT);
    MyIntReverseFFT(C[1], aux2_FFT);

    for(i=0; i<N0; i++)
    { 
        C[0][i] = (C[0][i] + e1[i]               + q0/2)%q0 - (q0/2);
        C[1][i] = (C[1][i] + e2[i] + (q0/2)*C[2][i] + q0/2)%q0 - (q0/2);
    } 

}


bool PEKS_Test( const long C[3][N0], const CC_t * const SKtd_FFT)
{
    unsigned int i;
    CC_t c0_FFT[N0], aux_FFT[N0];
    bool fout = false;
    long k[N0];
    MyIntFFT(c0_FFT, C[0]);

    for(i=0; i<N0; i++)
    {
        aux_FFT[i] = c0_FFT[i]*SKtd_FFT[i];
    }

    MyIntReverseFFT(k, aux_FFT);
      
    int z = 0; //Zähler
  
    for(i=0; i<N0; i++)
    {
        k[i] = C[1][i] - k[i];
        k[i] = ((unsigned long)(k[i] ))%q0;
        k[i] = (k[i] + (q0>>2) )/(q0>>1);
        k[i] %= 2;
        if (C[2][i] == k[i]) //==
        {
            z++;
            //cout << "1";
        }
        /*
        else
        {
            cout << "0";
        }
        */
     }
     
     if (z == N0)
     {
         fout = true;
     }
   /* if (fout)
        cout << endl<< endl<<endl <<"  TEST algorithm has been succesful"<< endl; */
    //cout << endl << endl;
    return fout;

}



//==============================================================================
//==============================================================================
//                             BENCHES AND TESTS
//                   FOR EXTRACTION AND ENCRYPTION/DECRYPTION
//==============================================================================
//==============================================================================


void Trapdoor_Bench(const unsigned int nb_extr, MSK_Data * MSKD)
{
    clock_t t1, t2;
    float diff;
    unsigned int i;
    vec_ZZ kw;
    ZZX SK_tr[2];

    t1 = clock();

    cout << "0%" << flush;
    for(i=0; i<nb_extr; i++)
    {
        kw = RandomVector();

        PEKS_Trapdoor(SK_tr, kw, MSKD);
        if((i+1)%(nb_extr/10)==0)
        {
            cout << "..." << (i+1)/(nb_extr/10) << "0%" << flush;
        }
    }

    t2 = clock();
    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << "\n\nIt took " << diff << " seconds to create  " << nb_extr << " trapdoors." << endl;
    cout << "That's " << (diff/nb_extr)*1000 << " milliseconds per trapdoor." << endl << endl;
}


void Encrypt_Bench(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD)
{
    clock_t te1, te2, td1, td2;
    float diffe = 0, diffd = 0;
    unsigned int i;//,j;
    vec_ZZ kw;
    ZZX SK_tr[2], w;
    CC_t SKid_FFT[N0];
   // long int message[N0], decrypted[N0];
    long int keyword[N0], Ciphertext[3][N0];
	//bool flag = true; 

    kw = RandomVector();
    PEKS_Trapdoor(SK_tr, kw, MSKD);
    PEKS_Verify_Trapdoor(SK_tr, kw, MSKD);
    ZZXToFFT(SKid_FFT, SK_tr[1]);

    for(i=0; i<N0; i++)
    {
        keyword[i] = conv<long int>(kw[i]);
    }

    cout << "0%" << flush ;
    for(i=0; i<nb_cryp; i++)
    {

   /*     for(j=0; j<N0; j++)
        {
            message[j] = (rand()%2);
        }*/
/*		for(j=0; j<N0; j++)
        {
             cout << " "<< j<<"th "<< message[j] << "\t"; ;
        }*/
	te1 = clock();	
		PEKS_Enc(Ciphertext, keyword, MPKD);
	te2 = clock();

	

	td1 = clock();
		if (!PEKS_Test(Ciphertext, SKid_FFT)){
			 cout << "TEST FAILED --- Exiting..."<<endl;
			break;
			}
			
	td2 = clock();
	
        if((i+1)%(nb_cryp/10)==0)
        {
            cout << "..." << (i+1)/(nb_cryp/10) << "0%" << flush;
        }

	diffe += ((float)te2 - (float)te1)/1000000.0l;
	diffd += ((float)td2 - (float)td1)/1000000.0l;

    }

    cout << "\n\nIt took " << diffe << " seconds to do " << nb_cryp << " encryptions." << endl;
    cout << "That's " << (diffe/nb_cryp)*1000 << " milliseconds per PEKS generation." << endl;
    cout << "That's " << (diffe/nb_cryp)*1000*1024/N0 << " milliseconds per PEKS per Kilobit." << endl << endl;

    cout << "\n\nIt took " << diffd << " seconds to do " << nb_cryp << " Tests." << endl;
    cout << "That's " << (diffd/nb_cryp)*1000 << " milliseconds per Tests." << endl;
    cout << "That's " << (diffd/nb_cryp)*1000*1024/N0 << " milliseconds per Tests per Kilobit." << endl << endl;

}


void Trapdoor_Test(const unsigned int nb_extr, MSK_Data * MSKD)
{
    unsigned int i, rep;
    vec_ZZ kw;
    ZZX SK_kw[2];

    rep = 0;

    cout << "0%" << flush;
    for(i=0; i<nb_extr; i++)
    {
        kw = RandomVector();

        PEKS_Trapdoor(SK_kw, kw, MSKD);
        rep += PEKS_Verify_Trapdoor(SK_kw, kw, MSKD);
        if((i+1)%(nb_extr/10)==0)
        {
            cout << "..." << (i+1)/(nb_extr/10) << "0%" << flush;
        }
    }

    cout << endl;
    if(rep == 0)
    {    cout << endl << nb_extr << " Trapdoor successfully performed!" << endl << endl;    }
    else
    {    cout << endl << rep << " out of " << nb_extr << " extractions failed miserabily!" << endl << endl;    }
}


void Encrypt_Test(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD)
{
    unsigned int i, rep;//j, rep;
    vec_ZZ kw;
    ZZX SK_td[2], m;
    CC_t SKtd_FFT[N0];
    long int kw0[N0], Ciphertext[2][N0];
   // long int message[N0], decrypted[N0];


    kw = RandomVector();
    PEKS_Trapdoor(SK_td, kw, MSKD);
    PEKS_Verify_Trapdoor(SK_td, kw, MSKD);
    ZZXToFFT(SKtd_FFT, SK_td[1]);

    rep = 0;

    for(i=0; i<N0; i++)
    {
        kw0[i] = conv<long int>(kw[i]);
	

    }

    cout << "0%" << flush;
    for(i=0; i<nb_cryp; i++)
    {

       /* for(j=0; j<N0; j++)
        {
            message[j] = (rand()%2);
        }*/

		PEKS_Enc(Ciphertext, kw0, MPKD);
        PEKS_Test(Ciphertext, SKtd_FFT);
        
/*        for(j=0; j<N0; j++)
        {
            if(message[j] != decrypted[j])
            {
                cout << "ERROR : Dec(Enc(m)) != m " << endl;
                rep++;
                break;
            }
        }*/

        if((i+1)%(nb_cryp/10)==0)
        {
            cout << "..." << (i+1)/(nb_cryp/10) << "0%" << flush;
        }
    }

    cout << endl;
    if(rep == 0)
    {    cout << endl << nb_cryp << " PEKS+TEST successfully performed!" << endl << endl;    }
    else
    {    cout << endl << rep << " out of " << nb_cryp << " PEKS+TEST failed miserabily!" << endl << endl;    }
}

//Setup keyword database
void Keyword_Database(CC_t (&keys)[24][N0], long int (&kwdb)[nb_kwcl][nb_sub][3][N0], MPK_Data * MPKD, MSK_Data * MSKD){
    unsigned int i, j;
    vec_ZZ kws[24];
    ZZX tempkeys[24][2];
    long int int_kws[24][N0], keyword_list0[nb_sub], keyword_list1[nb_sub], keyword_list2[nb_sub], keyword_list3[nb_sub];
    FILE * tempfile;
    char line [20];
    
    for(i=0; i<24; i++){
        kws[i] = RandomVector();
        PEKS_Trapdoor(tempkeys[i], kws[i], MSKD);
        PEKS_Verify_Trapdoor(tempkeys[i], kws[i], MSKD);
        ZZXToFFT(keys[i], tempkeys[i][1]);
        for(j=0; j<N0; j++){
            int_kws[i][j] = conv<long int>(kws[i][j]);
        }
    }
    
    //class 0 = gender
    tempfile = fopen("/Users/cecilie/Documents/Masterproject/NTRU-PEKS/frgc_gender_10000", "r");
    if (tempfile==NULL) perror ("Error opening file");
    else {
        for(j=0; j<nb_sub; j++){
            fgets(line, 20, tempfile);
            if(line[0] == 109){
                keyword_list0[j] = 0; //ASCII m = 109
            } else if(line[0] == 102){ //ASCII f = 102
                keyword_list0[j] = 1;
            } else{
                cout << "Error: Unknown keyword in frgc_gender!" << endl;
            }
        }
    fclose(tempfile);
    }
    
    //Setup kwdb[0]
    for(j=0; j<nb_sub; j++){
        if(keyword_list0[j] == 0) PEKS_Enc(kwdb[0][j], int_kws[0], MPKD);
        else PEKS_Enc(kwdb[0][j], int_kws[1], MPKD);
    }
    
    
    //class 1 = ethnicity
    tempfile = fopen("/Users/cecilie/Documents/Masterproject/NTRU-PEKS/frgc_ethnicity_10000", "r");
    if (tempfile==NULL) perror ("Error opening file");
    else {
        for(j=0; j<nb_sub; j++){
            fgets(line, 20, tempfile);
            if(line[0] == 99){
                keyword_list1[j] = 0; //caucasian
            } else if(line[0] == 104){
                keyword_list1[j] = 2; //hispanic
            } else if(line[0] == 105){
                keyword_list1[j] = 4; //indian
            } else if(line[0] == 98){
                keyword_list1[j] = 5; //black
            } else if(line[0] == 97){
                if(line[1] == 115){
                    keyword_list1[j] = 1; //asian
                } else if(line[1] == 114){
                    keyword_list1[j] = 3; //arabic
                }
            } else{
                cout << "Error: Unknown keyword in frgc_ethnicity!" << endl;
            }
        }
    fclose(tempfile);
    }
    
    //Setup kwdb[1]
    for(j=0; j<nb_sub; j++){
        if(keyword_list1[j] == 0) PEKS_Enc(kwdb[1][j], int_kws[2], MPKD);
        else if(keyword_list1[j] == 1) PEKS_Enc(kwdb[1][j], int_kws[3], MPKD);
        else if(keyword_list1[j] == 2) PEKS_Enc(kwdb[1][j], int_kws[4], MPKD);
        else if(keyword_list1[j] == 3) PEKS_Enc(kwdb[1][j], int_kws[5], MPKD);
        else if(keyword_list1[j] == 4) PEKS_Enc(kwdb[1][j], int_kws[6], MPKD);
        else if(keyword_list1[j] == 5) PEKS_Enc(kwdb[1][j], int_kws[7], MPKD);
    }
    
    //class 2 = skin type
    tempfile = fopen("/Users/cecilie/Documents/Code/NTRU-PEKS/frgc_skintype_10000", "r");
    if (tempfile==NULL) perror ("Error opening file");
    else {
        for(j=0; j<nb_sub; j++){
            fgets(line, 20, tempfile);
            if(line[0] == 49){
                if(line[1] == 10){
                    keyword_list2[j] = 0; //skin type 1 // ASCII 10 = LF
                }
                else if(line[1] == 48){
                    keyword_list2[j] = 9; //skin type 4-5
                }
                else if(line[1] == 49){
                    keyword_list2[j] = 10; //skin type 5-6
                }
                else cout << "Error: Unknown keyword in frgc_skintype!" << endl;
            } else if(line[0] == 50){
                keyword_list2[j] = 1; //skin type 2
            } else if(line[0] == 51){
                keyword_list2[j] = 2; //skin type 3
            } else if(line[0] == 52){
                keyword_list2[j] = 3; //skin type 4
            } else if(line[0] == 53){
                keyword_list2[j] = 4; //skin type 5
            } else if(line[0] == 54){
                keyword_list2[j] = 5; //skin type 6
            } else if(line[0] == 55){
                keyword_list2[j] = 6; //skin type 1-2
            } else if(line[0] == 56){
                keyword_list2[j] = 7; //skin type 2-3
            } else if(line[0] == 57){
                keyword_list2[j] = 8; //skin type 3-4
            } else{
                cout << "Error: Unknown keyword in frgc_skintype!" << endl;
            }
        }
    fclose(tempfile);
    }
    
    //Setup kwdb[2]
    for(j=0; j<nb_sub; j++){
        if(keyword_list2[j] == 0) PEKS_Enc(kwdb[2][j], int_kws[8], MPKD); //skin type 1
        else if(keyword_list2[j] == 1) PEKS_Enc(kwdb[2][j], int_kws[9], MPKD); //skin type 2
        else if(keyword_list2[j] == 2) PEKS_Enc(kwdb[2][j], int_kws[10], MPKD); //skin type 3
        else if(keyword_list2[j] == 3) PEKS_Enc(kwdb[2][j], int_kws[11], MPKD); //skin type 4
        else if(keyword_list2[j] == 4) PEKS_Enc(kwdb[2][j], int_kws[12], MPKD); //skin type 5
        else if(keyword_list2[j] == 5) PEKS_Enc(kwdb[2][j], int_kws[13], MPKD); //skin type 6
        else if(keyword_list2[j] == 6) PEKS_Enc(kwdb[2][j], int_kws[14], MPKD); //skin type 7
        else if(keyword_list2[j] == 7) PEKS_Enc(kwdb[2][j], int_kws[15], MPKD); //skin type 8
        else if(keyword_list2[j] == 8) PEKS_Enc(kwdb[2][j], int_kws[16], MPKD); //skin type 9
        else if(keyword_list2[j] == 9) PEKS_Enc(kwdb[2][j], int_kws[17], MPKD); //skin type 10
        else if(keyword_list2[j] == 10) PEKS_Enc(kwdb[2][j], int_kws[18], MPKD); //skin type 11
    }
    
    //class 3 = age group
    tempfile = fopen("/Users/cecilie/Documents/Code/NTRU-PEKS/frgc_agegroup_10000", "r");
    if (tempfile==NULL) perror ("Error opening file");
    else {
        for(j=0; j<nb_sub; j++){
            fgets(line, 20, tempfile);
            if(line[0] == 121){
                if(line[5] == 10){ //ASCII 10 = LF
                    keyword_list3[j] = 0; //age group young
                }
                else if(line[5] == 45){ //ASCII 45 = -
                    keyword_list3[j] = 1; //age group young-middle
                } else cout << "Error: Unknown keyword in frgc_agegroup!" << endl;
            } else if(line[0] == 109){
                if(line[6] == 10){
                    keyword_list3[j] = 2; //age group middle
                }
                else if(line[6] == 45){
                    keyword_list3[j] = 3; //age group middle-old
                } else cout << "Error: Unknown keyword in frgc_agegroup!" << endl;
            } else if(line[0] == 111){
                keyword_list3[j] = 4; //age group old
            } else{
                cout << "Error: Unknown keyword in frgc_agegroup!" << endl;
            }
        }
    fclose(tempfile);
    }
    
    //Setup kwdb[3]
    for(j=0; j<nb_sub; j++){
        if(keyword_list3[j] == 0) PEKS_Enc(kwdb[3][j], int_kws[19], MPKD);
        else if(keyword_list3[j] == 1) PEKS_Enc(kwdb[3][j], int_kws[20], MPKD);
        else if(keyword_list3[j] == 2) PEKS_Enc(kwdb[3][j], int_kws[21], MPKD);
        else if(keyword_list3[j] == 3) PEKS_Enc(kwdb[3][j], int_kws[22], MPKD);
        else if(keyword_list3[j] == 4) PEKS_Enc(kwdb[3][j], int_kws[23], MPKD);
    }
    
    cout << "Database setup complete" << endl << endl;
}

void Single_Keyword_Search(int kw_class, int keyword, long int (&kwdb)[nb_kwcl][nb_sub][3][N0], CC_t (&keys)[nb_kw][N0], long int (&result)[nb_sub], MPK_Data * MPKD, MSK_Data * MSKD){
    unsigned int i, j;
   
    for(j=0; j<nb_sub; j++){
        result[j] = 0;
    }
    
    i = 0;
    for(j=0; j<nb_sub; j++){
        if(PEKS_Test(kwdb[kw_class][j], keys[keyword]) == 1){ //kwdb[0] = gender, kwdb[1] = ethnicity
            result[i] = j;
            i++;
        }
    }
    
    
}

void Cascading_Keyword_Search(int (&search_for)[2*nb_kwcl_search], long int (&kwdb)[nb_kwcl][nb_sub][3][N0], CC_t (&keys)[nb_kw][N0], long int (&result)[nb_sub], MPK_Data * MPKD, MSK_Data * MSKD){
    unsigned int i, j, k, temp_nb_sub;
    long int tempresults[nb_kwcl_search+1][nb_sub];
    
    for(j=0; j<nb_sub; j++){
        tempresults[0][j] = j;
    }
    
    for(i=1; i<nb_kwcl_search+1; i++){
        for(j=0; j<nb_sub; j++){
            tempresults[i][j] = 0;
        }
    }
    
    temp_nb_sub = nb_sub;
    for(k=0; k<nb_kwcl_search; k++){
        i = 0;
        for(j=0; j<temp_nb_sub; j++){
            if(PEKS_Test(kwdb[search_for[2*k]][tempresults[k][j]], keys[search_for[(2*k)+1]]) == 1){
                tempresults[k+1][i] = tempresults[k][j];
                i++;
            }
        }
        temp_nb_sub = i;
    }
    
    for(j=0; j<nb_sub; j++){
        result[j] = tempresults[nb_kwcl_search][j];
    }
}

void Biometric_Test_Old(long int (&kwdb)[nb_kwcl][nb_sub][3][N0], MPK_Data * MPKD, MSK_Data * MSKD)
{    
    int i, rep = 0;
    CC_t keys[nb_kw][N0];
    long int result[nb_sub];
    
    enum kw_classes{gender=0, ethnicity=1, skintype=2, agegroup=3};
    enum keywords{male=0, female=1, caucasian=2, asian=3, hispanic=4, arabic=5, indian=6, black=7, skintype1=8, skintype2=9, skintype3=10, skintype4=11, skintype5=12, skintype6=13, skintype1_2=14, skintype2_3=15, skintype3_4=16, skintype4_5=17, skintype5_6=18, young=19, young_middle=20, middle=21, middle_old=22, old=23};
    
    int search_for[2*nb_kwcl_search] = {gender, female, ethnicity, asian};
    
    //Setup keyword database
    Keyword_Database(keys, kwdb, MPKD, MSKD);

    //Search for one keyword
    Single_Keyword_Search(gender, male, kwdb, keys, result, MPKD, MSKD);
    
    cout << "keyword search for <male>" << endl;
    for(i=0; i<nb_sub; i++){
        cout << result[i] << ", ";
    }
    
    //Search for multiple keywords
    Cascading_Keyword_Search(search_for, kwdb, keys, result, MPKD, MSKD);
    
    //0, 0, 1, 3 -> male, asian
    //0, 1, 1, 7 -> female, black
    //1, 2, 0, 0 -> caucasian, male
    //2, 11, 0, 1 -> skin type 4, female
    
    cout << endl << endl << "keyword search for <female>, <asian>" << endl;
    for(i=0; i<nb_sub; i++){
        cout << result[i] << ", ";
    }
    
    if(rep == 0)
    {    cout << endl << "Biometric keyword search successfully performed!" << endl << endl;    }
    else
    {    cout << endl << rep << "Biometric keyword search failed miserabily!" << endl << endl;    }
}

void Biometric_Bench_Old(long int (&kwdb)[nb_kwcl][nb_sub][3][N0], MPK_Data * MPKD, MSK_Data * MSKD){
    
    CC_t keys[nb_kw][N0];
    long int result[nb_sub];
    clock_t t1, t2, t5, t6;//, t8, t3, t4, t7;
    float diff1, diff3;// diff4; //diff2;
    
    enum kw_classes{gender=0, ethnicity=1, skintype=2, agegroup=3};
    enum keywords{male=0, female=1, caucasian=2, asian=3, hispanic=4, arabic=5, indian=6, black=7, skintype1=8, skintype2=9, skintype3=10, skintype4=11, skintype5=12, skintype6=13, skintype1_2=14, skintype2_3=15, skintype3_4=16, skintype4_5=17, skintype5_6=18, young=19, young_middle=20, middle=21, middle_old=22, old=23};
    
    //two keywords
    //int search_for_bc[2*nb_kwcl_search] = {agegroup, middle_old, ethnicity, asian};
    //int search_for_wc[2*nb_kwcl_search] = {ethnicity, caucasian, gender, male};

    //three keywords
    //int search_for_bc[2*nb_kwcl_search] = {agegroup, middle_old, ethnicity, asian, skintype, skintype4};
    //int search_for_wc[2*nb_kwcl_search] = {ethnicity, caucasian, gender, male, skintype, skintype2};

    //four keywords
    //int search_for_bc[2*nb_kwcl_search] = {agegroup, middle_old, ethnicity, asian, skintype, skintype4, gender, male};
    //int search_for_wc[2*nb_kwcl_search] = {ethnicity, caucasian, gender, male, skintype, skintype2, agegroup, middle};
    
    t5 = clock();
    //Setup keyword database
    Keyword_Database(keys, kwdb, MPKD, MSKD);
    t6 = clock();
    diff3 = ((float)t6 - (float)t5)/1000000.0l;
    
    cout << "\nKeyword_Database took " << diff3 << " seconds for " << nb_sub << " subjects" << endl;
    
    t1 = clock();
    //Search for one keyword
    Single_Keyword_Search(gender, male, kwdb, keys, result, MPKD, MSKD);
    t2 = clock();
    diff1 = ((float)t2 - (float)t1)/1000000.0l;
    
    cout << "\n\nSingle_Keyword_Search took " << diff1 << " seconds for " << nb_sub << " subjects" << endl;
    
    /*
    t3 = clock();
    //Search for multiple keywords best case
    Cascading_Keyword_Search(search_for_bc, kwdb, keys, result, MPKD, MSKD);
    t4 = clock();
    diff2 = ((float)t4 - (float)t3)/1000000.0l;
    
    cout << "\n\nBest-Case Cascading_Keyword_Search took " << diff2 << " seconds for " << nb_sub << " subjects" << endl;

    t7 = clock();
    //Search for multiple keywords worst case
    Cascading_Keyword_Search(search_for_wc, kwdb, keys, result, MPKD, MSKD);
    t8 = clock();
    diff4 = ((float)t8 - (float)t7)/1000000.0l;
    
    cout << "\n\nWorst-Case Cascading_Keyword_Search took " << diff4 << " seconds for " << nb_sub << " subjects" << endl;
     */
}

void Gen_Trapdoor_List(CC_t (&keys)[nb_classes][N0], long int (&int_kws)[nb_classes][N0], MSK_Data * MSKD){ //reverse PEKS(figure) , 64 random kw, generate trapdoor, 
    unsigned int i, j;
    vec_ZZ kws[nb_classes];
    ZZX tempkeys[nb_classes][2];
    
    for(i=0; i<nb_classes; i++){
        kws[i] = RandomVector();
        PEKS_Trapdoor(tempkeys[i], kws[i], MSKD);
        PEKS_Verify_Trapdoor(tempkeys[i], kws[i], MSKD);
        ZZXToFFT(keys[i], tempkeys[i][1]);
        for(j=0; j<N0; j++){
            int_kws[i][j] = conv<long int>(kws[i][j]);
        }
    }
}

void Gen_Lookup_Table(long int (&LookupTable)[nb_ne_classes]){ //binning
    for(int i = 0; i < 8; i++){
        LookupTable[i] = i;
    }
    LookupTable[8] = 2;
    LookupTable[9] = 3;
    LookupTable[10] = 4;
    LookupTable[11] = 5;
    LookupTable[12] = 6;
    LookupTable[13] = 7;
    LookupTable[14] = 7;
    LookupTable[15] = 6;
    LookupTable[16] = 8;
    LookupTable[17] = 8;
    LookupTable[18] = 8;
    LookupTable[19] = 5;
    LookupTable[20] = 4;
    LookupTable[21] = 8;
    LookupTable[22] = 8;
    LookupTable[23] = 8;
    for(int i = 24; i < 28; i++){
        LookupTable[i] = 9;
    }
    LookupTable[28] = 8;
    LookupTable[29] = 9;
    LookupTable[30] = 9;
    LookupTable[31] = 3;
    LookupTable[32] = 1;
    for(int i = 33; i < 37; i++){
        LookupTable[i] = 9;
    }
    LookupTable[37] = 10;
    LookupTable[38] = 10;
    LookupTable[39] = 2;
    LookupTable[40] = 7;
    for(int i = 41; i < 56; i++){
        LookupTable[i] = 10;
    }
    LookupTable[56] = 11;
    LookupTable[57] = 11;
    LookupTable[58] = 9;
    for(int i = 59; i < 73; i++){
        LookupTable[i] = 11;
    }
    LookupTable[73] = 10;
    for(int i = 74; i < 85; i++){
        LookupTable[i] = 11;
    }
    for(int i = 85; i < 155; i++){
        LookupTable[i] = 12;
    }
}

void Biometric_Test(CC_t (&keys)[nb_classes][N0], MPK_Data * MPKD, MSK_Data * MSKD){
    
    long int Probe[3][N0], int_kws[nb_classes][N0], LookupTable[nb_ne_classes];
    int found = 0; //int i
    
    Gen_Trapdoor_List(keys, int_kws, MSKD); 
    Gen_Lookup_Table(LookupTable);
    
    PEKS_Enc(Probe, int_kws[120], MPKD);
    
    for(int i = 0; i < nb_ne_classes; i++){
        if(PEKS_Test(Probe, keys[i]) == 1){
            found = 1;
            cout << "PEKS_Test == 1 for i = " << i << " and probe is in bin " << LookupTable[i] << endl;
        }
    }
    if(found == 0){
        cout << "No subjects with these keywords is enrolled" << endl;
    }
    
    cout << endl << "Biometric_Test successfully performed!" << endl << endl; 

}

void Biometric_Bench(CC_t (&keys)[nb_classes][N0], MPK_Data * MPKD, MSK_Data * MSKD){
    
    long int Probe[3][N0], int_kws[nb_classes][N0], LookupTable[nb_ne_classes];
    int found = 0; //int i
    clock_t t1, t2, t3, t4;
    float diff1, diff2;
    
    t1 = clock();
    //Setup keyword database
    Gen_Trapdoor_List(keys, int_kws, MSKD);
    //Gen_Lookup_Table(LookupTable);
    t2 = clock();
    diff1 = ((float)t2 - (float)t1)/1000000.0l;
    
    cout << "Trapdoor List and Lookup Table Setup took " << diff1 << " seconds for " << nb_ne_classes << " non-empty keyword classes" << endl;
    
    PEKS_Enc(Probe, int_kws[77], MPKD);
    
    t3 = clock();
    //Search for one keyword
    for(int i = 0; i < nb_ne_classes; i++){
        if(PEKS_Test(Probe, keys[i]) == 1){
            found = 1;
            cout << "\nPEKS_Test == 1 for i = " << i << " and probe is in bin " << LookupTable[i] << endl;
        }
    }
    t4 = clock();
    diff2 = ((float)t4 - (float)t3)/1000000.0l;

    if(found == 0){
        cout << "No subjects with these keywords is enrolled" << endl;
    }
    
    cout << "\nKeyword Search took " << diff2 << " seconds for " << nb_ne_classes << " non-empty keyword classes" << endl;
}




//==============================================================================
//My Implementation: Use Stable Hash Codes (Cluster Centers) As Keywords
//==============================================================================

void Generate_Cluster_Trapdoors(CC_t (&cluster_trapdoors)[nb_kw_clusters][N0], long int (&int_kws)[nb_kw_clusters][N0], MSK_Data * MSKD){
    // Generate 64 trapdoors from 64 random keyword vectors 
    // cluster_trapdoors: stores the trapdoor for each keyword cluster
    // int_kws: stores the plaintext keyword 
    
    unsigned int i, j; 
    vec_ZZ kws[nb_kw_clusters];
    ZZX temp_trapdoors[nb_kw_clusters][2]; 
    
    for (i=0; i<nb_kw_clusters; i++) {
        kws[i] = RandomVector(); // Initalize the kws array as a randon vector
        PEKS_Trapdoor(temp_trapdoors[i], kws[i], MSKD); // Produce a trapdoor and store the trapdoor for keyword cluster kws[i] in temp_trapdoors[i]
        PEKS_Verify_Trapdoor(temp_trapdoors[i], kws[i], MSKD); // Verify computed correct trapdoor 
        ZZXToFFT(cluster_trapdoors[i], temp_trapdoors[i][1]); // Transform the trapdoors from temp_trapdoors to cluster_trapdoors
        
        for (j=0; j<N0; j++) {
            int_kws[i][j] = conv<long int>(kws[i][j]); // Store the plain keyword (random vector) kws[i] in int_kws[i]
        }
    }
    
    //cout << "Generation of trapdoors for each keyword cluster are completed" << endl << endl; 
}

void Encrypted_Cluster_Keyword_Database(long int (&enc_kwdb)[nb_subjects][3][N0], long int (&int_kws)[nb_kw_clusters][N0], MPK_Data * MPKD, MSK_Data * MSKD){
    // Setup the encrypted keyword database
    // Generate searchable ciphertexts for each subject based on their keyword cluster
    // enc_kwdb: stores the searchable ciphertext for each subject 
    // int_kws: stores the plaintext keyword 
    
    unsigned int i; 
    string line; 
    fstream tempfile; 
    
    tempfile.open("/users/cecilie/documents/code/identification/features/frgc_arcface512_float/workload_k_means/frgc_clusters.txt", ios::in); //read the keyword cluster file of frgc--> corresponding keyword cluster for each subject
    //tempfile.open("/users/cecilie/documents/code/identification/features/feret_arcface512_float/workload_k_means/feret_clusters.txt", ios::in); //read the keyword cluster file of feret--> corresponding keyword cluster for each subject

    if (tempfile.fail()) perror("Error with opening the file"); 
    else if (tempfile.is_open()) {
        for (i=0; i<nb_subjects; i++) {
            getline(tempfile, line); 
            PEKS_Enc(enc_kwdb[i], int_kws[stoi(line)], MPKD); // Generate the searchable ciphertext for each subject based on the keyword int_kws (from the line in the file)
        }
        tempfile.close();
    }
    
    //cout << "Encrypted keyword datbase setup for the subjects are completed" << endl << endl; 
}

void Single_Cluster_Keyword_Search(long int (&candidate_list)[nb_subjects], long int(&enc_kwdb)[nb_subjects][3][N0], CC_t (&cluster_trapdoors)[nb_kw_clusters][N0], int keyword_cluster, MPK_Data * MPKD, MSK_Data * MSKD) {
    // Finds all the candidates from the database which have the same keywords as a probe keyword cluster
    // candidate_list: stores the candidates which have the same keywords as the identified keyword cluster (from a probe)
    // enc_kwdb: stores the searchable ciphertext for each subject
    // cluster_trapdoors: stores the trapdoor for each keyword cluster
    // keyword_cluster: the index of the identified keyword cluster (from a probe)
    // subjects_db: stores the subject id's for all the subjects --> in ascending order (same order as the keyword cluster file is defined) 
    
    unsigned int i, j, index;
    bool frgc = true; // If only frgc database is used 
    
    for (j=0; j<nb_subjects; j++) {
        candidate_list[j] = 0; // Initialize the candidate list as empty
    }
    
    index = 0;
    
    if (frgc) { // Only to return the correct identifier corresponding to the use of frgc database, i.e., return number between 529 and 1061
        i = 0;
        index = (nb_total_subjects-nb_subjects);
        for (j=0; j<nb_subjects; j++) {
            if (PEKS_Test(enc_kwdb[j], cluster_trapdoors[keyword_cluster]) == 1) { // Find all the subjects with the same keywords as the indicated keyword_cluster
                candidate_list[i] = index;
                i++;
            }
            index++;
        }
    }
    else { // if use feret database
        i = 0; 
        for (j=0; j<nb_subjects; j++) {
            if (PEKS_Test(enc_kwdb[j], cluster_trapdoors[keyword_cluster]) == 1) { // Find all the subjects with the same keywords as the indicated keyword_cluster
                candidate_list[i] = j;
                i++;
            }
        }
    }
}

void Generate_Candidate_List_File(long int(&enc_kwdb)[nb_subjects][3][N0], CC_t (&cluster_trapdoors)[nb_kw_clusters][N0], int keyword_cluster, MPK_Data * MPKD, MSK_Data * MSKD) {
    // Creates candidatelist file for the candidates associated with the specified keyword_cluster
    // enc_kwdb: stores the searchable ciphertext for each subject
    // cluster_trapdoors: stores the trapdoor for each keyword cluster
    // keyword_cluster: the index of the identified keyword cluster (from a probe)
    // subjects_db: stores the subject id's for all the subjects --> in ascending order (same order as the keyword cluster file is defined) 
    
    unsigned int i; 
    long int candidate_list[nb_subjects];
    ofstream CandidateFile("/users/cecilie/documents/code/identification-PEKS/frgc/candidatelist" + to_string(keyword_cluster) + ".txt"); 
    //ofstream CandidateFile("/users/cecilie/documents/code/identification-PEKS/feret/candidatelist" + to_string(keyword_cluster) + ".txt"); 

    Single_Cluster_Keyword_Search(candidate_list, enc_kwdb, cluster_trapdoors, keyword_cluster, MPKD, MSKD); // Finds all the candidates from the database which have the same keywords as the specified keyword_cluster

    for (i=0; i<nb_subjects; i++) {
        if (candidate_list[i] !=0) {
            CandidateFile << candidate_list[i] << endl; 
        }
    }
    
    CandidateFile.close();
    
    //cout << "\n\nGeneration of the candidate list file for keyword cluster " << keyword_cluster << " is successfully computed" << endl << endl; 
}

void Subjects_Database(long int (&subjects_db)[nb_total_subjects]) {
    // Create an array consisting of the subject id's from the features 
    // subjects_db: resulting array consisting of subject id's to nb_subjects subjects --> in ascending order (same order as the keyword cluster list is defined) 
    
    string frgc_dir_path = "/users/cecilie/documents/code/identification/features/frgc_arcface512_float/frgc_arcface512_npy";
    string feret_dir_path = "/users/cecilie/documents/code/identification/features/feret_arcface512_float/feret_arcface512_npy"; 
    unsigned int i, j;
    vector<string> frgc_files; 
    vector<string> feret_files;
    ofstream MappingFile("/users/cecilie/documents/code/identification-PEKS/MappingSubjects.txt"); 

    // FERET mapping 
    for (const auto & file: fs::directory_iterator(feret_dir_path)) {
        fs::path path(file.path()); 
        string subject = path.filename().string().substr(0, 5);
        feret_files.push_back(subject); 
    }
    
    sort(feret_files.begin(), feret_files.end()); // Sort the vector of subject id's in ascending order 
    feret_files.erase(unique(feret_files.begin(), feret_files.end()), feret_files.end()); // Remove duplicates from the vector 
    
    MappingFile << "Mapping of subject id (0-1061) to corresponding subject id in FERET and FRGC database" << endl;
    MappingFile << "FERET: 0-528, FRGC: 529-1061" << endl << endl;
    
    for (i = 0; i<529; i++) {
        subjects_db[i] = stoi(feret_files[i]); // Initialize the resulting subject id array
        MappingFile << i << " : " << subjects_db[i] << endl;
    }
    
    // FRGC mapping
    for (const auto & file: fs::directory_iterator(frgc_dir_path)) {
        fs::path path(file.path());
        string subject = path.filename().string().substr(0, 5); 
        frgc_files.push_back(subject);
    }
    
    sort(frgc_files.begin(), frgc_files.end()); // Sort the vector of subject id's in ascending order 
    frgc_files.erase(unique(frgc_files.begin(), frgc_files.end()), frgc_files.end()); // Remove duplicates from the vector 
    
    for (j = 0; j<533; j++) {
        subjects_db[i] = stoi(frgc_files[j]); // Initialize the resulting subject id array
        MappingFile << i << " : " << subjects_db[i] << endl;
        i++;
        
    }
    
    MappingFile.close();
    
    //cout << "Subject id array is successfully generated" << endl << endl;
}

void Biometric_Test_With_Cluster_Keywords(CC_t (&cluster_trapdoors)[nb_kw_clusters][N0], MPK_Data * MPKD, MSK_Data * MSKD){ 
    // Biometric test function with the keyword cluster approach 
    // cluster_trapdoors: stores the trapdoor for each keyword cluster
    
    unsigned int i, j, keyword_cluster; 
    long int Probe[3][N0], int_kws[nb_kw_clusters][N0], enc_kwdb[nb_subjects][3][N0], candidate_list[nb_subjects];
    long int subjects_db[nb_total_subjects];
    int found_match = 0, nb_candidates = 0; 
    
    Subjects_Database(subjects_db); 
    cout << "Subject id array is successfully generated" << endl << endl;

    Generate_Cluster_Trapdoors(cluster_trapdoors, int_kws, MSKD); 
    cout << "Generation of trapdoors for each keyword cluster are completed" << endl << endl; 
    
    Encrypted_Cluster_Keyword_Database(enc_kwdb, int_kws, MPKD, MSKD);
    cout << "Encrypted keyword datbase setup for the subjects are completed" << endl << endl;
    
    PEKS_Enc(Probe, int_kws[12], MPKD); // Generate the searchable ciphertext for the probe with the keyword cluster int_kws[] of the probe
    
    for (i=0; i<nb_kw_clusters; i++) {
        if (PEKS_Test(Probe, cluster_trapdoors[i]) == 1) { //Find the corresponding trapdoor to the probe 
            found_match = 1; 
            keyword_cluster = i; 
            cout << "\nFor the probe, PEKS_Test == 1 for subject i = " << i << endl << endl; // Index i is the trapdoor index corresponding to the probe trapdoor
        }
    }
    
    if (found_match == 0) {
        cout << "No subjects with this cluster center as keyword is enrolled" << endl; 
    }
    
    Single_Cluster_Keyword_Search(candidate_list, enc_kwdb, cluster_trapdoors, keyword_cluster, MPKD, MSKD); // Find all the candidate subjects with the same keywords as the probe

    
    cout << "Subjects in the candidate list are: ";
    for (j=0; j<nb_subjects; j++) {
        if (candidate_list[j] != 0) {
            nb_candidates ++; 
            cout << candidate_list[j] << ", ";
        }
    }
    
    /*// Generate all the candidate lists for every keyword cluster
    for (i=0; i<nb_kw_clusters; i++) {
        //Generate_Candidate_List_File(enc_kwdb, cluster_trapdoors, i, subjects_db, MPKD, MSKD);
        Generate_Candidate_List_File(enc_kwdb, cluster_trapdoors, i, MPKD, MSKD);
    }
    */
    
    Generate_Candidate_List_File(enc_kwdb, cluster_trapdoors, keyword_cluster, MPKD, MSKD); // Generate candidate list for the keyword cluster associated with the probe 
    cout << "\n\nGeneration of the candidate list file for keyword cluster " << keyword_cluster << " is successfully computed" << endl << endl; 


    cout << "There are " << nb_candidates << " subjects that are similar to the probe and added to the candidate list" << endl; 
    
    cout << "\n\nBiometric test with cluster centers are successfully performed" << endl; 
}

void Biometric_Bench_With_Cluster_Keywords(CC_t (&cluster_trapdoors)[nb_kw_clusters][N0], MPK_Data * MPKD, MSK_Data * MSKD){
    // Measure the time for the functionality in the Biometric test function with the keyword cluster approach 
    // cluster_trapdoors: stores the trapdoor for each keyword cluster
    
    unsigned int i, j, keyword_cluster;
    long int Probe[3][N0], int_kws[nb_kw_clusters][N0], enc_kwdb[nb_subjects][3][N0], candidate_list[nb_subjects]; 
    long int subjects_db[nb_total_subjects];
    int found_match = 0;
    clock_t time1, time2, time3, time4, time5, time6, time7, time8, time9, time10, time11, time12; 
    float diff1, diff2, diff3, diff4, diff5, diff6; 
    
    time1 = clock();
    Subjects_Database(subjects_db); 
    time2 = clock();
    diff1 = ((float)time2 - (float)time1)/1000000.01;
    cout << "Subject id array is successfully generated" << endl << endl;
    cout << "Creatimg the subject id array took " << diff1 << " seconds for " << nb_subjects << " subjects" << endl << endl << endl;
    
    time3 = clock(); 
    Generate_Cluster_Trapdoors(cluster_trapdoors, int_kws, MSKD);
    time4 = clock(); 
    diff2 = ((float)time4 - (float)time3)/1000000.01;
    cout << "Generation of trapdoors for each keyword cluster are completed" << endl << endl; 
    cout << "Generating the trapdoors took " << diff2 << " seconds for " << nb_kw_clusters << " keyword cluster classes" << endl << endl << endl; 
    
    time5 = clock(); 
    Encrypted_Cluster_Keyword_Database(enc_kwdb, int_kws, MPKD, MSKD);
    time6 = clock();
    diff3 = ((float)time6 - (float)time5)/1000000.01;
    cout << "Encrypted keyword datbase setup for the subjects are completed" << endl << endl;
    cout << "Generating the searchable ciphertexts for the subjects and setup the encrypted keyword database took " << diff3 << " seconds for " << nb_subjects << " subjects" << endl << endl;

    time7 = clock();
    PEKS_Enc(Probe, int_kws[44], MPKD); 
    time8 = clock();
    diff4 = ((float)time8 - (float)time7)/1000000.01;
    cout << "Generating the searchable ciphertexts for the probe took " << diff4 << " seconds" << endl << endl;

    time9 = clock(); 
    for (i=0; i<nb_kw_clusters; i++){
        if(PEKS_Test(Probe, cluster_trapdoors[i]) == 1){
            found_match = 1; 
            keyword_cluster = i;
            cout << "\nFor the probe, PEKS_Test == 1 for subject i = " << i << endl;

        }
    }
    time10 = clock();
    diff5 = ((float)time10 - (float)time9)/1000000.01;
    
    if (found_match == 0) {
        cout << "\nNo subjects with this cluster center as keyword is enrolled" << endl; 
    }
    
    cout << "\nKeyword search to find corresponding probe trapdoor took " << diff5 << " seconds for " << nb_kw_clusters << " keyword cluster classes" << endl << endl; 
    
    time11 = clock();
    Single_Cluster_Keyword_Search(candidate_list, enc_kwdb, cluster_trapdoors, keyword_cluster, MPKD, MSKD); 
    time12 = clock();
    diff6 = ((float)time12 - (float)time11)/1000000.01;
    
    cout << "Keyword cluster search to find the candidates with the same keywords as the probe took " << diff6 << " seconds" << endl << endl;
    
    cout << "Subjects in the candidate list are: ";
    for (j=0; j<nb_subjects; j++) {
        if (candidate_list[j] != 0) {
            cout << candidate_list[j] << ", ";
        }
    }
    cout << endl << endl;
}

void Cluster_Trapdoor_Bench(const unsigned int nb_executions, MSK_Data * MSKD) {
    clock_t time1, time2, time3, time4; 
    float diff1, diff2; 
    unsigned int i, j; 
    CC_t cluster_trapdoors[N0];
    long int int_kws[N0]; 
    vec_ZZ kws;
    ZZX temp_trapdoors[2];   
    float execution[nb_executions], median = 0;  
    ofstream ExecutionFile("/users/cecilie/documents/code/identification-PEKS/execution/frgc_" + to_string(nb_executions) + "_trapdoor_bench.txt"); 
    
    time1 = clock(); 
    cout << "0%" << flush;
    for (i=0; i<nb_executions; i++) {
        time2 = clock(); 
        kws = RandomVector(); 
        
        PEKS_Trapdoor(temp_trapdoors, kws, MSKD); 
        PEKS_Verify_Trapdoor(temp_trapdoors, kws, MSKD);  
        ZZXToFFT(cluster_trapdoors, temp_trapdoors[1]); 
        
        for (j=0; j<N0; j++) {
            int_kws[j] = conv<long int>(kws[j]);
        }
        
        if((i+1)%(nb_executions/10)==0){
            cout << "..." << (i+1)/(nb_executions/10) << "0%" << flush;
        }
        time3 = clock();
        diff1 = ((float)time3 - (float)time2)/1000000.01;
        execution[i] = diff1;
        ExecutionFile << diff1 << endl;
        
    }
    time4 = clock(); 
    diff2 = ((float)time4 - (float)time1)/1000000.01;
    
    sort(execution, execution + nb_executions);
    if (nb_executions % 2 != 0) {
        median = execution[nb_executions/2];
    }
    else {
        median = (execution[(nb_executions-1)/2] + execution[nb_executions/2])/2.0;
    }
    
    cout << "\n\nIt took " << diff2 << " seconds to create " << nb_executions << " cluster trapdoors." << endl;
    cout << "That's " << (diff2/nb_executions)*1000 << " milliseconds per cluster trapdoor." << endl << endl;
    cout << "The median time for one trapdoor generation is: " << median << " seconds" << endl << endl;
    
    ExecutionFile << "\nFor " << nb_executions << " executions of trapdoor generation, the median time is: ";
    ExecutionFile << median << " seconds" << endl << endl;

}

void Cluster_Trapdoor_List_Bench(const unsigned int nb_executions, MSK_Data * MSKD) {
    unsigned int i, j, k;
    vec_ZZ kws[nb_kw_clusters];
    ZZX temp_trapdoors[nb_kw_clusters][2]; 
    CC_t cluster_trapdoors[nb_kw_clusters][N0];
    long int int_kws[nb_kw_clusters][N0]; 
    
    clock_t time1, time2, time3, time4, time5, time6; 
    float diff1, diff2, diff3; 
    float execution_trapdoor_list[nb_executions], median_trapdoor_list = 0;
    float execution_trapdoors[nb_executions*nb_kw_clusters], median_trapdoors = 0;
    ofstream ExecutionFile("/users/cecilie/documents/code/identification-PEKS/execution/frgc_" + to_string(nb_executions) + "_trapdoor_list_bench.txt");  
    
    time1 = clock(); 
    cout << "0%" << flush;
    for (i=0; i<nb_executions; i++) {
        time2 = clock(); 
        for (j=0; j<nb_kw_clusters; j++) {
            time3 = clock();
            kws[j] = RandomVector(); 
            
            PEKS_Trapdoor(temp_trapdoors[j], kws[j], MSKD); 
            PEKS_Verify_Trapdoor(temp_trapdoors[j], kws[j], MSKD); 
            ZZXToFFT(cluster_trapdoors[j], temp_trapdoors[j][1]); 
            
            for (k=0; k<N0; k++) {
                int_kws[j][k] = conv<long int>(kws[j][k]); 
            }
            time4 = clock(); 
            diff1 = ((float)time4 - (float)time3)/1000000.01;
            execution_trapdoors[i] = diff1;
            
        }
        if((i+1)%(nb_executions/10)==0){
            cout << "..." << (i+1)/(nb_executions/10) << "0%" << flush;
        }
        time5 = clock();
        diff2 = ((float)time5 - (float)time2)/1000000.01;
        execution_trapdoor_list[i] = diff2;
        ExecutionFile << diff2 << endl;
    }
    time6 = clock(); 
    diff3 = ((float)time6 - (float)time1)/1000000.01;
    
    sort(execution_trapdoor_list, execution_trapdoor_list + nb_executions);
    sort(execution_trapdoors, execution_trapdoors + nb_executions);
    if (nb_executions % 2 != 0) {
        median_trapdoor_list = execution_trapdoor_list[nb_executions/2];
        median_trapdoors = execution_trapdoors[nb_executions/2];
    }
    else {
        median_trapdoor_list = (execution_trapdoor_list[(nb_executions-1)/2] + execution_trapdoor_list[nb_executions/2])/2.0;
        median_trapdoors = (execution_trapdoors[(nb_executions-1)/2] + execution_trapdoors[nb_executions/2])/2.0;
    }
    
    cout << "\n\nIt took " << diff3 << " seconds to create " << nb_executions << " cluster trapdoor list of " << nb_kw_clusters << " clusters" << endl;
    cout << "That's " << (diff3/nb_executions)*1000 << " milliseconds per cluster trapdoor list generation." << endl << endl;
    cout << "The median time per cluster trapdoor list generation is: " << median_trapdoor_list << " seconds" << endl;
    cout << "The median time for one cluster trapdoor generation is: " << median_trapdoors << " seconds" << endl << endl;
    
    ExecutionFile << "\nFor " << nb_executions << " executions of trapdoor list generation with " << nb_kw_clusters << " cluster keywords, the median time is: ";
    ExecutionFile << median_trapdoor_list << " seconds" << endl;
    ExecutionFile << "\nThe median time per one trapdoor generation is: " << median_trapdoors << " seconds" << endl << endl;
}

void Encrypted_Searchable_Ciphertext_Bench(const unsigned int nb_executions, MPK_Data * MPKD, MSK_Data * MSKD) {
    clock_t time1, time2, time3, time4; 
    float diff1, diff2; 
    unsigned int i, j; 
    CC_t cluster_trapdoors[N0];
    long int int_kws[N0]; 
    vec_ZZ kws;
    ZZX temp_trapdoors[2];   
    float execution[nb_executions], median = 0;
    long int searchable_ciphertext[3][N0];
    ofstream ExecutionFile("/users/cecilie/documents/code/identification-PEKS/execution/frgc_" + to_string(nb_executions) + "_searchable_ciphertext_bench.txt");  

   
    kws = RandomVector(); 
        
    PEKS_Trapdoor(temp_trapdoors, kws, MSKD); 
    PEKS_Verify_Trapdoor(temp_trapdoors, kws, MSKD); 
    ZZXToFFT(cluster_trapdoors, temp_trapdoors[1]); 
    
    for (j=0; j<N0; j++) {
        int_kws[j] = conv<long int>(kws[j]);
    }
    
    time1 = clock(); 
    cout << "0%" << flush;

    for (i=0; i<nb_executions; i++) {
        time2 = clock(); 
        
        PEKS_Enc(searchable_ciphertext, int_kws, MPKD); 

        if((i+1)%(nb_executions/10)==0){
            cout << "..." << (i+1)/(nb_executions/10) << "0%" << flush;
        }
        time3 = clock();
        diff1 = ((float)time3 - (float)time2)/1000000.01;
        execution[i] = diff1;
        ExecutionFile << diff1 << endl;
        
    }
    time4 = clock(); 
    diff2 = ((float)time4 - (float)time1)/1000000.01;
    
    sort(execution, execution + nb_executions);
    if (nb_executions % 2 != 0) {
        median = execution[nb_executions/2];
    }
    else {
        median = (execution[(nb_executions-1)/2] + execution[nb_executions/2])/2.0;
    }
    
    cout << "\n\nIt took " << diff2 << " seconds to create " << nb_executions << " searchable ciphertexts." << endl;
    cout << "That's " << (diff2/nb_executions)*1000 << " milliseconds per searchable ciphertext." << endl << endl;
    cout << "The median time for one searchable ciphertext generation is: " << median << " seconds" << endl << endl;
    
    ExecutionFile << "\nFor " << nb_executions << " executions of encryption of searchable ciphertext, the median time is: ";
    ExecutionFile << median << " seconds" << endl << endl;
}

void Reverse_PEKS_Bench(const unsigned int nb_executions, MPK_Data * MPKD, MSK_Data * MSKD) {
    clock_t time1, time2, time3, time4; 
    float diff1, diff2; 
    unsigned int i, j; 
    CC_t cluster_trapdoors[nb_kw_clusters][N0];
    long int int_kws[nb_kw_clusters][N0]; 
    vec_ZZ kws[nb_kw_clusters];
    ZZX temp_trapdoors[nb_kw_clusters][2];   
    float execution[nb_executions], median = 0;
    long int searchable_ciphertext[3][N0];
    bool found_match;
    ofstream ExecutionFile("/users/cecilie/documents/code/identification-PEKS/execution/frgc_" + to_string(nb_executions) + "_reverse_peks_search_bench.txt");  
    
    
    for (i=0; i<nb_kw_clusters; i++) {
        kws[i] = RandomVector(); 
        
        PEKS_Trapdoor(temp_trapdoors[i], kws[i], MSKD); 
        PEKS_Verify_Trapdoor(temp_trapdoors[i], kws[i], MSKD); 
        ZZXToFFT(cluster_trapdoors[i], temp_trapdoors[i][1]); 
        
        for (j=0; j<N0; j++) {
            int_kws[i][j] = conv<long int>(kws[i][j]);
        }
        
    }

    time1 = clock(); 
    
    cout << "0%" << flush;

    for (i=0; i<nb_executions; i++) {
        time2 = clock(); 
        
        int random_kw = rand() % 64;
        
        PEKS_Enc(searchable_ciphertext, int_kws[random_kw], MPKD); 
        
        for (j=0; j<nb_kw_clusters; j++){
            if(PEKS_Test(searchable_ciphertext, cluster_trapdoors[j]) == 1){
                found_match = true; 
            }
        }
        
        if (!found_match) {
            cout << "PEKS test failed!" << endl; 
            break; 
        }
    
        if((i+1)%(nb_executions/10)==0){
            cout << "..." << (i+1)/(nb_executions/10) << "0%" << flush;
        }
        time3 = clock();
        diff1 = ((float)time3 - (float)time2)/1000000.01;
        execution[i] = diff1;
        ExecutionFile << diff1 << endl;
    }
    time4 = clock(); 
    diff2 = ((float)time4 - (float)time1)/1000000.01;
    
    sort(execution, execution + nb_executions);
    if (nb_executions % 2 != 0) {
        median = execution[nb_executions/2];
    }
    else {
        median = (execution[(nb_executions-1)/2] + execution[nb_executions/2])/2.0;
    }
    
    cout << "\n\nIt took " << diff2 << " seconds to perform " << nb_executions << " reverse PEKS search." << endl;
    cout << "That's " << (diff2/nb_executions)*1000 << " milliseconds per PEKS test." << endl << endl;
    cout << "The median time for one PEKS test is: " << median << " seconds" << endl << endl;
    
    ExecutionFile << "\nFor " << nb_executions << " executions of reverse PEKS search, the median time is: ";
    ExecutionFile << median << " seconds" << endl << endl;
}