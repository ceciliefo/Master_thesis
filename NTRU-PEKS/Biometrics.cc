/*

Adopted from Thomas Prest's project https://github.com/tprest/Lattice-IBE and edited for this transformation. 
Therefore, please respect its license and requirements.

*/



#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>


#include "params.h"
#include "io.h"
#include "FFT.h"
#include "Sampling.h"
#include "Random.h"
#include "Algebra.h"
#include "Scheme.h"


using namespace std;
using namespace NTL;

long int kwdb[nb_kwcl][nb_sub][3][N0];
long int enc_kwdb[nb_sub][3][N0]; // implemented for stable hash 

CC_t keys[nb_classes][N0] = {};
CC_t cluster_trapdoors[nb_kw_clusters][N0] = {}; // implemented for stable hash


//==============================================================================
//==============================================================================
//                                  MAIN
//==============================================================================
//==============================================================================


int main()
{
    cout << "\n=======================================================================\n";
    cout << "This program is a proof-of concept soft-biometric keyword search using\n";
    cout << "efficient PEKS over lattices. It generates a NTRU lattice of dimension\n";
    cout << "2N and associated modulus q, and perform benches and tests.";
    cout << "\n=======================================================================\n\n";

    ZZX MSK[4];
    ZZ_pX phiq, MPK;
    unsigned int i;
    float diff;
    MSK_Data * MSKD = new MSK_Data;
    MPK_Data * MPKD = new MPK_Data;
    clock_t t1, t2;
    const ZZX phi = Cyclo();

    struct timespec t;
    srand(clock_gettime(CLOCK_MONOTONIC, &t));
    
    cout << "N = " << N0 << endl;
    cout << "q = " << q0 << endl;

    ZZ_p::init(q1);
    zz_p::init(q0);

    phiq = conv<ZZ_pX>(phi);
    ZZ_pXModulus PHI(phiq);


    cout << "\n===================================================================\n KEY GENERATION";
    cout << "\n===================================================================\n";
    t1 = clock();
    for(i=0; i<1; i++)
    {
        Keygen(MPK, MSK);
    }

    CompleteMSK(MSKD, MSK);
    CompleteMPK(MPKD, MPK);

    t2 = clock();
    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << "It took " << diff << " seconds to generate user keys" << endl;

    
    //==============================================================================
    //Key extraction bench and encryption/decryption bench
    //==============================================================================
    const unsigned int nb_trdb = 1000;
    const unsigned int nb_crypb = 1000;
    //const unsigned int nb_decrypb = 1000;

    cout << "\n===================================================================\n RUNNING PEKS BENCH FOR ";
    cout << nb_crypb << " DIFFERENT KEYWORDS\n===================================================================\n";
    //Encrypt_Bench(nb_crypb, MPKD, MSKD);


    cout << "\n===================================================================\n RUNNING TRAPDOOR BENCH FOR ";
    cout << nb_trdb << " DIFFERENT KEYWORDS\n===================================================================\n";
    //Trapdoor_Bench(nb_trdb, MSKD);

   
    ///==============================================================================
    //Key extraction test and encryption/decryption test
    //==============================================================================
    const unsigned int nb_trdt = 100;
    const unsigned int nb_crypt = 100;
    
    
    cout << "\n===================================================================\n CHECKING PEKS VALIDITY FOR ";
    cout << nb_crypt << " DIFFERENT KEYWORDS\n===================================================================\n";
    //Encrypt_Test(nb_crypt, MPKD, MSKD);
    
    cout << "\n===================================================================\n CHECKING TRAPDOOR VALIDITY FOR ";
    cout << nb_trdt << " DIFFERENT KEYWORDS\n===================================================================\n";
	//Trapdoor_Test(nb_trdt, MSKD);
     
    
    cout << "\n===================================================================\n BIOMETRIC KEYWORD SEARCH ";
    cout << "\n===================================================================\n";
    //Biometric_Test(keys, MPKD, MSKD);
    //Biometric_Test_Old(kwdb, MPKD, MSKD);
    
    
    cout << "\n===================================================================\n RUNNING BIOMETRIC KEYWORD SEARCH BENCH";
    cout << "\n===================================================================\n";
    //Biometric_Bench(keys, MPKD, MSKD);
    //Biometric_Bench_Old(kwdb, MPKD, MSKD);
    
    
///==============================================================================
//My Implementation: Use Stable Hash Codes (Cluster Centers) As Keywords
//===============================================================================
    
    //==============================================================================
    //Trapdoor Bench, Encryption Bench and PEKS Test Bench
    //==============================================================================
        
    const unsigned int nb_executions = 1000;
    const unsigned int nb_executions_full_scheme = 10;
    
    cout << "\n===================================================================\n RUNNING CLUSTER TRAPDOOR BENCH AND MEASURING THE MEDIAN TIME\n ";
    cout << "FOR " << nb_executions << " DIFFERENT KEYWORDS\n===================================================================\n";
    Cluster_Trapdoor_Bench(nb_executions, MSKD);
    
    cout << "\n===================================================================\n RUNNING CLUSTER TRAPDOOR LIST BENCH AND MEASURING THE MEDIAN TIME\n ";
    cout << "FOR " << nb_executions_full_scheme << " EXECUTIONS WITH " << nb_kw_clusters << " DIFFERENT KEYWORDS\n===================================================================\n";
    Cluster_Trapdoor_List_Bench(nb_executions_full_scheme, MSKD);

    cout << "\n===================================================================\n RUNNING SEARCHABLE CIPHERTEXT BENCH AND MEASURING THE MEDIAN TIME\n ";
    cout << "FOR " << nb_executions << " SUBJECT KEYWORDS\n===================================================================\n";
    Encrypted_Searchable_Ciphertext_Bench(nb_executions, MPKD, MSKD); 
    
    cout << "\n===================================================================\n RUNNING REVERSE PEKS SEARCH BENCH AND MEASURING THE MEDIAN TIME\n ";
    cout << "FOR " << nb_executions << " DIFFERENT KEYWORDS\n===================================================================\n";
    Reverse_PEKS_Bench(nb_executions, MPKD, MSKD); 
    
    cout << "\n===================================================================\n BIOMETRIC KEYWORD SEARCH WITH CLUSTER CENTERS";
    cout << "\n===================================================================\n";
    Biometric_Test_With_Cluster_Keywords(cluster_trapdoors, MPKD, MSKD);
    
    cout << "\n===================================================================\n RUNNING BIOMETRIC KEYWORD SEARCH BENCH WITH CLUSTER CENTERS";
    cout << "\n===================================================================\n";
    Biometric_Bench_With_Cluster_Keywords(cluster_trapdoors, MPKD, MSKD);
    
    
    free(MSKD);
    free(MPKD);
    return 0;
}