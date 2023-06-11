#ifndef LIBE_SCHEME_H
#define LIBE_SCHEME_H

#include "params.h"
#include "Sampling.h"
#include "FFT.h"
#include "Random.h"
#include "Algebra.h"

void Keygen(ZZ_pX& PublicKey, ZZX* PrivateKey);
void CompletePrivateKey(mat_ZZ& B, const ZZX * const PrivateKey);
void GPV(RR_t * v, const RR_t * const c, const RR_t s, const MSK_Data * const MSKD);
void CompleteMSK(MSK_Data * MSKD, ZZX * MSK);
void CompleteMPK(MPK_Data * MPKD, ZZ_pX MPK);
void IBE_Extract(ZZX SK_id[2], vec_ZZ id, const MSK_Data * const MSKD);
unsigned long IBE_Verify_Key(const ZZX SK_id[2], const vec_ZZ id, const MSK_Data * const MSKD);
void PEKS_Enc(long C[3][N0],  const long id0[N0],  const MPK_Data * const MPKD);
bool PEKS_Test( const long C[3][N0], const CC_t * const SKtd_FFT);

void Trapdoor_Bench(const unsigned int nb_extr, MSK_Data * MSKD);
void Encrypt_Bench(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD);
void Trapdoor_Test(const unsigned int nb_extr, MSK_Data * MSKD);
void Encrypt_Test(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD);
void Keyword_Database(CC_t (&keys)[8][N0], long int (&kwdb)[nb_kwcl][nb_sub][3][N0], MPK_Data * MPKD, MSK_Data * MSKD);
void Single_Keyword_Search(int kw_class, int keyword, long int (&kwdb)[nb_kwcl][nb_sub][3][N0], CC_t (&keys)[8][N0], long int (&result)[nb_sub], MPK_Data * MPKD, MSK_Data * MSKD);
void Cascading_Keyword_Search(int kw_class0, int keyword0, int kw_class1, int keyword1, int kw_class2, int keyword2, int kw_class3, int keyword3, long int (&kwdb)[nb_kwcl][nb_sub][3][N0], CC_t (&keys)[8][N0], long int (&result)[nb_sub], MPK_Data * MPKD, MSK_Data * MSKD);
void Biometric_Test_Old(long int (&kwdb)[nb_kwcl][nb_sub][3][N0], MPK_Data * MPKD, MSK_Data * MSKD);
void Biometric_Bench_Old(long int (&kwdb)[nb_kwcl][nb_sub][3][N0], MPK_Data * MPKD, MSK_Data * MSKD);
void Gen_Trapdoor_List(CC_t (&keys)[nb_classes][N0], long int (&int_kws)[nb_classes][N0], MSK_Data * MSKD);
void Gen_Lookup_Table(long int (&LookupTable)[nb_ne_classes]);
void Biometric_Test(CC_t (&keys)[nb_classes][N0], MPK_Data * MPKD, MSK_Data * MSKD);
void Biometric_Bench(CC_t (&keys)[nb_classes][N0], MPK_Data * MPKD, MSK_Data * MSKD);


//==============================================================================
//My Implementation: Use Stable Hash Codes (Cluster Centers) As Keywords
//==============================================================================
void Generate_Cluster_Trapdoors(CC_t (&cluster_trapdoors)[nb_kw_clusters][N0], long int (&int_kws)[nb_kw_clusters][N0], MSK_Data * MSKD);
void Encrypted_Cluster_Keyword_Database(long int (&enc_kwdb)[nb_sub][3][N0], long int (&int_kws)[nb_kw_clusters][N0], MPK_Data * MPKD, MSK_Data * MSKD);
void Single_Cluster_Keyword_Search(long int (&candidate_list)[nb_sub], long int(&enc_kwdb)[nb_sub][3][N0], CC_t (&cluster_trapdoors)[nb_kw_clusters][N0], int keyword_cluster, MPK_Data * MPKD, MSK_Data * MSKD);
void Generate_Candidate_List_File(long int(&enc_kwdb)[nb_sub][3][N0], CC_t (&cluster_trapdoors)[nb_kw_clusters][N0], int keyword_cluster, MPK_Data * MPKD, MSK_Data * MSKD);
void Subjects_Database(long int (&subjects_db)[nb_sub]);
void Biometric_Test_With_Cluster_Keywords(CC_t (&cluster_trapdoors)[nb_kw_clusters][N0], MPK_Data * MPKD, MSK_Data * MSKD);
void Biometric_Bench_With_Cluster_Keywords(CC_t (&cluster_trapdoors)[nb_kw_clusters][N0], MPK_Data * MPKD, MSK_Data * MSKD);

void Cluster_Trapdoor_Bench(const unsigned int nb_executions, MSK_Data * MSKD);
void Cluster_Trapdoor_List_Bench(const unsigned int nb_executions, MSK_Data * MSKD);
void Encrypted_Searchable_Ciphertext_Bench(const unsigned int nb_executions, MPK_Data * MPKD, MSK_Data * MSKD);
void Reverse_PEKS_Bench(const unsigned int nb_executions, MPK_Data * MPKD, MSK_Data * MSKD);

#endif
