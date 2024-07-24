/******************************************************************************
  openmx2xatu.c:
  =============

    openmx2xatu.x extract the Kohn-Sham Hamiltonian from filename.scfout, and 
    write it in xatu input format

    N. N. Batista
    nnardoto@gmail.com
    17/Jun./2024
    
    Adaptation of analysis_example.c from openmx@3.9 code
    https://www.openmx-square.org/	

  
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "read_scfout.h"

//#include "mpi.h"
#define ToEV 27.2114

void ExtractHamiltonian(FILE *fp, double ****RH);
void ExtractOverlap(FILE *fp, double ****Overlap);


void main(int argc, char *argv[]) 
{
  // read scfout from openmx@3.9
  read_scfout(argv);


  FILE *fHS = fopen("H.dat", "w");
  FILE *fSS = fopen("S.dat", "w");

  // Hks[0]: Hamiltonian nonSpin Case
  // OLP   : Overlap Matrices

  ExtractHamiltonian(fHS, Hks[0]);
  ExtractOverlap(fSS, OLP);

  fclose(fHS);
  fclose(fSS);

}


void ExtractOverlap(FILE *fp, double ****Overlap)
{
  static int i,j,tnoA,tnoB,Anum,Bnum,NUM,GA_AN,LB_AN,GB_AN;
  static int l1,l2,l3,Rn;
  int MP[atomnum];

  Anum = 1;
  for (i=1; i<=atomnum; i++){
    MP[i] = Anum;
    Anum += Total_NumOrbs[i];
  }
  NUM = Anum - 1;

  /****************************************************
                       set overlap
  ****************************************************/
  for (GA_AN=1; GA_AN<=atomnum; GA_AN++)
  {
    tnoA = Total_NumOrbs[GA_AN];
    Anum = MP[GA_AN];

    for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++)
    {
      GB_AN = natn[GA_AN][LB_AN];
      Rn = ncn[GA_AN][LB_AN];
      tnoB = Total_NumOrbs[GB_AN];

      l1 = atv_ijk[Rn][1];
      l2 = atv_ijk[Rn][2];
      l3 = atv_ijk[Rn][3];

      Bnum = MP[GB_AN];
      for (i=0; i<tnoA; i++)
      {
   	   for (j=0; j<tnoB; j++)
   	   {
	       fprintf(fp, "% 2i  % 2i  % 2i  % 4i  % 4i  % 10.7f\n", l1, l2, l3, Anum + i, Bnum + j, OLP[GA_AN][LB_AN][i][j]);
   	   }
      }
    }
  }
}

void ExtractHamiltonian(FILE *fp, double ****RH)
{
  static int i,j,tnoA,tnoB,Anum,Bnum,NUM,GA_AN,LB_AN,GB_AN;
  static int l1,l2,l3,Rn;
  int MP[atomnum];

  Anum = 1;
  for (i=1; i<=atomnum; i++)
  {
    MP[i] = Anum;
    Anum += Total_NumOrbs[i];
  }
  NUM = Anum - 1;


  /****************************************************
                    set Hamiltonian
  ****************************************************/

  for (GA_AN=1; GA_AN<=atomnum; GA_AN++)
  {
    tnoA = Total_NumOrbs[GA_AN];
    Anum = MP[GA_AN];

    for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++)
    {
      GB_AN = natn[GA_AN][LB_AN];
      Rn = ncn[GA_AN][LB_AN];
      tnoB = Total_NumOrbs[GB_AN];

      l1 = atv_ijk[Rn][1];
      l2 = atv_ijk[Rn][2];
      l3 = atv_ijk[Rn][3];

      Bnum = MP[GB_AN];
      for (i=0; i<tnoA; i++)
      {
	      for (j=0; j<tnoB; j++)
	      {
	        fprintf(fp, "% 2i  % 2i  % 2i  % 4i  % 4i  % 10.7f\n", l1, l2, l3, Anum + i, Bnum + j, ToEV*RH[GA_AN][LB_AN][i][j]);
	      } 
      } 
    }
  }
}
