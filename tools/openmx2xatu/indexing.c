/******************************************************************************
  indexing.c:
  =============

    reindexing of Kohn-Sham Hamiltonian from filename.scfout, for make interfaces
    e.g Xatu code.

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
#include "indexing.h"


int** aIndex;
int Index(int l, int m, int n)
{
    return Map[N + l][N + m] - 1;
}

void Clear()
{
    // verificar se é só isso mesmo
    for(int i = 0; i < 2*N + 1; i++)
    {
        free(Map[i]);
    }
    free(Map);
}

void MakeMap()
{
    N = 0; 
    for (int GlobalA = 1; GlobalA <= atomnum; GlobalA++)
    {
        for (int LocalB = 0; LocalB <= FNAN[GlobalA]; LocalB++)
        {
            int Rn = ncn[GlobalA][LocalB];
      
            for(int i = 0; i < 3; i++)
            {
                if(N < abs(atv_ijk[Rn][i + 1]))
                {
                    N = abs(atv_ijk[Rn][i + 1]);
                }
            }
        }
    } 

    // Alloca memoria para indexação
    Map = malloc((2*N + 1) * sizeof(int*));
    for(int i = 0; i < 2*N + 1; i++)
    {
        Map[i] = malloc((2*N + 1) * sizeof(int));
        for(int j = 0; j < 2*N + 1; j++)
        {
            Map[i][j] = 0;
        }
    }

    // Verifica quais matrizes existem
    nMat = 0;
    for (int GlobalA = 1; GlobalA <= atomnum; GlobalA++)
    {
        for (int LocalB = 0; LocalB <= FNAN[GlobalA]; LocalB++)
        {
            int Rn = ncn[GlobalA][LocalB];
            int l, m, n;
            l = atv_ijk[Rn][1];
            m = atv_ijk[Rn][2];
            
            Map[N + l][N + m] = 1;
        }
    } 

    // indexa corretamente
    for(int i = -N; i <= N; i++)
    {
        for(int j = -N; j <= N; j++)
        {
            if(Map[N + i][N + j])
            {
                nMat += 1;
                Map[N + i][N + j] = nMat;
            }
        }
    }

    // preenche a indexação reversa
    aIndex = malloc(nMat * sizeof(int*));
    for(int i = 0; i < nMat; i++)
    {
        aIndex[i] = malloc(3 * sizeof(int));
        aIndex[i][0] = 0;
        aIndex[i][1] = 0;
        aIndex[i][2] = 0;
    }
}

void Export2Xatu(char* SystemName)
{
    static int i,j,tnoA,tnoB,Anum,Bnum,NUM,GA_AN,LB_AN,GB_AN;
    static int l1,l2,l3,Rn;
    int MP[atomnum];
    int Dimension = 2;

    double ***H, ***S;
    
    Anum = 1;
    for (i=1; i<=atomnum; i++){
        MP[i] = Anum;
        Anum += Total_NumOrbs[i];
    }
    NUM = Anum - 1;
    MSize = NUM;

    // Aloca memoria
    H = malloc(nMat * sizeof(double*));
    S = malloc(nMat * sizeof(double*));
    for(int nn = 0; nn < nMat; nn++)
    {
        H[nn] = malloc(NUM * sizeof(double*));
        S[nn] = malloc(NUM * sizeof(double*));    
        for(int i = 0; i < NUM; i++)
        {
            H[nn][i] = malloc(NUM * sizeof(double));
            S[nn][i] = malloc(NUM * sizeof(double));
        }
    }

    // inicializa com zeros 
    for(int nn = 0; nn < nMat; nn++)
    {
        for(int i = 0; i < NUM; i++)
        {
            for(int j = 0; j < NUM; j++)
            {
                H[nn][i][j] = 0.0;
                S[nn][i][j] = 0.0;
            }
        }
    }


    // Reindex Hamiltonian
    ExtractHamiltonian(H, Hks[0]);
    ExtractOverlap(S, OLP);

    // create xatu file
    FILE* fp = fopen("SystemName.model", "w");

    // Set head of xatu input model
    fprintf(fp, "# dimension\n");
    fprintf(fp, "2\n");

    fprintf(fp, "# norbitals\n");
    for(int k = 1; k < atomnum + 1; k++)
    {
        fprintf(fp, "%d ", Total_NumOrbs[k]);
    }
    fprintf(fp, "\n");

    fprintf(fp, "# bravaislattice\n");
    for(int k = 1; k < Dimension + 1; k++)
    {
        fprintf(fp, "% 12.8lf    % 12.8lf    % 12.8lf\n", ToAng*tv[k][1], ToAng*tv[k][2], ToAng*tv[k][3]);
    }

    fprintf(fp, "# motif\n");
    for(int k = 1; k < atomnum + 1; k++)
    {
        fprintf(fp, "% 12.8lf    % 12.8lf    %12.8lf    %12.8lf\n", ToAng*Gxyz[k][1], ToAng*Gxyz[k][2], ToAng*Gxyz[k][3], (double) k - 1);
    }

    // Acesss the bravais vector 
    fprintf(fp, "# bravaisvectors\n");
    for(int nn = 0; nn < nMat; nn++)
    {
        int i = aIndex[nn][0];
        int j = aIndex[nn][1];

        double xx = i * tv[1][1] + j*tv[2][1];
        double yy = i * tv[1][2] + j*tv[2][2];

        fprintf(fp, "% 12.8lf    % 12.8lf    % 12.8lf\n", ToAng*xx, ToAng*yy, 0.0);
    }

    // Access at same order of bravais vectors for constructions of hamiltonian
    fprintf(fp, "# hamiltonian\n");
    for(int nn = 0; nn < nMat; nn++)
    {
        for(int i = 0; i < NUM; i++)
        {
            for(int j = 0; j < NUM; j++)
            {
                fprintf(fp, "% 12.8lf  % 12.8lfj    ", H[nn][i][j], 0.0);
            }
        
            fprintf(fp, "\n");
        }
        fprintf(fp, "&\n");
    }    

    fprintf(fp, "# overlap\n");
    for(int nn = 0; nn < nMat; nn++)
    {
        for(int i = 0; i < NUM; i++)
        {
            for(int j = 0; j < NUM; j++)
            {
                fprintf(fp, "% 12.8lf  % 12.8lfj    ", S[nn][i][j], 0.0);
            }
        
            fprintf(fp, "\n");
        }
        fprintf(fp, "&\n");
    }

    fprintf(fp, "# filling\n");
    fprintf(fp, "%i\n", Valence_Electrons/2);
    fprintf(fp, "#");

    fclose(fp);
}

void ExtractOverlap(double ***NewOverlap, double ****Overlap)
{
    static int i,j,tnoA,tnoB,Anum,Bnum,NUM,GA_AN,LB_AN,GB_AN;
    static int l1,l2,l3,Rn;
    int MP[atomnum];

    // Size of Matrices
    Anum = 1;
    for (i = 1; i <= atomnum; i++)
    {
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
                    NewOverlap[Index(l1, l2, l3)][Anum + i - 1][Bnum + j - 1] = OLP[GA_AN][LB_AN][i][j];
   	            }
            }
        }
    }
}

void ExtractHamiltonian(double ***NewHH, double ****RH)
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
                    NewHH[Index(l1, l2, l3)][Anum + i - 1][Bnum + j - 1] = ToEV*RH[GA_AN][LB_AN][i][j];
                    aIndex[Index(l1, l2, l3)][0] = l1;
                    aIndex[Index(l1, l2, l3)][1] = l2;
                    aIndex[Index(l1, l2, l3)][2] = l3;
                } 
            } 
        }
    }
}
