#ifndef TMalign_HwRMSD_h
#define TMalign_HwRMSD_h 1
#include "NWalign.h"
#include "se.h"

double Kabsch_Superpose(double **r1, double **r2, double **xt,
    double **xa, double **ya, int xlen, int ylen, int invmap[],
    int& L_ali, double t[3], double u[3][3], const int mol_type)
{
    L_ali = 0;
    int i, j;
    for (j = 0; j<ylen; j++)
    {
        i = invmap[j];
        if (i >= 0)
        {
            r1[L_ali][0]  = xa[i][0];
            r1[L_ali][1]  = xa[i][1];
            r1[L_ali][2]  = xa[i][2];

            r2[L_ali][0]  = ya[j][0];
            r2[L_ali][1]  = ya[j][1];
            r2[L_ali][2]  = ya[j][2];

            L_ali++;
        }
        else if (i != -1) PrintErrorAndQuit("Wrong map!\n");
    }

    double RMSD = 0;
    Kabsch(r1, r2, L_ali, 1, &RMSD, t, u);
    RMSD = sqrt( RMSD/(1.0*L_ali) );

    for (i=0; i<xlen; i++)
    {
        xt[i][0] = xa[i][0];
        xt[i][1] = xa[i][1];
        xt[i][2] = xa[i][2];
    }
    do_rotation(xa, xt, xlen, t,u);
    return RMSD;
}

/* outfmt_opt>=2 should not parse the sequence alignment */
int HwRMSD_main(double **xa, double **ya, const char *seqx, const char *seqy,
    const char *secx, const char *secy, double t0[3], double u0[3][3],
    double &TM1, double &TM2, double &TM3, double &TM4, double &TM5,
    double &d0_0, double &TM_0, double &d0A, double &d0B, double &d0u,
    double &d0a, double &d0_out, string &seqM, string &seqxA, string &seqyA,
    double &rmsd0, int &L_ali, double &Liden, double &TM_ali,
    double &rmsd_ali, int &n_ali, int &n_ali8, const int xlen, const int ylen,
    const vector<string>&sequence, const double Lnorm_ass,
    const double d0_scale, const int i_opt,
    const int a_opt, const bool u_opt, const bool d_opt, const int mol_type,
    const int outfmt_opt, int *invmap, const int glocal=0, const int iter_opt=1)
{
    /***********************/
    /* allocate memory     */
    /***********************/
    double t[3], u[3][3]; //Kabsch translation vector and rotation matrix
    double **xt;          //for saving the superposed version of r_1 or xtm
    double **r1, **r2;    // for Kabsch rotation
    int minlen = min(xlen, ylen);
    NewArray(&xt, xlen, 3);
    NewArray(&r1, minlen, 3);
    NewArray(&r2, minlen, 3);
    int *invmap_tmp = new int[ylen+1];

    int i, j, i1, i2, L;
    double TM1_tmp,TM2_tmp,TM3_tmp,TM4_tmp,TM5_tmp,TM_ali_tmp;
    string seqxA_tmp,seqyA_tmp,seqM_tmp;
    double rmsd0_tmp;
    int L_ali_tmp,n_ali_tmp,n_ali8_tmp;
    double Liden_tmp;
    double rmsd_ali_tmp;

    /* initialize alignment */
    TM1=TM2=TM1_tmp=TM2_tmp=L_ali=-1;

    if (i_opt)
    {
        seqxA_tmp=sequence[0];
        seqyA_tmp=sequence[1];
    }
    else NWalign_main(seqx, seqy, xlen, ylen, seqxA_tmp, seqyA_tmp,
            mol_type, invmap_tmp, 1, glocal);
    int total_iter=(i_opt==3 || iter_opt<1)?1:iter_opt;

    /*******************************/
    /* perform iterative alignment */
    /*******************************/
    for (int iter=0;iter<total_iter;iter++)
    {
        n_ali_tmp=n_ali8_tmp=0;
        /* get ss alignment for the second iteration */
        if (iter==1 && !i_opt) NWalign_main(secx, secy, xlen, ylen,
            seqxA_tmp, seqyA_tmp, mol_type, invmap_tmp, 1, glocal);

        /* parse initial alignment */
        if (seqxA_tmp.size())
        {
            for (j = 0; j < ylen; j++) invmap_tmp[j] = -1;
            i1 = -1;
            i2 = -1;
            L = min(seqxA_tmp.size(), seqyA_tmp.size());
            for (j = 0; j<L; j++)
            {
                if (seqxA_tmp[j] != '-') i1++;
                if (seqyA_tmp[j] != '-')
                {
                    i2++;
                    if (i2 >= ylen || i1 >= xlen) j = L;
                    else if (seqxA_tmp[j] != '-') invmap_tmp[i2] = i1;
                }
            }
        }

        /* superpose */
        Kabsch_Superpose(r1, r2, xt, xa, ya, xlen, ylen, invmap_tmp,
            L_ali, t, u, mol_type);

        /* derive new alignment */
        se_main(xt, ya, seqx, seqy, TM1_tmp, TM2_tmp, TM3_tmp, TM4_tmp,
            TM5_tmp, d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
            seqM_tmp, seqxA_tmp, seqyA_tmp, rmsd0_tmp, L_ali_tmp, Liden_tmp,
            TM_ali_tmp, rmsd_ali_tmp, n_ali_tmp, n_ali8_tmp, xlen, ylen,
            sequence, Lnorm_ass, d0_scale, i_opt==3, a_opt, u_opt, d_opt,
            mol_type, outfmt_opt, invmap_tmp);

        /* accept new alignment */
        if (TM1_tmp>TM1 && TM2_tmp>TM2)
        {
            /* return values */
            for (i=0; i<3; i++)
            {
                t0[i]=t[i];
                for (j=0;j<3;j++) u0[i][j]=u[i][j];
            }
            TM1=TM1_tmp;
            TM2=TM2_tmp;
            TM3=TM3_tmp;
            TM4=TM4_tmp;
            TM5=TM5_tmp;

            TM_0=TM1;
            if (a_opt>0) TM_0=TM3;
            else if (u_opt) TM_0=TM4;
            else if (d_opt) TM_0=TM5;

            if (outfmt_opt<2)
            {
                seqxA =seqxA_tmp;
                seqM  =seqM_tmp;
                seqyA =seqyA_tmp;
            }
            else for (j=0; j<ylen; j++) invmap[j]=invmap_tmp[j];

            rmsd0 =rmsd0_tmp;
            Liden =Liden_tmp;
            n_ali =n_ali_tmp;
            n_ali8=n_ali8_tmp;

            /* user specified initial alignment parameters */
            if (i_opt && L_ali==-1)
            {
                L_ali=L_ali_tmp;
                TM_ali=TM_ali_tmp;
                rmsd_ali=rmsd_ali_tmp;
            }
        }
        else
        {
            if (iter>=2) break;
            if (outfmt_opt<2)
            {
                seqxA_tmp  = seqxA;
                seqyA_tmp  = seqyA;
            }
            else for (j=0; j<ylen; j++) invmap_tmp[j]=invmap[j];
            rmsd0_tmp  = 0;
            Liden_tmp  = 0;
            n_ali_tmp  = 0;
            n_ali8_tmp = 0;
        }
    }

    /************/
    /* clean up */
    /************/
    if (outfmt_opt<2)
    {
        seqxA_tmp.clear();
        seqM_tmp.clear();
        seqyA_tmp.clear();
    }
    delete [] invmap_tmp;
    DeleteArray(&xt, xlen);
    DeleteArray(&r1, minlen);
    DeleteArray(&r2, minlen);
    return 0;
}
#endif
