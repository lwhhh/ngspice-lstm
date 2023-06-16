/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"

#define ccap (qcap+1)


//add dyc gai we could use the one plus order to calculate instead using the equal order to use beacuse of the extra pre x and dt
double diff_calculate(int qcap, CKTcircuit* ckt)
{
    double diff[9];
    double deltmp[9];
    int i, j;
    diff[0] = ckt->QState[ckt->Qindex-1];
    for (i = ckt->CKTorder + 1; i >= 0; i--) {
        diff[i+1] = ckt->CKTstates[i][qcap];
    }
    deltmp[0] = ckt->pre_timestep_Q[ckt->Qindex - 1];
    for (i = 0; i <= ckt->CKTorder; i++) {
        deltmp[i+1] = ckt->CKTdeltaOld[i];
    }
    j = ckt->CKTorder + 1;
    for (;;) {
        for (i = 0; i <= j; i++) {
            diff[i] = (diff[i] - diff[i + 1]) / deltmp[i];
        }
        if (--j < 0) break;
        for (i = 0; i <= j; i++) {
            deltmp[i] = deltmp[i + 1] + ckt->CKTdeltaOld[i];
        }
    }
    return diff[1];

}

void
CKTterr(int qcap, CKTcircuit *ckt, double *timeStep)
{ 
    double volttol;
    double chargetol;
    double tol;
    double del;
    double diff[8];
    double deltmp[8];
    double factor=0;
    int i;
    int j;
    int  shit = 0;
    static double gearCoeff[] = {
        .5,
        .2222222222,
        .1363636364,
        .096,
        .07299270073,
        .05830903790
    };
    static double trapCoeff[] = {
        .5,
        .08333333333
    };
    ckt->QState[ckt->Qindex] = ckt->CKTstate0[qcap];
    ckt->Qcap[ckt->Qindex] = qcap;
    shit = ckt->Qindex;
    ckt->Qindex = shit + 1;// add dyc store the Qcap and state
    if (ckt->Qindex >= 2000)
    {
        int fuck = 1;
    }
    volttol = ckt->CKTabstol + ckt->CKTreltol * 
            MAX( fabs(ckt->CKTstate0[ccap]), fabs(ckt->CKTstate1[ccap]));
            
    chargetol = MAX(fabs(ckt->CKTstate0[qcap]),fabs(ckt->CKTstate1[qcap]));
    chargetol = ckt->CKTreltol * MAX(chargetol,ckt->CKTchgtol)/ckt->CKTdelta;
    tol = MAX(volttol,chargetol);
    /* now divided differences */
    for(i=ckt->CKTorder+1;i>=0;i--) {
        diff[i] = ckt->CKTstates[i][qcap];
    }
    for(i=0 ; i <= ckt->CKTorder ; i++) {
        deltmp[i] = ckt->CKTdeltaOld[i];
    }
    j = ckt->CKTorder;
    for (;;) {
        for(i=0;i <= j;i++) {
            diff[i] = (diff[i] - diff[i+1])/deltmp[i];
        }
        if (--j < 0) break;
        for(i=0;i <= j;i++) {
            deltmp[i] = deltmp[i+1] + ckt->CKTdeltaOld[i];
        }
    }
    if (ckt->use_Qdiff_flag == 1 && ckt->pre_timestep_Q[ckt->Qindex - 1]>1e-12)
    {
        int dada = 0;
        diff[0] = diff_calculate(qcap, ckt); //add dyc use the changed diff
    }
    switch(ckt->CKTintegrateMethod) {
        case GEAR:
            factor = gearCoeff[ckt->CKTorder-1];
            break;

        case TRAPEZOIDAL:
            factor = trapCoeff[ckt->CKTorder - 1] ;
            break;
    }
    del = ckt->CKTtrtol * tol/MAX(ckt->CKTabstol,factor * fabs(diff[0]));
    if(ckt->CKTorder == 2) {
        del = sqrt(del);
    } else if (ckt->CKTorder > 2) {
        del = exp(log(del)/ckt->CKTorder);
    }
    ckt->deltatimeinQ[ckt->Qindex - 1] = del;
    *timeStep = MIN(*timeStep,del);
    return;
}
