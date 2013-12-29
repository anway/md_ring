
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846

int i, j;
FILE* fin = fopen("in.txt", "r");
FILE* fout = fopen("out.txt", "w");

// input vars
float zpos, ang_0, pot_lj_s, pot_lj_e, pot_hs_k, pot_hs_r, delt, tmax, temp, tequil, ztoler;
int iseed, nsamp;

// ring positions
float ring1[3][6], ring2[3][6];

// velocities
float vring1[3][6], vring2[3][6];
float vcm[3];

// other global vars
float v1, v2, v3, fxring1[6], fxring2[6], fyring1[6], fyring2[6], fzring1[6], fzring2[6];

/*
 Generate a random number between 0 and 1
 */
float randfrac()
{
    return (double)rand() / (double)RAND_MAX;
}


/*
 Rescale velocities to right temp
 */
temprescale()
{
    int sumv2=0, scale;
    for (i=0; i<=6; i++)
    {
        for (j=0; j<3; j++)
        {
            sumv2 += pow(vring1[j][i], 2.0) + pow(vring2[j][i], 2.0);
        }
    }
    sumv2 /= 14.0;
    scale = sqrt(3.0 * temp / sumv2);
    for (i=0; i<=6; i++)
    {
        for (j=0; j<3; j++)
        {
            vring1[j][i] = (vring1[j][i] - vcm[j]) * scale;
            vring2[j][i] = (vring2[j][i] - vcm[j]) * scale;
        }
    }
}

/*
 Read in values and initialize positions and velocities
 */
void init()
{
    // input
    fscanf(fin, "%f %f", &zpos, &ang_0);
    fscanf(fin, "%f %f", &pot_lj_s, &pot_lj_e);
    fscanf(fin, "%f %f", &pot_hs_k, &pot_hs_r);
    fscanf(fin, "%d %d", &iseed, &nsamp);
    fscanf(fin, "%f %f %f %f", &delt, &tmax, &temp, &tequil, &ztoler);
    
    // more initializing
    srand(iseed);
    
    float alat = 1.0;
    ang_0 = ang_0*PI/180.0;
    float ang_rot = PI/3.0;
    
    // initialize positions
    ring1[0][0] = 0.0;
    ring1[1][0] = 0.0;
    ring1[2][0] = 0.0;
    ring2[0][0] = 0.0;
    ring2[1][0] = 0.0;
    ring2[2][0] = zpos;
    
    for (i=1; i<=6; i++)
    {
        ring1[0][i] = alat * cos(ang_rot * i);
        ring1[1][i] = alat * sin(ang_rot * i);
        ring1[2][i] = 0.0;
        ring2[0][i] = alat * cos(ang_0 + ang_rot * i);
        ring2[1][i] = alat * sin(ang_0 + ang_rot * i);
        ring2[2][i] = zpos;
    }
    
    // initialize velocities
    for (i=0; i<=6; i++)
    {
        for (j=0; j<3; j++)
        {
            vring1[j][i] = randfrac() - 0.5;
            vring2[j][i] = randfrac() - 0.5;
        }
    }
    temprescale();
}

/*
 Calculate forces
 */
void force()
{
    int i, j, k;
    for (i=0; i<=6; i++)
    {
        fxring1[i]=0.0;
        fxring2[i]=0.0;
        fyring1[i]=0.0;
        fyring2[i]=0.0;
        fzring1[i]=0.0;
        fzring2[i]=0.0;
    }
    v1=0.0;
    v2=0.0;
    v3=0.0;
    // interactions in first ring
    float d10, d11, v1curr, v2curr, v3curr;
    for (i=1; i<=6; i++)
    {
        d10 = 0.0;
        d11 = 0.0;
        for (j=0; j<3; j++)
        {
            d10 += sqrt(pow(ring1[j][i]-ring1[j][0], 2.0));
        }
        v1 += pot_hs_k * pow(pot_hs_r - d10, 2.0);
        fxring1[i] += 2.0 * pot_hs_k * abs(d10 - pot_hs_r) * (ring1[0][0] - ring1[0][i]) / d10;
        fxring1[0] += 2.0 * pot_hs_k * abs(d10 - pot_hs_r) * (ring1[0][i] - ring1[0][0]) / d10;
        fyring1[i] += 2.0 * pot_hs_k * abs(d10 - pot_hs_r) * (ring1[1][0] - ring1[1][i]) / d10;
        fyring1[0] += 2.0 * pot_hs_k * abs(d10 - pot_hs_r) * (ring1[1][i] - ring1[1][0]) / d10;
        fzring1[i] += 2.0 * pot_hs_k * abs(d10 - pot_hs_r) * (ring1[2][0] - ring1[2][i]) / d10;
        fzring1[0] += 2.0 * pot_hs_k * abs(d10 - pot_hs_r) * (ring1[2][i] - ring1[2][0]) / d10;
        if (i==6)
        {
            d11 = sqrt(pow(ring1[j][i]-ring1[j][1], 2.0));
            v1 += pot_hs_k * pow(pot_hs_r - d11, 2.0);
            fxring1[i] += 2.0 * pot_hs_k * abs(d11 - pot_hs_r) * (ring1[0][1] - ring1[0][i]) / d11;
            fxring1[1] += 2.0 * pot_hs_k * abs(d11 - pot_hs_r) * (ring1[0][i] - ring1[0][1]) / d11;
            fyring1[i] += 2.0 * pot_hs_k * abs(d11 - pot_hs_r) * (ring1[1][1] - ring1[1][i]) / d11;
            fyring1[1] += 2.0 * pot_hs_k * abs(d11 - pot_hs_r) * (ring1[1][i] - ring1[1][1]) / d11;
            fzring1[i] += 2.0 * pot_hs_k * abs(d11 - pot_hs_r) * (ring1[2][1] - ring1[2][i]) / d11;
            fzring1[1] += 2.0 * pot_hs_k * abs(d11 - pot_hs_r) * (ring1[2][i] - ring1[2][1]) / d11;
        }
        else
        {
            d11 = sqrt(pow(ring1[j][i]-ring1[j][i+1], 2.0));
            v1 += pot_hs_k * pow(pot_hs_r - d11, 2.0);
            fxring1[i] += 2.0 * pot_hs_k * abs(d11 - pot_hs_r) * (ring1[0][i+1] - ring1[0][i]) / d11;
            fxring1[i+1] += 2.0 * pot_hs_k * abs(d11 - pot_hs_r) * (ring1[0][i] - ring1[0][i+1]) / d11;
            fyring1[i] += 2.0 * pot_hs_k * abs(d11 - pot_hs_r) * (ring1[1][i+1] - ring1[1][i]) / d11;
            fyring1[i+1] += 2.0 * pot_hs_k * abs(d11 - pot_hs_r) * (ring1[1][i] - ring1[1][i+1]) / d11;
            fzring1[i] += 2.0 * pot_hs_k * abs(d11 - pot_hs_r) * (ring1[2][i+1] - ring1[2][i]) / d11;
            fzring1[i+1] += 2.0 * pot_hs_k * abs(d11 - pot_hs_r) * (ring1[2][i] - ring1[2][i+1]) / d11;
        }
    }
    // interactions in second ring
    float d20, d21;
    for (i=1; i<=6; i++)
    {
        d20 = 0.0;
        d21 = 0.0;
        for (j=0; j<3; j++)
        {
            d20 += sqrt(pow(ring2[j][i]-ring2[j][0], 2.0));
        }
        v2 += pot_hs_k * pow(pot_hs_r - d20, 2.0);
        fxring2[i] += 2.0 * pot_hs_k * abs(d20 - pot_hs_r) * (ring2[0][0] - ring2[0][i]) / d20;
        fxring2[0] += 2.0 * pot_hs_k * abs(d20 - pot_hs_r) * (ring2[0][i] - ring2[0][0]) / d20;
        fyring2[i] += 2.0 * pot_hs_k * abs(d20 - pot_hs_r) * (ring2[1][0] - ring2[1][i]) / d20;
        fyring2[0] += 2.0 * pot_hs_k * abs(d20 - pot_hs_r) * (ring2[1][i] - ring2[1][0]) / d20;
        fzring2[i] += 2.0 * pot_hs_k * abs(d20 - pot_hs_r) * (ring2[2][0] - ring2[2][i]) / d20;
        fzring2[0] += 2.0 * pot_hs_k * abs(d20 - pot_hs_r) * (ring2[2][i] - ring2[2][0]) / d20;
        if (i==6)
        {
            d21 = sqrt(pow(ring2[j][i]-ring2[j][1], 2.0));
            v2 += pot_hs_k * pow(pot_hs_r - d21, 2.0);
            fxring2[i] += 2.0 * pot_hs_k * abs(d21 - pot_hs_r) * (ring2[0][1] - ring2[0][i]) / d21;
            fxring2[1] += 2.0 * pot_hs_k * abs(d21 - pot_hs_r) * (ring2[0][i] - ring2[0][1]) / d21;
            fyring2[i] += 2.0 * pot_hs_k * abs(d21 - pot_hs_r) * (ring2[1][1] - ring2[1][i]) / d21;
            fyring2[1] += 2.0 * pot_hs_k * abs(d21 - pot_hs_r) * (ring2[1][i] - ring2[1][1]) / d21;
            fzring2[i] += 2.0 * pot_hs_k * abs(d21 - pot_hs_r) * (ring2[2][1] - ring2[2][i]) / d21;
            fzring2[1] += 2.0 * pot_hs_k * abs(d21 - pot_hs_r) * (ring2[2][i] - ring2[2][1]) / d21;
        }
        else
        {
            d21 = sqrt(pow(ring2[j][i]-ring2[j][i+1], 2.0));
            v2 += pot_hs_k * pow(pot_hs_r - d21, 2.0);
            fxring2[i] += 2.0 * pot_hs_k * abs(d21 - pot_hs_r) * (ring2[0][i+1] - ring2[0][i]) / d21;
            fxring2[i+1] += 2.0 * pot_hs_k * abs(d21 - pot_hs_r) * (ring2[0][i] - ring2[0][i+1]) / d21;
            fyring2[i] += 2.0 * pot_hs_k * abs(d21 - pot_hs_r) * (ring2[1][i+1] - ring2[1][i]) / d21;
            fyring2[i+1] += 2.0 * pot_hs_k * abs(d21 - pot_hs_r) * (ring2[1][i] - ring2[1][i+1]) / d21;
            fzring2[i] += 2.0 * pot_hs_k * abs(d21 - pot_hs_r) * (ring2[2][i+1] - ring2[2][i]) / d21;
            fzring2[i+1] += 2.0 * pot_hs_k * abs(d21 - pot_hs_r) * (ring2[2][i] - ring2[2][i+1]) / d21;
        }
    }
    // interactions between rings
    float d3, r3_6, r3_12;
    for (i=0; i<=6; i++)
    {
        for (j=0; j<=6; j++)
        {
            d3 = 0.0;
            for (k=0; k<3; k++)
            {
                d3 += sqrt(pow(ring1[k][i]-ring2[k][j], 2.0));
            }
            r3_6 = pow(pot_lj_s / d3, 6.0);
            r3_12 = r3_6 * r3_6;
            v3 += pot_lj_e * (r3_12 - r3_6);
            r3_8 = pow(pot_lj_s / d3, 8.0);
            fxring1[i] += 12 * pot_lj_e * (r3_6 * r3_8 - 0.5 * r3_8) * (ring1[0][i] - ring2[0][j]);
            fxring2[j] += 12 * pot_lj_e * (r3_6 * r3_8 - 0.5 * r3_8) * (ring2[0][j] - ring1[0][i]);
            fyring1[i] += 12 * pot_lj_e * (r3_6 * r3_8 - 0.5 * r3_8) * (ring1[1][i] - ring2[1][j]);
            fyring2[j] += 12 * pot_lj_e * (r3_6 * r3_8 - 0.5 * r3_8) * (ring2[1][j] - ring1[1][i]);
            fzring1[i] += 12 * pot_lj_e * (r3_6 * r3_8 - 0.5 * r3_8) * (ring1[2][i] - ring2[2][j]);
            fzring2[j] += 12 * pot_lj_e * (r3_6 * r3_8 - 0.5 * r3_8) * (ring2[2][j] - ring1[2][i]);
        }
    }

}

/*
 Find next positions
 */
void solve()
{
    int i, j;
    for (i=0; i<=6; i++)
    {
        vring1[0][i] += delt * fxring1[i];
        vring1[1][i] += delt * fyring1[i];
        vring1[2][i] += delt * fzring1[i];
        for (j=0; j<3; j++)
        {
            ring1[j][i] += delt * vring1[j][i];
        }
        vring2[0][i] += delt * fxring2[i];
        vring2[1][i] += delt * fyring2[i];
        vring2[2][i] += delt * fzring2[i];
        for (j=0; j<3; j++)
        {
            ring2[j][i] += delt * vring2[j][i];
        }
    }
}

/*
 Sample parameters of interest
 */
void sample()
{

}

int main(void)
{
    init();
    
    int step=0, time=0, nstep = (int) tmax / delt;
               
    do {
        force();
        solve();
        time += delt;
        if (time < tequil)
        {
            if (time / delt % 20 == 0)
                temprescale();
        }
        else if ( time / delt % 14 == 0)
            sample();
    } while (time < tmax);
    
    // write things of interest to file
    
    return 0;
}
