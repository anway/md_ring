
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846

// input vars
float ang_0, pot_lj_s, pot_lj_e, pot_hs_k, pot_hs_r, k2, z0, delt, tmax, temp, tequil;
int iseed;

// ring positions
float ring1[3][7], ring2[3][7];

// velocities
float vring1[3][7], vring2[3][7];
float vcm[3];

// other global vars
float v1=0.0, v2=0.0, v3=0.0, v12=0.0, fxring1[7], fxring2[7], fyring1[7], fyring2[7], fzring1[7], fzring2[7];

float vtot=0.0;

/*
 Generate a random number between 0 and 1
 */
float randfrac()
{
    return (double)rand() / (double)RAND_MAX;
}

/*
 Generate random number from Gaussian distribution
 */
float randgauss()
{
    float x1, x2, w, y1, y2;
    
    do {
        x1 = 2.0 * randfrac() - 1.0;
        x2 = 2.0 * randfrac() - 1.0;
        w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    return x1 * sqrt( (-2.0 * log( w ) ) / w );
}


/*
 Rescale velocities to right temp
 */
void temprescale()
{
    int i, j;
    float sumv2=0, scale;
    for (i=0; i<3; i++)
        vcm[i] = 0;
    for (i=0; i<=6; i++)
    {
        for (j=0; j<3; j++)
        {
            sumv2 += pow(vring1[j][i], 2.0) + pow(vring2[j][i], 2.0);
            vcm[j] += vring1[j][i] + vring2[j][i];
        }
    }
    sumv2 /= 14.0;
    scale = sqrt(3.0 * temp / sumv2);
    for (i=0; i<3; i++)
        vcm[i] /= 14;
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
void init(FILE* fin)
{
    // input
    fscanf(fin, "%f", &ang_0);
    fscanf(fin, "%f %f", &pot_lj_s, &pot_lj_e);
    fscanf(fin, "%f %f %f %f", &pot_hs_k, &pot_hs_r, &k2, &z0);
    fscanf(fin, "%d", &iseed);
    fscanf(fin, "%f %f %f %f", &delt, &tmax, &temp, &tequil);
    
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
    ring2[2][0] = z0;
    
    int i, j;
    for (i=1; i<=6; i++)
    {
        ring1[0][i] = alat * cos(ang_rot * (i-1));
        ring1[1][i] = alat * sin(ang_rot * (i-1));
        ring1[2][i] = 0.0;
        ring2[0][i] = alat * cos(ang_0 + ang_rot * (i-1));
        ring2[1][i] = alat * sin(ang_0 + ang_rot * (i-1));
        ring2[2][i] = z0;
    }
    
    // initialize velocities
    for (i=0; i<=6; i++)
    {
        for (j=0; j<3; j++)
        {
            vring1[j][i] = randgauss();
            vring2[j][i] = randgauss();
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
            d10 += pow(ring1[j][i]-ring1[j][0], 2.0);
        }
        d10 = sqrt(d10);
        v1 += pot_hs_k * pow(pot_hs_r - d10, 2.0);
        fxring1[i] += 2.0 * pot_hs_k * (d10 - pot_hs_r) * (ring1[0][0] - ring1[0][i]) / d10;
        fxring1[0] += 2.0 * pot_hs_k * (d10 - pot_hs_r) * (ring1[0][i] - ring1[0][0]) / d10;
        fyring1[i] += 2.0 * pot_hs_k * (d10 - pot_hs_r) * (ring1[1][0] - ring1[1][i]) / d10;
        fyring1[0] += 2.0 * pot_hs_k * (d10 - pot_hs_r) * (ring1[1][i] - ring1[1][0]) / d10;
        fzring1[i] += 2.0 * pot_hs_k * (d10 - pot_hs_r) * (ring1[2][0] - ring1[2][i]) / d10;
        fzring1[0] += 2.0 * pot_hs_k * (d10 - pot_hs_r) * (ring1[2][i] - ring1[2][0]) / d10;
        
        if (i==6)
        {
            for (j=0; j<3; j++)
            {
                d11 += pow(ring1[j][i]-ring1[j][1], 2.0);
            }
            d11 = sqrt(d11);
            v1 += pot_hs_k * pow(pot_hs_r - d11, 2.0);
            fxring1[i] += 2.0 * pot_hs_k * (d11 - pot_hs_r) * (ring1[0][1] - ring1[0][i]) / d11;
            fxring1[1] += 2.0 * pot_hs_k * (d11 - pot_hs_r) * (ring1[0][i] - ring1[0][1]) / d11;
            fyring1[i] += 2.0 * pot_hs_k * (d11 - pot_hs_r) * (ring1[1][1] - ring1[1][i]) / d11;
            fyring1[1] += 2.0 * pot_hs_k * (d11 - pot_hs_r) * (ring1[1][i] - ring1[1][1]) / d11;
            fzring1[i] += 2.0 * pot_hs_k * (d11 - pot_hs_r) * (ring1[2][1] - ring1[2][i]) / d11;
            fzring1[1] += 2.0 * pot_hs_k * (d11 - pot_hs_r) * (ring1[2][i] - ring1[2][1]) / d11;
        }
        else
        {
            for (j=0; j<3; j++)
            {
                d11 += pow(ring1[j][i]-ring1[j][i+1], 2.0);
            }
            d11 = sqrt(d11);
            v1 += pot_hs_k * pow(pot_hs_r - d11, 2.0);
            fxring1[i] += 2.0 * pot_hs_k * (d11 - pot_hs_r) * (ring1[0][i+1] - ring1[0][i]) / d11;
            fxring1[i+1] += 2.0 * pot_hs_k * (d11 - pot_hs_r) * (ring1[0][i] - ring1[0][i+1]) / d11;
            fyring1[i] += 2.0 * pot_hs_k * (d11 - pot_hs_r) * (ring1[1][i+1] - ring1[1][i]) / d11;
            fyring1[i+1] += 2.0 * pot_hs_k * (d11 - pot_hs_r) * (ring1[1][i] - ring1[1][i+1]) / d11;
            fzring1[i] += 2.0 * pot_hs_k * (d11 - pot_hs_r) * (ring1[2][i+1] - ring1[2][i]) / d11;
            fzring1[i+1] += 2.0 * pot_hs_k * (d11 - pot_hs_r) * (ring1[2][i] - ring1[2][i+1]) / d11;
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
            d20 += pow(ring2[j][i]-ring2[j][0], 2.0);
        }
        d20 = sqrt(d20);
        v2 += pot_hs_k * pow(pot_hs_r - d20, 2.0);
        fxring2[i] += 2.0 * pot_hs_k * (d20 - pot_hs_r) * (ring2[0][0] - ring2[0][i]) / d20;
        fxring2[0] += 2.0 * pot_hs_k * (d20 - pot_hs_r) * (ring2[0][i] - ring2[0][0]) / d20;
        fyring2[i] += 2.0 * pot_hs_k * (d20 - pot_hs_r) * (ring2[1][0] - ring2[1][i]) / d20;
        fyring2[0] += 2.0 * pot_hs_k * (d20 - pot_hs_r) * (ring2[1][i] - ring2[1][0]) / d20;
        fzring2[i] += 2.0 * pot_hs_k * (d20 - pot_hs_r) * (ring2[2][0] - ring2[2][i]) / d20;
        fzring2[0] += 2.0 * pot_hs_k * (d20 - pot_hs_r) * (ring2[2][i] - ring2[2][0]) / d20;
        if (i==6)
        {
            for (j=0; j<3; j++)
            {
                d21 += pow(ring2[j][i]-ring2[j][1], 2.0);
            }
            d21 = sqrt(d21);
            v2 += pot_hs_k * pow(pot_hs_r - d21, 2.0);
            fxring2[i] += 2.0 * pot_hs_k * (d21 - pot_hs_r) * (ring2[0][1] - ring2[0][i]) / d21;
            fxring2[1] += 2.0 * pot_hs_k * (d21 - pot_hs_r) * (ring2[0][i] - ring2[0][1]) / d21;
            fyring2[i] += 2.0 * pot_hs_k * (d21 - pot_hs_r) * (ring2[1][1] - ring2[1][i]) / d21;
            fyring2[1] += 2.0 * pot_hs_k * (d21 - pot_hs_r) * (ring2[1][i] - ring2[1][1]) / d21;
            fzring2[i] += 2.0 * pot_hs_k * (d21 - pot_hs_r) * (ring2[2][1] - ring2[2][i]) / d21;
            fzring2[1] += 2.0 * pot_hs_k * (d21 - pot_hs_r) * (ring2[2][i] - ring2[2][1]) / d21;
        }
        else
        {
            for (j=0; j<3; j++)
            {
                d21 += pow(ring2[j][i]-ring2[j][i+1], 2.0);
            }
            d21 = sqrt(d21);
            v2 += pot_hs_k * pow(pot_hs_r - d21, 2.0);
            fxring2[i] += 2.0 * pot_hs_k * (d21 - pot_hs_r) * (ring2[0][i+1] - ring2[0][i]) / d21;
            fxring2[i+1] += 2.0 * pot_hs_k * (d21 - pot_hs_r) * (ring2[0][i] - ring2[0][i+1]) / d21;
            fyring2[i] += 2.0 * pot_hs_k * (d21 - pot_hs_r) * (ring2[1][i+1] - ring2[1][i]) / d21;
            fyring2[i+1] += 2.0 * pot_hs_k * (d21 - pot_hs_r) * (ring2[1][i] - ring2[1][i+1]) / d21;
            fzring2[i] += 2.0 * pot_hs_k * (d21 - pot_hs_r) * (ring2[2][i+1] - ring2[2][i]) / d21;
            fzring2[i+1] += 2.0 * pot_hs_k * (d21 - pot_hs_r) * (ring2[2][i] - ring2[2][i+1]) / d21;
        }
    }
    
    // central spring
    float d12=0.0;
    for (j=0; j<3; j++)
    {
        d12 += pow(ring1[j][0]-ring2[j][0], 2.0);
    }
    d12 = sqrt(d12);
    v12 = k2 * pow(z0 - d12, 2.0);
    fxring1[0] += 2.0 * k2 * (d12 - z0) * (ring2[0][0] - ring1[0][0]) / d12;
    fxring2[0] += 2.0 * k2 * (d12 - z0) * (ring1[0][0] - ring2[0][0]) / d12;
    fyring1[0] += 2.0 * k2 * (d12 - z0) * (ring2[1][0] - ring1[1][0]) / d12;
    fyring2[0] += 2.0 * k2 * (d12 - z0) * (ring1[1][0] - ring2[1][0]) / d12;
    fzring1[0] += 2.0 * k2 * (d12 - z0) * (ring2[2][0] - ring1[2][0]) / d12;
    fzring2[0] += 2.0 * k2 * (d12 - z0) * (ring1[2][0] - ring2[2][0]) / d12;
    
    // interactions between rings
    float d3, r3_6, r3_12, r3_8;
    for (i=0; i<=6; i++)
    {
        for (j=0; j<=6; j++)
        {
            d3 = 0.0;
            for (k=0; k<3; k++)
            {
                d3 += pow(ring1[k][i]-ring2[k][j], 2.0);
            }
            d3 = sqrt(d3);
            r3_6 = pow(pot_lj_s / d3, 6.0);
            r3_12 = r3_6 * r3_6;
            v3 += pot_lj_e * (r3_12 - r3_6);
            //fprintf(fout, "%d, %d: %f %f\n", i, j, r3_12, r3_6);
            fxring1[i] += 12.0 * pot_lj_e * (r3_12 / d3 - 0.5 * r3_6 / d3) * (ring1[0][i] - ring2[0][j]) / d3;
            fxring2[j] += 12.0 * pot_lj_e * (r3_12 / d3 - 0.5 * r3_6 / d3) * (ring2[0][j] - ring1[0][i]) / d3;
            fyring1[i] += 12.0 * pot_lj_e * (r3_12 / d3 - 0.5 * r3_6 / d3) * (ring1[1][i] - ring2[1][j]) / d3;
            fyring2[j] += 12.0 * pot_lj_e * (r3_12 / d3 - 0.5 * r3_6 / d3) * (ring2[1][j] - ring1[1][i]) / d3;
            fzring1[i] += 12.0 * pot_lj_e * (r3_12 / d3 - 0.5 * r3_6 / d3) * (ring1[2][i] - ring2[2][j]) / d3;
            fzring2[j] += 12.0 * pot_lj_e * (r3_12 / d3 - 0.5 * r3_6 / d3) * (ring2[2][j] - ring1[2][i]) / d3;
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
    vtot += v1 + v2 + v3 + v12;
}

/*
 Write energy
 */
void writeener(FILE* fener)
{
    float sumv2 = 0.0;
    int i;
    for (i=0; i<=6; i++)
    {
        sumv2 += pow(vring1[0][i], 2.0) + pow(vring1[1][i], 2.0) + pow(vring1[2][i], 2.0) + pow(vring2[0][i], 2.0) + pow(vring2[1][i], 2.0) + pow(vring2[2][i], 2.0);
    }
    sumv2 /= 2.0;
    fprintf(fener, "ring1: %f, ring2: %f, lj: %f, kinetic: %f, sum: %f\n", v1, v2, v3, sumv2, v1+v2+v3+sumv2);
}

int main(void)
{
    
    int i, j;
    FILE* fin = fopen("in.txt", "r");
    FILE* fout = fopen("out.txt", "w");
    FILE* fener = fopen("ener.txt", "w");
    FILE* ftraj = fopen("md_ring.xyz", "w");
    FILE* fp = fopen("p.txt", "w");
    
    init(fin);
    
    int step=0, nstep = (int) tmax / delt;
    float time;
    
    do {
        force();
        solve();
        
        //write(ftraj);
        //write_py();
        writeener(fener);
        //writep(fp);
        
        time += delt;
        if (time < tequil)
        {
            if ((int)(time / delt) % 20 == 0)
                temprescale();
        }
        //else if ((int)(time / delt) % 14 == 0)
            sample();
    } while (time < tmax);
    
    float epsilon0 = -10.0 * temp;
    vtot /= (int)(tmax/delt);
    vtot += epsilon0;
    
    FILE* data = fopen("plot3d.py", "a");
    fprintf(data, "[%f, %f, %f]", temp, z0, vtot);

    return 0;
}
