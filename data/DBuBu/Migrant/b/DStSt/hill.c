//
//
/*
 individual-based simulations of hill-climbing
 island model
 each locus is a bit
 in rec gens, loci are only unlinked, mixing is not perfect
 */
/**********************************************/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <gsl/gsl_randist.h>
#include <time.h>

typedef struct
{
    unsigned long popsize, numdemes; //size of one deme, number of demes
    unsigned int numloci, tgap;
    double frac, mutrate, s, recrate; //note that mutrate is genomic mutation rate
    char *name;                       //the prefix for output files
} parameters;

char **pop_setup(parameters *params, gsl_rng *RNDM);

double *landscape(parameters *params, gsl_rng *RNDM);

int generation(parameters *params, char **population, double **W, unsigned int **mutcount, double *sels, unsigned long t, gsl_rng *RNDM);

int add_data(unsigned int **arr1, int value1, unsigned int **arr2, int value2, unsigned int **cap_size, unsigned int parity, unsigned int init_cap);

int main(int argc, char *argv[])
{
    unsigned long SEED, N;
    unsigned long t, tfinal, i, j, k, D;
    unsigned int L, d, parity, Tgap, **mutcount, run, fittgap = 10, **max_fit_ind, **max_fit_deme, **cap_size, dist, relateness;
    unsigned int ind1, d1, ind2, d2, init_cap;
    char **population;
    double *sels, **W, sbar, *max_fit, *freq_max_fit, fit_temp, **gene_freq, het;
    clock_t start, end;
    //    double count;
    parameters *params = malloc(sizeof(parameters));
    params->name = malloc(100 * sizeof(char));
    char *outfile, *mutn, *recn;
    FILE *paramfile, *fitfile, *tfitfile, *hetfile, *genfile, *relfile; //, *genfile, *divfile

    start = clock();
    //Initialize variables:
    //argc should be one more than the number of required parameters:
    if (argc != 10)
    {
        printf("Usage: s tgap U V D N L tfinal ranseed outfile\n");
        return EXIT_FAILURE;
    }
    j = 1;
    params->s = atof(argv[j++]); //characteristic selection coefficient (log fitness)
    Tgap = atof(argv[j++]);
    params->tgap = Tgap;
    mutn = argv[j];
    params->mutrate = atof(argv[j++]); //genomic mutation rate
    recn = argv[j];
    params->recrate = atof(argv[j++]);   //genomic recombination rate
    D = (unsigned long)atoll(argv[j++]); //number of demes
    params->numdemes = D;
    N = (unsigned long)atoll(argv[j++]); //population size of each deme
    params->popsize = N;
    L = (unsigned int)atol(argv[j++]); //number of "loci" (each of which is a char)
    params->numloci = L;
    tfinal = (unsigned long)atoll(argv[j++]); //stopping time
    SEED = (unsigned long)atol(argv[j]);      //RNG seed
    outfile = argv[j];                        //name for output files
    strcpy(params->name, outfile);

    //gsl random setup:
    const gsl_rng_type *T;
    gsl_rng *RNDM;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    RNDM = gsl_rng_alloc(T);

    for (; params->tgap <= 300; params->tgap += 1000)
    {
        Tgap = params->tgap;
        params->frac = 0.0025;
        //Print out variables to paramfile:
        sprintf(outfile, "%d%c%ld%c%ld%c%s%c%s", params->tgap, 'D', params->numdemes, 'N', params->popsize, 'M', mutn, 'R', recn);
        strcpy(params->name, outfile);
        fitfile = fopen(strcat(outfile, ".fit"), "w");
        strcpy(outfile, params->name);
        paramfile = fopen(strcat(outfile, ".params"), "w");
        strcpy(outfile, params->name);
        tfitfile = fopen(strcat(outfile, ".tfit"), "w");
        strcpy(outfile, params->name);
        hetfile = fopen(strcat(outfile, ".het"), "w");
        strcpy(outfile, params->name);
        genfile = fopen(strcat(outfile, ".gen"), "w");
        strcpy(outfile, params->name);
        relfile = fopen(strcat(outfile, ".rel"), "w");

        fprintf(paramfile, "ranseed = %ld,frac = %lf,tgap = %u,s = %lf,U = %lf,V = %lf,D = %lu,N = %lu,L = %u,tfinal = %lu\n", SEED, params->frac, params->tgap, params->s, params->mutrate, params->recrate, params->numdemes, N, L, tfinal);
        //initialize random number generator
        gsl_rng_set(RNDM, SEED);
        for (run = 0; run < 5; run++)
        {
            sbar = 0;
            mutcount = malloc((D) * sizeof(unsigned int *));
            for (d = 0; d < D; d++)
            {
                mutcount[d] = calloc(2 * N, sizeof(unsigned int));
            }
            gene_freq = malloc(L * sizeof(double *));
            for (i = 0; i < L; i++)
            {
                gene_freq[i] = calloc(D + 1, sizeof(double));
            }
            cap_size = malloc(2 * sizeof(unsigned int *));     // capaticy and current size
            max_fit_ind = malloc(2 * sizeof(unsigned int *));  // maxfit individual
            max_fit_deme = malloc(2 * sizeof(unsigned int *)); // demes of maxfit individual
            init_cap = (unsigned int)(0.1 * N);
            for (i = 0; i < 2; i++) //2 means old and new generations
            {
                cap_size[i] = calloc(2, sizeof(unsigned int));
                cap_size[i][0] = init_cap;                                      //capacity
                cap_size[i][1] = 1;                                             //initially, the first cell is the best-fit(they are the same)
                max_fit_ind[i] = calloc(cap_size[i][0], sizeof(unsigned int));  // position of the individual
                max_fit_deme[i] = calloc(cap_size[i][0], sizeof(unsigned int)); // deme of that individual
            }

            freq_max_fit = calloc(D + 1, sizeof(double));
            max_fit = calloc(D + 1, sizeof(double));
            //initialize population as numdeme demes of N individuals, each with L "loci"
            //population is an unsigned char **
            population = pop_setup(params, RNDM);

            //initialize selective coefficients at each locus:
            sels = landscape(params, RNDM);

            //initialize individual fitnesses:
            W = malloc(D * sizeof(double *));
            for (d = 0; d < D; d++)
            {
                W[d] = calloc(N, sizeof(double));
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < L; j++)
                    {
                        mutcount[d][i] += population[d][i * L + j];
                    }
                    W[d][i] = exp(sels[mutcount[d][i]]);
                }
            }

            //Evolution:
            for (t = 0; t < tfinal; t++)
            {
                parity = t % 2;
                //evolve one generation:
                if (generation(params, population, W, mutcount, sels, t, RNDM) != EXIT_SUCCESS)
                {
                    printf("Problem in generation %lu\n", t);
                    return EXIT_FAILURE;
                }

                if (t % fittgap == fittgap - 1 || (t + 1) % fittgap == fittgap - 1 || t % Tgap == Tgap - 1 || (t + 1) % Tgap == Tgap - 1 || t % Tgap == 0) // at the output gen or right brfore the output gen, find max-fit individule
                {
                    sbar = 0; // calculate mean fitness and max fitness
                    max_fit[D] = 0;
                    for (d = 0; d < D; d++)
                    {
                        max_fit[d] = 0;
                        for (i = 0; i < N; i++)
                        {
                            fit_temp = log(W[d][i]);
                            sbar += fit_temp / (N * (D)); //mean fitness
                            if (fit_temp > max_fit[d])
                            {
                                max_fit[d] = fit_temp;
                            }
                            if (fit_temp > max_fit[D])
                            {
                                max_fit[D] = fit_temp;
                            }
                        }
                    }

                    //record the individual with highest genotype
                    cap_size[1 - parity][1] = 0; //No genotype recorded
                    for (d = 0; d < D; d++)      //find the one with highest fitness
                    {
                        for (i = 0; i < N; i++)
                        {
                            fit_temp = log(W[d][i]);
                            if (fit_temp == max_fit[D])
                            {
                                add_data(max_fit_deme, d, max_fit_ind, i + (1 - parity) * N, cap_size, parity, init_cap);
                            }
                        }
                    }
                }

                if (t % fittgap == fittgap - 1 || t % Tgap == Tgap - 1 || t % Tgap == 0)
                {

                    fprintf(tfitfile, "%lu\t %lf\t", t, sbar);

                    freq_max_fit[D] = 0; //output frequencies of most fit ones
                    fprintf(fitfile, "%lu\t", t);
                    for (d = 0; d < D; d++)
                    {
                        freq_max_fit[d] = 0;
                        for (i = 0; i < N; i++)
                        {
                            fit_temp = log(W[d][i]);
                            if (fit_temp >= max_fit[d] - params->s)//allow a difference of 1 mutation 
                            {
                                freq_max_fit[d] += 1.0 / N;
                            }
                            if (fit_temp >= max_fit[D] - params->s)//allow a difference of 1 mutation 
                            {
                                freq_max_fit[D] += 1.0 / (N * D);
                            }
                        }
                        fprintf(fitfile, "%f\t%f\t,\t", max_fit[d], freq_max_fit[d]);
                    }
                    fprintf(fitfile, "%f\t%f\t;\t", max_fit[D], freq_max_fit[D]);

                    for (i = 0; i < L; i++) //update gene frequency
                    {
                        gene_freq[i][D] = 0; //total gene frequency
                        for (d = 0; d < D; d++)
                        {
                            gene_freq[i][d] = 0; // gene frequency in a deme
                            for (j = (1 - parity) * N; j < ((2 - parity) * N); j++)
                            {
                                gene_freq[i][d] += 1.0 * population[d][j * L + i] / N;
                            }
                            gene_freq[i][D] += 1.0 * gene_freq[i][d] / D;
                        }
                    }
                    fprintf(hetfile, "%lu\t", t); //output heterozygosity
                    for (d = 0; d <= D; d++)      //last one is population-wide heterozygosity
                    {
                        het = 0;
                        for (i = 0; i < L; i++)
                        {
                            het += 2.0 * gene_freq[i][d] * (1.0 - gene_freq[i][d]);
                        }
                        fprintf(hetfile, "%lf\t", het);
                    }
                    fprintf(hetfile, ",\t");

                    relateness = 0; //calculate the smallest relateness between best-fit individuals of two generation
                    for (i = 0; i < cap_size[1 - parity][1]; i++)
                    {
                        ind1 = max_fit_ind[1 - parity][i];
                        d1 = max_fit_deme[1 - parity][i];
                        for (j = 0; j < cap_size[parity][1]; j++)
                        {
                            dist = 0;
                            ind2 = max_fit_ind[parity][j];
                            d2 = max_fit_deme[parity][j];
                            for (k = 0; k < L; k++) //calculate the distance between two individual;
                            {
                                if (population[d1][ind1 * L + k] == population[d2][ind2 * L + k])
                                {
                                    dist++;
                                }
                            }
                            if (dist > relateness)
                            {
                                relateness = dist;
                            }
                        }
                    }
                    fprintf(relfile, "%lu\t%u\t,\t", t, relateness);
                }

                if (t % Tgap == Tgap - 1)
                {
                    fprintf(genfile, "%lu\t,\t", t);
                    for ( i = 0; i < L; i++)
                    {
                        fprintf(genfile, "%lf\t", gene_freq[i][D]);
                    }
                    fprintf(genfile, ";\t");
                }
                
            }

            for (i = 0; i < L; i++) //output gene profile
            {
                fprintf(genfile, "%lu\t,\t", t-1);
                    for ( i = 0; i < L; i++)
                    {
                        fprintf(genfile, "%lf\t", gene_freq[i][D]);
                    }
                    fprintf(genfile, ";\t");
            }

            for (d = 0; d < D; d++)
            {
                free(mutcount[d]);
                free(population[d]);
                free(W[d]);
            }
            for (i = 0; i < L; i++)
            {
                free(gene_freq[i]);
            }
            for (i = 0; i < 2; i++)
            {
                free(max_fit_deme[i]);
                free(max_fit_ind[i]);
                free(cap_size[i]);
            }
            free(max_fit);
            free(freq_max_fit);
            free(mutcount);
            free(population);
            free(W);
            free(sels);
            free(gene_freq);
            free(cap_size);
            free(max_fit_ind);
            free(max_fit_deme);
            fprintf(tfitfile, "\n");
            fprintf(fitfile, "\n");
            fprintf(hetfile, "\n");
            fprintf(relfile, "\n");
            fprintf(genfile, "\n");
        }
        fclose(paramfile);
        fclose(fitfile);
        fclose(hetfile);
        fclose(tfitfile);
        fclose(genfile);
        fclose(relfile);
    }
    end = clock();
    printf("Time: %lf\n", (double)(end - start) / CLOCKS_PER_SEC);
    return EXIT_SUCCESS;
}

/********************************************/
//initializes population
//thanks to J Kelleher for advice
char **pop_setup(parameters *params, gsl_rng *RNDM)
{
    unsigned int d, L = (params->numloci);
    unsigned long N = (params->popsize), D = (params->numdemes);
    char **population = malloc(D * sizeof(char *));
    if (population == NULL)
    {
        printf("Cannot allocate memory\n");
        exit(EXIT_FAILURE);
    }
    for (d = 0; d < D; d++)
    {
        population[d] = calloc(2 * N * L, sizeof(char));
        if (population[d] == NULL)
        {
            printf("Cannot allocate memory\n");
            exit(EXIT_FAILURE);
        }
    }
    //need 2N spots, for parental and offspring generations
    //population is initialized with no mutations
    return population;
}

/********************************************/
//initializes fitness landscape
//fitness is always just a function of the number of 1s in your genotype
double *landscape(parameters *params, gsl_rng *RNDM)
{
    int l, L = (params->numloci);
    double *sels = calloc(L + 1, sizeof(double));

    //generate the selective coefficient for each locus
    if (sels == NULL)
    {
        printf("Cannot allocate memory\n");
        exit(EXIT_FAILURE);
    }
    for (l = 0; l <= L; l++)
    {
        //hill:
        sels[l] = l * params->s;
    }
    return sels;
}

/********************************************/
//does one generation of evolution
//selection -> (migration/recombination) ->mutation
int generation(parameters *params, char **population, double **W, unsigned int **mutcount, double *sels, unsigned long t, gsl_rng *RNDM)
{
    unsigned long i, d, dpar[2], j, parent[2], N = params->popsize, D = params->numdemes;
    unsigned int y, L = (params->numloci), parity = t % 2;

    //define per-locus mutation rate:
    double m = (params->mutrate) / (double)L;
    int vargen = 1, migrants = 0, recombine = 0;
    double rec = params->recrate;

    //initialize lookup table for population sampling:
    gsl_ran_discrete_t *gametes[D];
    for (d = 0; d < D; d++)
    {
        gametes[d] = gsl_ran_discrete_preproc(N, W[d]);
    }
    //check if it's a migrant/recombinant generation:
    // if (t % params->tgap == params->tgap - 1)
    // {
    //     vargen = 1;
    //     rec = params->recrate * params->tgap;
    // }

    //Generate new demes:
    for (d = 0; d < D; d++)
    {
        //Generate N new individuals:
        //figure out how many should be migrants/recombinants:
        if (vargen)
        {
            migrants = gsl_ran_binomial(RNDM, params->frac, N);
        }
        for (i = (1 - parity) * N; i < ((2 - parity) * N); i++)
        {
            //reset mutcount:
            mutcount[d][i] = 0;

            if (i - (1 - parity) * N < migrants)
            { //the individual is a migrant
                recombine = gsl_ran_bernoulli(RNDM, rec);
                if (recombine)
                { //if it is a recombinant
                    //draw from parents from random demes:
                    for (j = 0; j < 2; j++)
                    {
                        dpar[j] = gsl_rng_uniform_int(RNDM, D);
                        parent[j] = parity * N + gsl_ran_discrete(RNDM, gametes[dpar[j]]);
                    }
                    //loci are unlinked:
                    for (y = 0; y < L; y++)
                    {
                        if (population[dpar[0]][parent[0] * L + y] == population[dpar[1]][parent[1] * L + y])
                        { //the parents are identical at the locus
                            population[d][i * L + y] = population[dpar[0]][parent[0] * L + y];
                        }
                        else
                        { //the parents differ at the locus:
                            population[d][i * L + y] = gsl_ran_bernoulli(RNDM, 0.5);
                        }
                        mutcount[d][i] += population[d][i * L + y];
                    }
                }
                else
                { //if it is not a recombinant
                    dpar[0] = gsl_rng_uniform_int(RNDM, D);
                    parent[0] = parity * N + gsl_ran_discrete(RNDM, gametes[dpar[0]]);
                    memcpy(population[d] + i * L, population[dpar[0]] + parent[0] * L, L * sizeof(unsigned char));
                    mutcount[d][i] = mutcount[dpar[0]][parent[0]];
                }
            }
            else
            { // the individual is not a migrant
                recombine = gsl_ran_bernoulli(RNDM, rec);
                if (recombine)
                { //if it is a recombinant
                    //draw from parents from the same deme:
                    for (j = 0; j < 2; j++)
                    {
                        parent[j] = parity * N + gsl_ran_discrete(RNDM, gametes[d]);
                    }
                    //loci are unlinked:
                    for (y = 0; y < L; y++)
                    {
                        if (population[d][parent[0] * L + y] == population[d][parent[1] * L + y])
                        { //the parents are identical at the locus
                            population[d][i * L + y] = population[d][parent[0] * L + y];
                        }
                        else
                        { //the parents differ at the locus:
                            population[d][i * L + y] = gsl_ran_bernoulli(RNDM, 0.5);
                        }
                        mutcount[d][i] += population[d][i * L + y];
                    }
                }
                else
                { //if it is not a recombinant
                    parent[0] = parity * N + gsl_ran_discrete(RNDM, gametes[d]);
                    memcpy(population[d] + i * L, population[d] + parent[0] * L, L * sizeof(unsigned char));
                    mutcount[d][i] = mutcount[d][parent[0]];
                }
            }

            //now add mutations:
            //find position of first mutation:
            y = gsl_ran_geometric(RNDM, m) - 1; //gsl_ran_geometric has min value 1
            while (y < L)
            {
                population[d][i * L + y] = 1 - population[d][i * L + y]; //flip the value
                mutcount[d][i] += 2 * population[d][i * L + y] - 1;
                y += gsl_ran_geometric(RNDM, m); //find position of next mutation
            }

            //update fitness:
            W[d][i - (1 - parity) * N] = exp(sels[mutcount[d][i]]);
        }
    }
    //Delete used gamete pool:
    for (d = 0; d < D; d++)
    {
        gsl_ran_discrete_free(gametes[d]);
    }

    return EXIT_SUCCESS;
}

/********************************************/
//adding data to a dynamic array whose size might not be enough
int add_data(unsigned int **deme, int d, unsigned int **ind, int i, unsigned int **cap_size, unsigned int parity, unsigned int init_cap)
{
    unsigned int *deme_new, *ind_new;
    if (cap_size[1 - parity][0] == cap_size[1 - parity][1])
    {
        deme_new = realloc(deme[1 - parity], (cap_size[1 - parity][0] + init_cap) * sizeof(int)); //increase the capsity of arr1
        ind_new = realloc(ind[1 - parity], (cap_size[1 - parity][0] + init_cap) * sizeof(int));   //increase the capsity of arr2
        if (deme_new && ind_new)
        {
            deme[1 - parity] = deme_new;
            deme[1 - parity][cap_size[1 - parity][1]] = d;
            ind[1 - parity] = ind_new;
            ind[1 - parity][cap_size[1 - parity][1]] = i;
            cap_size[1 - parity][0] += init_cap;
            cap_size[1 - parity][1]++;
        }
        else
        {
            printf("Realloc memory failed\n");
            return EXIT_FAILURE;
        }
    }
    else
    {
        deme[1 - parity][cap_size[1 - parity][1]] = d;
        ind[1 - parity][cap_size[1 - parity][1]] = i;
        cap_size[1 - parity][1]++;
    }
    return EXIT_SUCCESS;
}