/*
* memetic programming
* Zhangli Wang, 20028336
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <sys/time.h>

// size of population
#define POPULATION_SIZE 10

// probability of performing crossover
#define CROSSOVER_PROBABILITY 0.8

// probability of performing uniform crossover
#define UNICROSSOVER_PROBABILITY 0.5

// probability of mutation
#define MUTATION_PROBABILITY 0.01

// maximum time of iteration of generation
#define MAX_ITER 1000000

// max iteration of vns local search
#define VNS_MAX_ITER 1

/* global parameters */
struct solution_struct best_sln;
int MAX_TIME = 300; // max amount of time for each iteration (sec)
int num_of_problems;
int RAND_SEED[] = {1,20,30,40,50,60,70,80,90,100,110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
int NUM_OF_RUNS = 1;
int CROSSOVER_MARK = 0; // determine crossover method
int K = 3; // VNS local search

/* main data structure */
struct item_struct {
    int dimN; // number of dimentions
    int* size; // number of item in all dimensions
    int p; // profit of item
    int indx; // index of the item in a problem
    double u; // profit per size, calculated by p/total size
};

struct problem_struct {
    int itemN; // number of items
    struct item_struct* items; // problem items
    int dimN; // number of dimensions
    int* capacities; // knapsack capacities
};

struct solution_struct {
    struct problem_struct* prob; // maintain a shallow copy of problem structure
    int* x; // chromosome
    int* cap_left; // capacities left in all dimensions
    int objective; // objective value
    int feasibility; // the feasibility of the solution
    double selectProbability; // probability of selected by roulette wheel
};

/* assistant functions */
//return a random number between 0 and 1
float rand_01()
{
    float number;
    number = (float) rand();
    number = number/RAND_MAX;
    //printf("rand01=%d\n", number);
    return number;
}

//return a random nunber ranging from min to max (inclusive)
int rand_int(int min, int max)
{
    int div = max-min+1;
    int val =rand() % div + min;
    //printf("rand_range= %d \n", val);
    return val;
}

/* problem assistant functions */
// init problem
void init_problem(int n, int dim, struct problem_struct** my_prob)
{
    struct problem_struct* new_prob = malloc(sizeof(struct problem_struct));
    new_prob->itemN=n; new_prob->dimN=dim;
    new_prob->items=malloc(sizeof(struct item_struct)*n);
    for(int j=0; j<n; j++)
        new_prob->items[j].size= malloc(sizeof(int)*dim);
    new_prob->capacities = malloc(sizeof(int)*dim);
    *my_prob = new_prob;
}

// free problem
void free_problem(struct problem_struct* prob)
{
    if(prob!=NULL)
    {
        if(prob->capacities !=NULL) free(prob->capacities);
        if(prob->items!=NULL)
        {
            for(int j=0; j<prob->itemN; j++)
            {
                if(prob->items[j].size != NULL)
                    free(prob->items[j].size);
            }
            free(prob->items);
        }
        free(prob);
    }
}

// example to create problem instances, actual date should come from file
struct problem_struct** load_problems(char* data_file)
{
    int i,j,k;
    //int num_of_probs;
    FILE* pfile = fopen(data_file, "r");
    if(pfile==NULL)
        {printf("Data file %s does not exist. Please check!\n", data_file); exit(2); }
    fscanf (pfile, "%d", &num_of_problems);
 
    struct problem_struct** my_problems = malloc(sizeof(struct problem_struct*)*num_of_problems);
    for(k=0; k<num_of_problems; k++)
    {
        int n, dim, obj_opt;
        fscanf (pfile, "%d", &n);
        fscanf (pfile, "%d", &dim); fscanf (pfile, "%d", &obj_opt);
        
        init_problem(n, dim, &my_problems[k]);  //allocate data memory
        for(j=0; j<n; j++)
        {
            my_problems[k]->items[j].dimN=dim;
            my_problems[k]->items[j].indx=j;
            fscanf(pfile, "%d", &(my_problems[k]->items[j].p)); //read profit data
            //printf("item[j].p=%d\n",my_problems[k]->items[j].p);
        }
        for(i=0; i<dim; i++)
        {
            for(j=0; j<n; j++)
            {
                fscanf(pfile, "%d", &(my_problems[k]->items[j].size[i])); //read size data
                //printf("my_problems[%i]->items[%i].size[%i]=%d\n",k,j,i,my_problems[k]->items[j].size[i]);
            }
        }
        for(i=0; i<dim; i++){
            fscanf(pfile, "%d", &(my_problems[k]->capacities[i]));
            //printf("capacities[i]=%d\n",my_problems[k]->capacities[i] );
        }
    }
    fclose(pfile); //close file
    return my_problems;
}

/* solution assistant functions */
// free solution
void free_solution(struct solution_struct* sln)
{
    if(sln!=NULL)
    {
        free(sln->x);
        free(sln->cap_left);
        sln->objective=0;
        sln->prob=NULL;
        sln->feasibility=false;
    }
    free(sln);
}

// copy a solution from another solution
bool copy_solution(struct solution_struct* dest_sln, struct solution_struct* source_sln)
{
    if(source_sln ==NULL) return false;
    if(dest_sln==NULL)
    {
        dest_sln = malloc(sizeof(struct solution_struct));
    }
    else{
        free(dest_sln->cap_left);
        free(dest_sln->x);
    }
    int n = source_sln->prob->itemN;
    int m =source_sln->prob->dimN;
    dest_sln->x = malloc(sizeof(int)*n);
    dest_sln->cap_left=malloc(sizeof(int)*m);
    for(int i=0; i<m; i++)
        dest_sln->cap_left[i]= source_sln->cap_left[i];
    for(int j=0; j<n; j++)
        dest_sln->x[j] = source_sln->x[j];
    dest_sln->prob= source_sln->prob;
    dest_sln->feasibility=source_sln->feasibility;
    dest_sln->objective=source_sln->objective;
    return true;
}

//output a given solution to a file
void output_solution(struct solution_struct* sln, char* out_file)
{
    if(out_file !=NULL){
        FILE* pfile = fopen(out_file, "a"); //append solution data
        fprintf(pfile, "%i\n", (int)sln->objective);
        for(int i=0; i<sln->prob->itemN; i++)
        {
            fprintf(pfile, "%i ", sln->x[i]);
        }
        fprintf(pfile, "\n");
        /*for(int j=0; j<sln->prob->n; j++)
            fprintf(pfile, "%i ", sln->prob->items[j].p);
        fprintf(pfile, "\n");*/
        fclose(pfile);
    }
    else
        printf("sln.feas=%d, sln.obj=%f\n", sln->feasibility, sln->objective);
}

// compare unit profit of two items
int cmpfunc(const void* a, const void* b) {
    const struct item_struct* item1 = a;
    const struct item_struct* item2 = b;
    if(item1->u<item2->u) return -1;
    if(item1->u>item2->u) return 1;
    return 0;
}

//update global best solution from sln
void update_best_solution(struct solution_struct* sln)
{
    if(best_sln.objective < sln->objective)
    copy_solution(&best_sln, sln);
}

// init a empty solution
struct solution_struct* initEmptySolution(struct problem_struct* prob) {
    struct solution_struct* initialS = malloc(sizeof(struct solution_struct));
    initialS->cap_left = malloc(sizeof(int)*prob->dimN);
    initialS->x = malloc(sizeof(int)*prob->itemN);
    initialS->feasibility = 1;
    initialS->objective = 0;
    initialS->prob = prob;
    return initialS;
}

// initialize a solution
struct solution_struct* initSolution(struct problem_struct* prob) {
    // memory allocation and initialization
    struct solution_struct* initialS = malloc(sizeof(struct solution_struct));
    initialS->prob = prob;
    initialS->objective = 0;
    initialS->x = malloc(sizeof(int)*prob->itemN);
    initialS->cap_left = malloc(sizeof(int)*prob->dimN);
    initialS->feasibility = 1;

    // assign values
    int i = 0;
    int j = 0;
    bool mark = true;
    int* cap = malloc(sizeof(int)*prob->dimN);
    for (j = 0; j < prob->dimN; j++)
        cap[j] = 0;

    // unpack all items
    for (i = 0; i < prob->itemN; i++)
    {
        initialS->x[i] = 0;
    }

    // try to pack an item
    for (i = 0; i < prob->itemN; i++)
    {
        // avoid repeated selection
        int select = rand_int(0, prob->itemN-1);
        if (initialS->x[select] == 1)
        {
            i--;
            continue;
        }
        else
        {
            // check whether the item can be packed
            struct item_struct* item_select = &prob->items[select]; 
            for (j = 0; j < prob->dimN; j++)
            {
                if (cap[j]+item_select->size[j] > prob->capacities[j])
                {
                    mark = false;
                    break;
                }
            }

            if (mark)
            {
                initialS->x[select] = 1;
                // update capacity
                for (j = 0; j < prob->dimN; j++)
                {
                    cap[j] += item_select->size[j];
                }
                // update objective value
                initialS->objective += item_select->p;
            }
        }
    }
    

    // update left capacity
    for (j = 0; j < prob->dimN; j++)
    {
        initialS->cap_left[j] = prob->capacities[j]-cap[j];
    }

    // free memory
    free(cap);

    return initialS;
}

/* GA assistant functions */
// initialize a population
struct solution_struct** initPopulation(struct problem_struct* prob) {
    struct solution_struct** initialPop = malloc(sizeof(struct solution_struct*)*POPULATION_SIZE);
    int i;
    for (i = 0; i < POPULATION_SIZE; i++)
    {
        initialPop[i] = initSolution(prob);
        update_best_solution(initialPop[i]);
    }
    return initialPop;
}


// roulette wheel selection
struct solution_struct* wheelSelection(struct solution_struct** pop) {
    double offset = 0.0;
    int pick = 0;

    int totalFitness = 0;
    int i;
    // calculate total fitness
    for (i = 0; i < POPULATION_SIZE; i++)
    {
        totalFitness += pop[i]->objective;
    }

    // calculate select probability of each solution
    for (i = 0; i < POPULATION_SIZE; i++)
    {
        pop[i]->selectProbability = (double) (pop[i]->objective)/totalFitness;
    }
    
    double rndNumber = rand()/(RAND_MAX+1.0);

    // pick the solution
    for (i = 0; i < POPULATION_SIZE; i++)
    {
        offset += pop[i]->selectProbability;
        if (rndNumber < offset)
        {
            pick = i;
            break;
        }
    }
    return pop[pick];
}

// modify the solutions that violate the capacity constraints
void feasibility_repair(struct solution_struct* sln)
{
    // initialize a list of pointers to all items in the problem
    struct item_struct** items_sorted = malloc(sizeof(struct item_struct*)*sln->prob->itemN);
    for (int i = 0; i < sln->prob->itemN; i++)
    {
        items_sorted[i] = &(sln->prob->items[i]);
    }
    
    // calculate total size and update u of items
    for (int i = 0; i < sln->prob->itemN; i++)
    {
        // calculate total size of an item
        int totalSize = 0;
        for (int j = 0; j < sln->prob->dimN; j++)
        {
            totalSize += sln->prob->items[i].size[j];
        }

        // calculate unit profit of an item
        sln->prob->items[i].u = (double)(sln->prob->items[i].p/totalSize);
    }

    // sort items according to u in desent order
    qsort(items_sorted, sln->prob->itemN, sizeof(struct item_struct*), cmpfunc);
    
    // drop items packed with lowest u
    for (int i = sln->prob->itemN-1; i >= 0; i--)
    {
        if (sln->x[items_sorted[i]->indx] == 1)
        {
            // check feasibility
            for (int j = 0; j < sln->prob->dimN; j++)
            {
                if (sln->cap_left[j] < 0)
                {
                    // unpack item, update solution
                    sln->x[items_sorted[i]->indx] = 0;
                    for (int z = 0; z < sln->prob->dimN; z++)
                    {
                        sln->cap_left[z] += sln->prob->items[items_sorted[i]->indx].size[z];
                    }
                    sln->objective -= sln->prob->items[items_sorted[i]->indx].p;
                    break;
                }
            }
        }
    }
     
    // add items unpacked with highest u
    for (int i = 0; i < sln->prob->itemN; i++)
    {
        if (sln->x[items_sorted[i]->indx] == 0)
        {
            // check feasibility
            bool count = true;
            for (int j = 0; j < sln->prob->dimN; j++)
            {
                if (sln->cap_left[j]-sln->prob->items[items_sorted[i]->indx].size[j] < 0)
                {
                    count = false;
                }
            }

            // pack item, update solution
            if (count)
            {
                sln->x[items_sorted[i]->indx] = 1;
                for (int z = 0; z < sln->prob->dimN; z++)
                {
                    sln->cap_left[z] -= sln->prob->items[items_sorted[i]->indx].size[z];
                }
                sln->objective += sln->prob->items[items_sorted[i]->indx].p;
            }
        }
        
    }
    
    sln->feasibility = 1;

    // free memory
    free(items_sorted);
}

// crossover two offsprings
void crossover_1p(struct solution_struct* offspring1, struct solution_struct* offspring2) {
    // determine whether to cut
    float rndCrossover = rand_01();
    if (rndCrossover > CROSSOVER_PROBABILITY)
    {
        return;
    }
    
    // determine cut point
    int rndCutpoint = rand_int(0, offspring1->prob->dimN-1);

    // crossover
    int i;
    for (i = 0; i < offspring1->prob->itemN; i++)
    {
        if (i > rndCutpoint)
        {
            int z;
            struct item_struct* item_i = &offspring1->prob->items[i];
            if (offspring1->x[i] == 0)
            {
                if (offspring2->x[i] == 1)
                {
                    offspring1->objective += item_i->p;
                    offspring2->objective -= item_i->p;
                    for (z = 0; z < offspring1->prob->dimN; z++)
                    {
                        offspring1->cap_left[z] -= item_i->size[z];
                        offspring2->cap_left[z] += item_i->size[z];
                    }
                }
            }
            else
            {
                if (offspring2->x[i] == 0)
                {
                    offspring1->objective -= item_i->p;
                    offspring2->objective += item_i->p;
                    for (z = 0; z < offspring1->prob->dimN; z++)
                    {
                        offspring1->cap_left[z] += item_i->size[z];
                        offspring2->cap_left[z] -= item_i->size[z];
                    }
                }
                
            }
            int tempx = offspring1->x[i];
            offspring1->x[i] = offspring2->x[i];
            offspring2->x[i] = tempx;
        }
        
    }

    // check feasibility
    for (i = 0; i < offspring1->prob->dimN; i++)
    {
        if (offspring1->cap_left[i] < 0)
        {
            offspring1->feasibility = 0;
            feasibility_repair(offspring1);
            //offspring1->objective = 0;
            break;
        }
    }

    for (i = 0; i < offspring1->prob->dimN; i++)
    {
        if (offspring2->cap_left[i] < 0)
        {
            offspring2->feasibility = 0;
            feasibility_repair(offspring2);
            //offspring2->objective = 0;
            break;
        }
    }
}

// uniform crossover function
void crossover_uniform(struct solution_struct* offspring1, struct solution_struct* offspring2) {

    // determine whether to cut
    float rndCrossover = rand_01();
    if (rndCrossover > CROSSOVER_PROBABILITY)
    {
        return;
    }

    for (int i = 0; i < offspring1->prob->itemN; i++)
    {
        // determine whether to crossover for a point
        float rndUniCrossover = rand_01();
        if (rndUniCrossover < UNICROSSOVER_PROBABILITY)
        {
            int z;
            // perform crossover
            struct item_struct* item_i = &offspring1->prob->items[i];
            if (offspring1->x[i] == 0)
            {
                if (offspring2->x[i] == 1)
                {
                    offspring1->objective += item_i->p;
                    offspring2->objective -= item_i->p;
                    for (z = 0; z < offspring1->prob->dimN; z++)
                    {
                        offspring1->cap_left[z] -= item_i->size[z]; 
                        offspring2->cap_left[z] += item_i->size[z]; 
                    }
                }
            }
            else
            {
                if (offspring2->x[i] == 0)
                {
                    offspring1->objective -= item_i->p;
                    offspring2->objective += item_i->p;
                    for (z = 0; z < offspring1->prob->dimN; z++)
                    {
                        offspring1->cap_left[z] += item_i->size[z]; 
                        offspring2->cap_left[z] -= item_i->size[z]; 
                    }
                }
                
            }
            int tempx = offspring1->x[i];
            offspring1->x[i] = offspring2->x[i];
            offspring2->x[i] = tempx;            
        }
        
    }

    // check feasibility
    for (int i = 0; i < offspring1->prob->dimN; i++)
    {
        if (offspring1->cap_left[i] < 0)
        {
            offspring1->feasibility = 0;
            feasibility_repair(offspring1);
            //offspring1->objective = 0;
            break;
        }
    }

    for (int i = 0; i < offspring1->prob->dimN; i++)
    {
        if (offspring2->cap_left[i] < 0)
        {
            offspring2->feasibility = 0;
            feasibility_repair(offspring2);
            //offspring2->objective = 0;
            break;
        }
    }
}


// mutation function
void mutate(struct solution_struct* chromosome) {
    int i;
    for (i = 0; i < chromosome->prob->itemN; i++)
    {
        float rndNumber = rand_01();
        if (rndNumber < MUTATION_PROBABILITY)
        {
            struct item_struct* item_i = &chromosome->prob->items[i];
            int j;
            if (chromosome->x[i] == 0)
            {
                // update solution
                chromosome->x[i] = 1;
                chromosome->objective += item_i->p;
                for (j = 0; j < chromosome->prob->dimN; j++)
                {
                    chromosome->cap_left[j] -= item_i->size[j];
                    // feasibility check
                    if (chromosome->cap_left[j] < 0)
                    {
                        chromosome->feasibility = 0;
                    }
                }

                if (chromosome->feasibility == 0)
                {
                    feasibility_repair(chromosome);
                    //chromosome->objective = 0;
                }
                
            }
            else
            {
                // update solution
                chromosome->x[i] = 0;
                chromosome->objective -= item_i->p;
                for (j = 0; j < chromosome->prob->dimN; j++)
                {
                    chromosome->cap_left[j] += item_i->size[j];
                }                
            }
        }
        
    }
    
}


/* VNS assistant functions*/
bool can_swap(struct solution_struct* sln, int out, int in)
{
    for(int d =0; d<sln->prob->dimN; d++)
    {
        if(sln->cap_left[d]+sln->prob->items[out].size[d] < sln->prob->items[in].size[d])
            return false;
    }
    return true;
}

bool can_move(int nb_indx, int* move, struct solution_struct* curt_sln ){
    bool ret=true;
    if(nb_indx==1)
    {
        int i = move[0];
        if(i<0) return false;
        for(int d=0; d<curt_sln->prob->dimN; d++){
            if(curt_sln->cap_left[d] < curt_sln->prob->items[i].size[d])
                return false;
        }
    }
    else if(nb_indx==2){
        ret=can_swap(curt_sln, move[0], move[1]);
    }
    else if(nb_indx==3){//3-item swap
        int i= move[0], j= move[1], k= move[2];
        if(i<0 || j<0 || k<0) return false;
        if(curt_sln->x[j]>0) {//2-1 swap
            for(int d=0; d<curt_sln->prob->dimN; d++){
                if(curt_sln->cap_left[d] + curt_sln->prob->items[i].size[d] +
                   curt_sln->prob->items[j].size[d] < curt_sln->prob->items[k].size[d])
                    return false;
            }
        }
        else {//1-2 swap
            for(int d=0; d<curt_sln->prob->dimN; d++){
                if(curt_sln->cap_left[d] + curt_sln->prob->items[i].size[d] <
                   curt_sln->prob->items[j].size[d] + curt_sln->prob->items[k].size[d])
                    return false;
            }
        }
        
    }
    else ret=false;
    return ret;
}

bool apply_move(int nb_indx, int* move, struct solution_struct* sln ){
    bool ret=true;
    if(nb_indx==1) // add items
    {
        int i = move[0];
        if(i<0) return false;
        for(int d=0; d<sln->prob->dimN; d++){
            sln->cap_left[d] -= sln->prob->items[i].size[d];
        }
        sln->objective += sln->prob->items[i].p;
        sln->x[i]=1;
    }
    else if(nb_indx==2){ // 1-1 swap
        for(int d=0; d<sln->prob->dimN; d++){
            sln->cap_left[d] = sln->cap_left[d] + sln->prob->items[move[0]].size[d]-
                sln->prob->items[move[1]].size[d];
        }
        sln->objective += sln->prob->items[move[1]].p-sln->prob->items[move[0]].p;
        sln->x[move[0]]=0; sln->x[move[1]]=1;
    }
    else if(nb_indx==3){//3-item swap
        int i= move[0], j= move[1], k= move[2];
        if(i<0 || j<0 || k<0) return false;
        if(sln->x[j]>0) {//2-1 swap
            for(int d=0; d<sln->prob->dimN; d++){
                sln->cap_left[d] = sln->cap_left[d]+sln->prob->items[i].size[d] +
                    sln->prob->items[j].size[d] - sln->prob->items[k].size[d];
            }
            sln->objective += sln->prob->items[k].p-sln->prob->items[i].p-sln->prob->items[j].p;
            sln->x[i]=0; sln->x[j]=0; sln->x[k]=1;
        }
        else {//1-2 swap
            for(int d=0; d<sln->prob->dimN; d++){
                sln->cap_left[d] = sln->cap_left[d]+sln->prob->items[i].size[d] -
                    sln->prob->items[j].size[d] - sln->prob->items[k].size[d];
            }
            sln->objective += sln->prob->items[j].p+sln->prob->items[k].p-sln->prob->items[i].p;
            sln->x[i]=0; sln->x[j]=1; sln->x[k]=1;
        }
        
    }
    else ret=false;
    return ret;
}

/* VNS local search*/
// best descent vns local search
struct solution_struct* best_descent_vns(int nb_indx, struct solution_struct* curt_sln)
{
    bool mark = true; // first descent controller
    struct solution_struct* best_neighb = malloc(sizeof(struct solution_struct));
    best_neighb->cap_left = malloc(sizeof(int)*curt_sln->prob->dimN);
    best_neighb->x = malloc(sizeof(int)*curt_sln->prob->itemN);
    copy_solution(best_neighb, curt_sln);
    int n=curt_sln->prob->itemN;
    int curt_move[] ={-1,-1,-1}, best_move []={-1,-1,-1}, delta=0, best_delta=0;  //storing best neighbourhood moves
    switch (nb_indx)
    {
        
        case 1: 
            // insert an new item
            for(int i=0; i<n && mark; i++){
                if(curt_sln->x[i]>0) continue;
                curt_move[0]=i;
                if(can_move(nb_indx, &curt_move[0], best_neighb)){
                    delta = curt_sln->prob->items[i].p;
                    if(delta> best_delta) {
                        best_delta = delta; best_move[0] = i;
                        //mark = false;
                        //break;
                    }
                }
            }
            if(best_delta>0) {    apply_move(nb_indx, &best_move[0], best_neighb);}
            break;
            
        case 2:
            // 1-1 swap
            for(int i=0; i<n && mark; i++){
                if(curt_sln->x[i]<=0) continue;
                for(int j=0; j<n && mark; j++){
                    if(curt_sln->x[j]==0)
                    {
                        curt_move[0]= i; curt_move[1]= j; curt_move[2]=-1;
                        if(can_move(nb_indx, &curt_move[0], best_neighb)){
                            delta = curt_sln->prob->items[j].p -curt_sln->prob->items[i].p;
                            if(delta > best_delta){
                                best_delta = delta; best_move[0] = i; best_move[1] = j; best_move[2]=-1;
                                //mark = false;
                                //break;
                            }
                        }
                    }
                }
            }
            if(best_delta>0) { apply_move(nb_indx, &best_move[0], best_neighb);}
            break;
            
           
        case 3:

            //2-1 swap
            for(int i=0; i<n&&mark; i++){
                
                if(curt_sln->x[i]==0) continue;
                for(int j=0; j!=i&&j<n&&mark; j++){
                    //printf("mark %d\n", mark);
                    if(curt_sln->x[j]==0) continue;
                    for(int k=0;k<n&&mark;k++){
                        
                        if(curt_sln->x[k] == 0)
                        {
                            curt_move[0]=i; curt_move[1]=j; curt_move[2]=k;
                            if(can_move(nb_indx, &curt_move[0], best_neighb)){
                                delta = curt_sln->prob->items[k].p -curt_sln->prob->items[i].p-curt_sln->prob->items[j].p;
                                if(delta > best_delta){
                                    best_delta = delta; best_move[0] = i; best_move[1] = j; best_move[2]=k;
                                    //mark = false;
                                    //break;
                                }
                            }
                        }
                    }
                }
            }
            
            // reset first descent controller
            //mark = true;

            //1-2 swap
            for(int i=0; i<n&&mark; i++){
                if(curt_sln->x[i]==0) continue;
                for(int j=0; j<n&&mark; j++){
                    if(curt_sln->x[j]>0) continue;
                    for(int k=0;k!=j&&k<n&&mark;k++){
                        if(curt_sln->x[k] == 0)
                        {
                            curt_move[0]=i; curt_move[1]=j; curt_move[2]=k;
                            if(can_move(nb_indx, &curt_move[0], curt_sln)){
                                delta = curt_sln->prob->items[k].p +curt_sln->prob->items[j].p-curt_sln->prob->items[i].p;
                                if(delta > best_delta){
                                    best_delta = delta; best_move[0] = i; best_move[1] = j; best_move[2]=k;
                                    //mark = false;
                                    //break;
                                }
                            }
                        }
                    }
                }
            }
            if(best_delta>0) { apply_move(nb_indx, &best_move[0], best_neighb);}
            break;
            
        default:
            printf("Neighbourhood index is out of the bounds, nothing is done!\n");
    }
                
    return best_neighb;
}

//VNS
struct solution_struct* varaible_neighbourhood_search(struct solution_struct* curt_sln, double VNS_max_time){
    int iter = 0;
    int nb_indx = 0; //neighbourhood index
    
    while(VNS_max_time > 0.0 && iter < VNS_MAX_ITER) //note that final computational time can be way beyond the MAX_TIME if best_descent is time consuming
    {
        while(nb_indx<K){
            struct solution_struct* neighb_s=best_descent_vns(nb_indx+1, curt_sln); //best solution in neighbourhood nb_indx
            if(neighb_s->objective > curt_sln->objective){
                copy_solution(curt_sln, neighb_s);
                nb_indx=1;
            }
            else nb_indx++;
            free_solution(neighb_s);
        }
        update_best_solution(curt_sln);
        nb_indx=0;
        
        iter++;
    }
    return curt_sln;
}

// genetic algorithm
int MemeticAlgorithm(struct problem_struct* prob) 
{
    // start time count
    clock_t time_start, time_fin;
    time_start = clock();

    best_sln.prob = prob;

    // initialize population
    struct solution_struct** curtPop = initPopulation(prob);

    int t = 0; // number of iteration
    double time_spent = 0; // time spent for each iteration
    while (t <= MAX_ITER && time_spent < MAX_TIME)
    {
        int i;

        // perform selection
        struct solution_struct** parents = malloc(sizeof(struct solution_struct*)*POPULATION_SIZE);
        for (i = 0; i < POPULATION_SIZE; i++) {
            parents[i] = wheelSelection(curtPop);
        }

        // initialize offsprings by copying parents
        struct solution_struct** offsprings = malloc(sizeof(struct solution_struct*)*POPULATION_SIZE);
        for (i = 0; i < POPULATION_SIZE; i++)
        {
            offsprings[i] = initEmptySolution(prob);
        }
        
        for (i = 0; i < POPULATION_SIZE; i++) {
            if (!copy_solution(offsprings[i], parents[i]))
            {
                printf("selection: failed copy solution %d\n", i);
            }
        }

        
        // do crossovers for all offsprings
        for (i = 0; i < POPULATION_SIZE; i += 2)
        {   
            if (CROSSOVER_MARK == 0)
            {
                crossover_uniform(offsprings[i], offsprings[i+1]);
            }
            if (CROSSOVER_MARK == 1)
            {
                crossover_1p(offsprings[i], offsprings[i+1]);
            }
        }
        
        
        // do mutations for all offsprings
        for (i = 0; i < POPULATION_SIZE; i++)
        {
            mutate(offsprings[i]);
        }
        

        // replacement
        for (i = 0; i < POPULATION_SIZE; i++)
        {
            if (!copy_solution(curtPop[i], offsprings[i]))
            {
                printf("replacement: failed copy solution %d\n", i);
            }
            update_best_solution(curtPop[i]);
        }

        
        // local search, execute this step every three generations
        for (i = 0; i < POPULATION_SIZE; i++) {
            time_fin = clock();
            double remained_time = (double)MAX_TIME-(time_fin-time_start)/CLOCKS_PER_SEC;
                //printf("remainedtime %f\n", remained_time);
                curtPop[i] = varaible_neighbourhood_search(curtPop[i], remained_time);
                update_best_solution(curtPop[i]);
        }
        

        // free in-cycle memories
        free(parents);
        for (i = 0; i < POPULATION_SIZE; i++)
        {
            free_solution(offsprings[i]);
        }
        free(offsprings);        

        // print results
        //printf("iter: %d, bestsln INSIDE objective: %d\n", t, best_sln.objective);

        // update termination criterias
        t++;     
        time_fin = clock();
        time_spent = (double)(time_fin-time_start)/CLOCKS_PER_SEC;
        //printf("timespent %f\n", time_spent);
    }

    printf("bestsln objective: %d\n", best_sln.objective);
    
    // free out-cycle memories
    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        free_solution(curtPop[i]);
    }
    free(curtPop);    

    return 0;
}

int main(int argc, const char* argv[]) {
    clock_t main_time_start, main_time_fin;
    main_time_start = clock();
    printf("Starting the run!\n");
    char data_file[50]={"somefile"}, out_file[50]={}, solution_file[50]={};  //max 50 problem instances per run
    if(argc<3)
    {
        printf("Insufficient arguments. Please use the following options:\n   -s data_file (compulsory)\n   -o out_file (default my_solutions.txt)\n   -c solution_file_to_check\n   -t max_time (in sec)\n");
        return 1;
    }
    else if(argc>9)
    {
        printf("Too many arguments.\n");
        return 2;
    }
    else
    {
        for(int i=1; i<argc; i=i+2)
        {
            if(strcmp(argv[i],"-s")==0)
                strcpy(data_file, argv[i+1]);
            else if(strcmp(argv[i],"-o")==0)
                strcpy(out_file, argv[i+1]);
            else if(strcmp(argv[i],"-c")==0)
                strcpy(solution_file, argv[i+1]);
            else if(strcmp(argv[i],"-t")==0)
                MAX_TIME = atoi(argv[i+1]);
        }
        //printf("data_file= %s, output_file= %s, sln_file=%s, max_time=%d", data_file, out_file, solution_file, MAX_TIME);
    }
    struct problem_struct** my_problems = load_problems(data_file);
    
    if(strlen(solution_file)<=0)
    {
        if(strcmp(out_file,"")==0) strcpy(out_file, "my_solutions.txt"); //default output
        FILE* pfile = fopen(out_file, "w"); //open a new file
        fprintf(pfile, "%d\n", num_of_problems); fclose(pfile);
        for(int k=0; k<num_of_problems; k++)
        {
            printf("problem %d, ", k);
            best_sln.objective=0; best_sln.feasibility=0;
            for(int run=0; run<NUM_OF_RUNS; run++)
            {
                srand(RAND_SEED[run]);
                
                MemeticAlgorithm(my_problems[k]); // call GA method
            }
            output_solution(&best_sln,out_file);
        }
    }
    for(int k=0; k<num_of_problems; k++)
    {
       free_problem(my_problems[k]); //free problem data memory
    }
    free(my_problems); //free problems array
    if(best_sln.x!=NULL && best_sln.cap_left!=NULL){ free(best_sln.cap_left); free(best_sln.x);} //free global
    
    // time count
    main_time_fin = clock();
    printf("time spent: %f seconds\n", (double)(main_time_fin-main_time_start)/CLOCKS_PER_SEC);
    return 0;
}
