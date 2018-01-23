using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace PathFinder
{
    class GeneticAlgorithm
    {

        // variables
        int[,] population;  // one row for each chromosome, 
        double[] fitness;
        int populationSize = 100;  // default values
        int maxSteps = 35;
        int tnSize = 2;
        Random randomNum;
        int[] child;
        double mutationRate;
        int length;


        public GeneticAlgorithm(int p, int t, double m, int ms)
        {
            // constructor sets up variables
            // population size, tournament size, mutation rate,
            // chromosome length (= 2 times max-steps as 2 bits for each step)

            populationSize = p;  // population size
            tnSize = t;          // tournament size
            maxSteps = ms;       // maximum steps 
            mutationRate = m;      // mutation rate   
            length = maxSteps * 2;   // chromosome length
            population = new int[populationSize, length];
            fitness = new double[populationSize];
            child = new int[length];  // placeholder for child
            randomNum = new Random();

        }


        public void initialise()
        {
            // initialise each string with 0,1s at random
            int j, k;
            for (j = 0; j < populationSize; j++)
                for (k = 0; k < length; k++)
                    population[j, k] = randomNum.Next(2);

        }

        public int select()
        {
            // tournament selection picks x random chromosomes,
            // and returns the fittest one

            int i, picked, bestIndex;
            double bestFit;

            // pick one
            picked = randomNum.Next(populationSize);
            bestFit = fitness[picked];
            bestIndex = picked;

            // pick tnSize-1 more and then see which is the best
            for (i = 0; i < tnSize - 1; i++)
            {
                picked = randomNum.Next(populationSize);
                if (fitness[picked] > bestFit)
                {
                    bestFit = fitness[picked];
                    bestIndex = picked;
                }
            }
            return (picked);	// return index of best one

        }

        // crossover 
        public int[] onept_crossover(int p1, int p2)
        {
            // one point
            int crossPoint, j;

            // pick crosspoint at random
            crossPoint = randomNum.Next(length);

            // first part of child is copied from p1, and 2nd
            // part from p2
            for (j = 0; j < crossPoint; j++)
                child[j] = population[p1, j];
            for (j = crossPoint; j < (length); j++)
                child[j] = population[p2, j];
            return (child);
        }

        public int[] twopt_crossover(int p1, int p2)
        {
            int j;
            // not coded yet ... just returns p1
            for (j = 0; j < length; j++)
                child[j] = population[p1, j];
            return (child);
        }

        public int[] uniform_crossover(int p1, int p2)
        {
            int j;
            // not coded yet ... just returns p1
            for (j = 0; j < length; j++)
                child[j] = population[p1, j];
            return (child);
        }

        // mutate
        public int[] mutate(int[] child)
        {
            // mutate each gene with a small probability
            for (int j = 0; j < length; j++)
            {
                if (randomNum.NextDouble() < mutationRate)
                {
                    if (child[j] == 0)
                    {
                        child[j] = 1;
                    }
                    else
                    {
                        child[j] = 0;
                    }
                }
            }
            return (child);
        }

        // myMutate:  add your own function here
        public int[] myMutate(int[] child)
        {

            for (int j = 0; j < length; j += 2)
            {
                if (randomNum.NextDouble() < mutationRate)
                {
                    // mutate pairs of genes
                    //child[j] = randomNum.Next(2);
                    //child[j + 1] = randomNum.Next(2);
                    if (j < length - 2)
                    {
                        if (((child[j] == 0 && child[j + 1] == 1) && (child[j + 2] == 0 && child[j + 3] == 0)) //south - north
                           || ((child[j] == 0 && child[j + 1] == 0) && (child[j + 2] == 0 && child[j + 3] == 1))) // north - south
                        {
                            switch (randomNum.Next(2))
                            {
                                case 0:
                                    child[j] = randomNum.Next(2);
                                    child[j + 1] = randomNum.Next(2);
                                    break;
                                case 1:
                                    child[j + 2] = randomNum.Next(2);
                                    child[j + 3] = randomNum.Next(2);
                                    break;
                            }
                        }
                        else if (((child[j] == 1 && child[j + 1] == 0) && (child[j + 2] == 1 && child[j + 3] == 1)) //east - west
                                 || ((child[j] == 1 && child[j + 1] == 1) && (child[j + 2] == 1 && child[j + 3] == 0))) // west - east
                        {
                            switch (randomNum.Next(2))
                            {
                                case 0:
                                    child[j] = randomNum.Next(2);
                                    child[j + 1] = randomNum.Next(2);
                                    break;
                                case 1:
                                    child[j + 2] = randomNum.Next(2);
                                    child[j + 3] = randomNum.Next(2);
                                    break;
                            }
                        }
                    }
                }
            }
            return (child);
        }

        // replace
        public void replace(int[] child, double childFitness)
        {

            // find the worst in the population
            // replace the worst with the child
            // (as long as the child is better than the worst
            // (if not, child is just ignored)

            double worst = fitness[0];
            int worstIndex = 0;

            for (int j = 0; j < populationSize; j++)
                if (fitness[j] < worst)
                {
                    worstIndex = j;
                    worst = fitness[j];
                }

            // now replace worst with child if child is better
            if (childFitness > worst)
            {
                for (int j = 0; j < length; j++)
                    population[worstIndex, j] = child[j];
                fitness[worstIndex] = childFitness;
            }

        }

        public void myReplace(int[] child, double childFitness)
        {
            // empty method for you to add your own mutation
            // function
            // for now, this just calls replace, you will need
            // to remove this line
            replace(child, childFitness);
        }


        public double evaluate(int x, int y, int endX, int endY)
        {
            // fitness is the sum of the difference between the
            // current X and Y position, and the X,Y of the exit
            // max fitness is 1 when X's and Y's are equal
            int DiffX, DiffY;

            DiffX = Math.Abs(x - endX);
            DiffY = Math.Abs(y - endY);
            return (1 / (double)(DiffX + DiffY + 1));
        }


        //******************************************************
        //  The method below are simply used to pass information
        //  between the main application and the EA.
        //  You should NOT need to modify any of these methods
        //******************************************************


        public int[] getRoute(int index)
        {
            // this returns an integer array containing 
            // the bit string at the index specified by "index" in the population
            int[] route = new int[length];
            for (int j = 0; j < length; j++)
                route[j] = population[index, j];
            return (route);
        }

        public void setFitness(int index, double f)
        {
            // set the fitness of the chromosome with index "index"
            fitness[index] = f;
        }
        public double returnFitness(int index)
        {
            // return the fitness of a given chromosome to the main program
            return (fitness[index]);
        }
        public int getBest()
        {
            // returns the INDEX in the population of the best chromosome
            int bestIndex = 0;
            double best = fitness[0];
            for (int j = 0; j < populationSize; j++)
                if (fitness[j] > best)
                {
                    bestIndex = j;
                    best = fitness[j];
                }
            return (bestIndex);
        }

        public double getBestFitness()
        {
            // returns the FITNESS of the best chromosome
            int bestIndex = 0;
            double best = fitness[0];
            for (int j = 0; j < populationSize; j++)
                if (fitness[j] > best)
                {
                    bestIndex = j;
                    best = fitness[j];
                }
            return (best);
        }
    }
}
