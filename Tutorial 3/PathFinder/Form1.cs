using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;



namespace PathFinder
{
    public partial class Form1 : Form
    {
        Boolean finishedGA = false;
        // MAZE DEFINITIONS      
        int mazeWidth;
        int mazeHeight;
        int rows = 10;
        int columns = 15;
       
        Graphics paper;
        int rectHeight;
        int rectWidth;
        internal int ExitX=14;
        internal int ExitY=7;
        int StartX=0;
        int StartY=2;

        // change the array if you want a new maze
        // make sure you use the codes in the variables below
        int WALL = 1;
        int EXIT = 5;
        int ENTRANCE = 8;
        int OPENSPACE = 0;
        int[,] theMaze = {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},{1,0,1,0,0,0,0,0,1,1,1,0,0,0,1},{8,0,0,0,0,0,0,0,1,1,1,0,0,0,1},{1,0,0,0,1,1,1,0,0,1,0,0,0,0,1},{1,0,0,0,1,1,1,0,0,0,0,0,1,0,1},{1,1,0,0,1,1,1,0,0,0,0,0,1,0,1},{1,0,0,0,0,1,0,0,0,0,1,1,1,0,1},{1,0,1,1,0,0,0,1,0,0,0,0,0,0,5},{1,0,1,1,0,0,0,1,0,0,0,0,0,0,1},{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}};

        // directions
        int NORTH = 0;
        int SOUTH = 1;
        int EAST = 2;
        int WEST = 3;

        // GA parameters
        GeneticAlgorithm myGA;
        int maxIterations = 200;
        int length;
        int populationSize;
    
        // types of crossover
        int ONE = 0;
        int TWO = 1;
        int UNIFORM = 2;

        // types of mutation
        int STANDARD = 0;
        int USER_M = 1;

        // types of replacement
        int WORST = 0;
        int USER_R = 1;

        // Player parameters
        int playerX;
        int playerY;
        Boolean finished=false;
        int[] bestPath; // keep a record of the best path

        public Form1()
        {
            
            InitializeComponent();
            // set up the Maze
            mazeHeight = mazePanel.Height;
            mazeWidth = mazePanel.Width;
            // take off a few pixels for around the edge
            // then divide by rows/columns to get rectangle size

            rectWidth = (int)(mazeWidth / columns);
            rectHeight = (int)(mazeHeight / rows);

            crossCombo.SelectedIndex = 0;
            mutateCB.SelectedIndex = 0;
            replaceCB.SelectedIndex = 0;

            paper = mazePanel.CreateGraphics();
            drawMaze();  // draws blank maze

        }
            

        protected override void OnPaint(PaintEventArgs e)
        {
            base.OnPaint(e);
            drawMaze();
            if (finishedGA == true)
                applyMoves(bestPath);
        }
        private void drawMaze()
        {
            int x, y,k,j;
            SolidBrush myBrush  = new SolidBrush(Color.Black);

            x=0;  // offset
            y=0;

            for (j = 0; j < rows; j++)
            {
                y = (j * rectHeight);
                x = 0;
                for (k = 0; k < columns; k++)
                {
                    x = k*rectWidth;
                  
                    if (theMaze[j, k] == WALL)
                        myBrush.Color = Color.Black;
                    else if (theMaze[j, k] == OPENSPACE)
                        myBrush.Color = Color.Yellow;
                    else if (theMaze[j, k] == ENTRANCE)
                        myBrush.Color = Color.Red;
                    else if (theMaze[j, k] == EXIT)
                        myBrush.Color = Color.Red;

                    paper.FillRectangle(myBrush, x, y, rectWidth, rectHeight);

                }
               
            }
           
           
        }

        private void updateScreen(){
             SolidBrush myBrush  = new SolidBrush(Color.Green);
             int x, y;
            // update new position 
       
            x=playerX*rectWidth;
            y=playerY*rectHeight;
            paper.FillRectangle(myBrush, x, y, rectWidth, rectHeight);
    
        }
        private void goButton_Click(object sender, EventArgs e)
        {
            // this is the code to set up and run the genetic algorithm
            int crossoverType;   // what crossover
            int mutationType;
            int replaceType;

            int j,p1,p2;
            double routeFitness;  
            int best,iters;         
            int [] child=new int[length];
            int []route=new int[length];
            int p,ms,t;
            double m;
            double fitness;

            // get parameter values from GUI
            p = (int) popSizeUpDown.Value;
            t = (int) tnSizeUpDown.Value;
            ms = (int) maxStepsUpDown.Value;
            maxIterations = (int)itersUpDown.Value;
            m = Convert.ToDouble(mutRateText.Text);
            length = ms * 2;
            crossoverType = crossCombo.SelectedIndex;
            populationSize = p;
            mutationType = mutateCB.SelectedIndex;
            replaceType = replaceCB.SelectedIndex;

            bestPath = new int[length];
            finished = false;

            // create a new instance of the GA
            myGA = new GeneticAlgorithm(p,t,m,ms);

            // initialise GA
            myGA.initialise();   // create random population

            
            for (j=0;j<populationSize;j++){  
                route = myGA.getRoute(j);  // get a chromosome
                applyMoves(route);  // apply the moves
                routeFitness = myGA.evaluate(playerX, playerY, ExitX, ExitY); // evaluate the initial population
                myGA.setFitness(j,routeFitness);  // set the fitness
            }

            // MAIN LOOP WHICH RUNS THE EA       
            iters = 0;
            while (iters<maxIterations && finished==false)
            {
                p1 = myGA.select();
                p2 = myGA.select();

                // crossover
                if (crossoverType == ONE)
                    child = myGA.onept_crossover(p1, p2);
                else if (crossoverType == TWO)
                    child = myGA.twopt_crossover(p1, p2);
                else
                    child = myGA.uniform_crossover(p1, p2);

                if (mutationType == STANDARD)
                    myGA.mutate(child);
                else
                    myGA.myMutate(child);

                // apply moves determined by child
                applyMoves(child);
               
                // evaluate the end position
                fitness = myGA.evaluate(playerX,playerY,ExitX,ExitY);

                if (replaceType == WORST)
                    myGA.replace(child, fitness);
                else
                    myGA.myReplace(child, fitness);
       
                iters++;
            }
            // at end, update screen with best
            best = myGA.getBest();
            bestPath = myGA.getRoute(best);
            finishedGA = true;
            string path = decodeChromo(bestPath);
            bestChromo.Text = path;
            applyMoves(bestPath);
            fitness = myGA.getBestFitness();
            fitnessLbl.Text = String.Format("{0:0.00}", fitness); 
            iterationsLbl.Text = iters.ToString();
        }


        private void applyMoves(int[] route){
            int direction;
            string code;
            int j;
            // START AT THE ENTRANCE
            playerX = StartX;
            playerY = StartY;
            drawMaze();  // start with a clear maze
            updateScreen();
                                                                                                                 
            // apply moves and see where player ends up
            for (j = 0; j < length; j += 2)
            {
                code = route[j].ToString() + route[j + 1].ToString();
                direction = decode(code);
                if (direction == NORTH)
                {
                    // check its OK first)
                    if (playerY>=1)
                        if (permittedMove(playerY-1,playerX))
                              playerY--;
                }
                else if (direction == SOUTH)
                {
                    if (playerY <= rows-2)
                        if (permittedMove(playerY + 1, playerX))
                            playerY++;                            
                }
                else if (direction == EAST)
                {
                    if (playerX<=(columns-2))
                        if (permittedMove(playerY, playerX + 1))
                             playerX++;
                }
                else if (direction == WEST)
                {
                    if (playerX>=1)
                        if (permittedMove(playerY, playerX - 1))
                            playerX--;
                }
               
                // update screen with the path
                updateScreen();

                // check if this is the exit
                if (theMaze[playerY, playerX] == EXIT)
                    finished = true;


              }
        }

        public Boolean permittedMove(int a, int b)
        {
            if (theMaze[a, b] == ENTRANCE || theMaze[a, b] == OPENSPACE || theMaze[a,b] == EXIT)
                return (true);
            else
                return (false);

        }


     public int decode(string code)
        {
            if (code.Equals("00"))
                return (NORTH);
            else if (code.Equals("01"))
                return (SOUTH);
            else if (code.Equals("10"))
                return (EAST);
            else 
                return (WEST);  // default
        }


     public string decodeChromo(int[] route)
     {
         string code;
         string path = "";
         int j, direction;
         for (j = 0; j < length; j += 2)
         {
             code = route[j].ToString() + route[j + 1].ToString();
             direction = decode(code);
             if (direction == NORTH) path += "N";
             else if (direction == SOUTH) path += "S";
             else if (direction == EAST) path += "E";
             else path += "W";
         }
         return (path);
     }

        private void mazePanel_Paint(object sender, PaintEventArgs e)
        {

        }

        private void resetButton_Click(object sender, EventArgs e)
        {
            finishedGA = false;
            drawMaze();
        }

        private void buttonShowPop_Click(object sender, EventArgs e)
        {
            int[] pathBits = new int[length];
            double pathFit;
            int i;
            List<string> myPop = new List<string>();
            
            // copy the current population to the list box 
            // (decode first)
            for (i = 0; i < populationSize; i++)
            {
                pathBits = myGA.getRoute(i);
                pathFit = myGA.returnFitness(i);
                string path = decodeChromo(pathBits);
                path += "\t";
                path += String.Format("{0:0.00}", pathFit); 
                //add this string to a list of items
                myPop.Add(path);
              
            }
            // create the pop up form
            // pass it a set of text to populate its list box
            // passed a List containing each member of the population
            // and a separate string with the best route
            int best = myGA.getBest();
            int[] route = myGA.getRoute(best);
            string bestRoute = decodeChromo(route);
            Population thePop = new Population(myPop, bestRoute);
            thePop.Show();
            
          }

        private void crossCombo_SelectedIndexChanged(object sender, EventArgs e)
        {

        }

        private void label12_Click(object sender, EventArgs e)
        {

        }

        
     
    }
}
