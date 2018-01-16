using System;
using System.Collections;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;


namespace Goofspiel
{
    
    public partial class Goofspiel : Form
    {

        Random computerRandom;
        ArrayList Deck;
        ArrayList permutation;
        int turn;
        int[] score = new int[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 };
        int upturnedCardValue;
        int computerPlayerValue;
        int humanPlayerValue;
        int humanScore;
        int computerScore;
        int humanWins = 0;
        int computerWins = 0;
        int cardsTurned = 0;
        int strategyUsed;
        Boolean[] humanPlayedCards;
        Boolean[] playedCards; 

        /**********************************************************************************
        * HOW TO MODIFY THE GAME TO UPDATE THE COMPUTER'S STRATEGY 
        *   In order to modify the game, you should only need to update the method below, 
        *   computerPlayerMove(); The comments in the code should give sufficient direction.
        *
        * VARIABLES YOU MIGHT NEED TO CHANGE THE STRATEGY
        * There are two pieces of information you might wish to use in designing a new strategy:
        * 
        * a) the list of cards already played by the other player:
        *      this is held in the array Boolean humanPlayedCards[]
        *      this INDEX of the array specifies the card (0=Ace, ..., 12=King)
        *      Each item in the array has a boolean value true/false depending on 
        *      whether the card has been played. 
        * 
        *  b) The value of the current upturned card:
        *       this is stored in a variable: int upturnedCardValue
        *       (0=ace, .... 13=king)
         * * * * */


        private void computerPlayerMove()
        {

            // picks a card for the computer player
            // the computer player can be given any strategy

            // STRATEGY 1: play a random card (that hasn't already been played)


            Boolean cardPlayed = true;
            int chosen = 0;  // this variable holds in the index of the card chosen

            /********************************************************************************
             * make your modifications here
             * ******************************************************************************/

            // This is the code to randomly choose a card
            // replace this if you want to change the strategy
            // make sure that the variable CHOSEN records the index of the card played,
            // however you decide to choose that card

            // The loop stops when the random card selected is one that has not already been chosen
            // (its not very efficient but speed doesn't really matter here

            //----------------------------------STRATEGIES-----------------------------------//
            switch (strategyUsed)
            {
                case 1:
                    while (cardPlayed)
                    {
                        chosen = computerRandom.Next(13);   // pick a random card
                        cardPlayed = playedCards[chosen];
                    }
                    break;
                case 2:
                    chosen = upturnedCardValue - 1;
                    cardPlayed = playedCards[chosen];
                    break;
                case 3:
                    chosen = upturnedCardValue - 1;
                    cardPlayed = playedCards[chosen];
                    break;
            }





            /******************************************************************************
             * *  END modifications
             * ****************************************************************************/

            // this code updates the array with the card chosen by the computer so it can't be chosen again, and
            // then  copies the chosen image to the computer player picture-box. It also records the value of the
            // card played so the updateScore method can determine the winning bidder
            // You will need this code whatever strategy you use

            playedCards[chosen] = true;
            computerPlayerValue = score[chosen];
            computerPB.Image = clubs.Images[chosen];

           


        }



        /************************************************************************
         *  None of the code below needs to be updated in order to modify the
         *  computer player's strategy
         * 
         * 
         * **********************************************************************/

        public Goofspiel()
        {
            InitializeComponent();
            permutation = new ArrayList();
            computerRandom = new Random();
            Deck = new ArrayList();
            disableCards();
            nextBtn.Enabled = false;
            humanPlayedCards = new Boolean[13];
            playedCards = new Boolean[13];     // keeps track of cards played this round
            
            // add images to players hand
            pb1.Image = spades.Images[0];
            pb2.Image = spades.Images[1];
            pb3.Image = spades.Images[2];
            pb4.Image = spades.Images[3];
            pb5.Image = spades.Images[4];
            pb6.Image = spades.Images[5];
            pb7.Image = spades.Images[6];
            pb8.Image = spades.Images[7];
            pb9.Image = spades.Images[8];
            pb10.Image = spades.Images[9];
            pb11.Image = spades.Images[10];
            pb12.Image = spades.Images[11];
            pb13.Image = spades.Images[12];

        }

        private void updateScores()
        {
            if (humanPlayerValue > computerPlayerValue)
                humanScore += upturnedCardValue;
            else if (humanPlayerValue < computerPlayerValue)
                computerScore += upturnedCardValue;

            compScoreTxt.Text = computerScore.ToString();
            humanScoreTxt.Text = humanScore.ToString();
        }

        private void getPermutation()
        {
            // returns a permutation of the numbers 0-12
            int size;
           
            
            // delete any current permutation
            permutation.Clear();
            
            
            // create  a list with 10 numbers in it
            ArrayList myArray = new ArrayList();
            for (int i = 0; i < 13; i++)
                myArray.Add(i);

            for (int i=0;i<13;i++){
             size = myArray.Count;
             int key = computerRandom.Next(size);
             permutation.Add(myArray[key]);
             myArray.RemoveAt(key);
            }
                
            
        }

        private void startButton_Click(object sender, EventArgs e)
        {
            //initializing random strategy to be used by computer
            Random rnd = new Random();
            strategyUsed = rnd.Next(1, 4); // creates a number between 1 and 3

            // create a random order to the cards to be revealed
            getPermutation();

            Deck.Clear();

            // use this permutation to create a new shuffled deck
            for (int i = 0; i < 13; i++)
                Deck.Add(permutation[i]);


            // initialise turn to 0
            turn = 0;
           

            // initialise scores
            humanScore = computerScore = 0;

            // reset players cards
            resetCards();

            for (int i = 0; i < 13; i++)
            {
                humanPlayedCards[i] = false;
                playedCards[i] = false;
            }

            // reset scores
            humanScoreTxt.Text = "";
            compScoreTxt.Text = "";

            // and turn card
            int index = (int)Deck[turn];
            deckPB.Image = hearts.Images[index];
            cardsTurned = 0;
            
            // keep record of value of upturned card
            upturnedCardValue = (int)Deck[turn] + 1;

            nextBtn.Enabled = false;
            enableCards();
            
        }

        private void nextBtn_Click(object sender, EventArgs e)
        {
            // turn over  a new card#
            turn++;
            int index = (int)Deck[turn];
            deckPB.Image = hearts.Images[index];
            upturnedCardValue = (int)Deck[turn]+1;
          
            // and blank out the previous cards

            computerPB.Image = cardBacks.Images[0];
            humanPB.Image = cardBacks.Images[0];
            
            // now disable the button until someone has played
            nextBtn.Enabled = false;
            nextBtn.BackColor = Color.Gray;
            enableCards();
        }

      
        // these methods deal with clicking on each of the cards in the human players hand

        private void pb1_Click(object sender, EventArgs e)
        {

            // copy to played PB
            humanPB.Image = pb1.Image;
            humanPlayerValue = 1;
            humanPlayedCards[0] = true;

            // turn over computer player
            computerPlayerMove();

            // calculate score
            updateScores();

            // remove card from pile
            pb1.Hide();
            nextBtn.Enabled = true;
            disableCards();
            cardsTurned++;
            checkFinishedGame();
        }
        private void pb2_Click(object sender, EventArgs e)
        {
            // copy to played PB
            humanPB.Image = pb2.Image;
            humanPlayerValue = 2;
            humanPlayedCards[1] = true;

            // turn over computer player
            computerPlayerMove();

            // calculate score
            updateScores();

            // remove card from pile
            pb2.Hide();
            nextBtn.Enabled = true;
            disableCards();
            cardsTurned++;
            checkFinishedGame();
        }

       
        private void pb3_Click(object sender, EventArgs e)
        {
            // copy to played PB
            humanPB.Image = pb3.Image;
            humanPlayerValue = 3;
            humanPlayedCards[2] = true;

            // turn over computer player
            computerPlayerMove();

            // calculate score
            updateScores();

            // remove card from pile
            pb3.Hide();
            nextBtn.Enabled = true;
            disableCards();
            cardsTurned++;
            checkFinishedGame();
        }

        private void pb4_Click(object sender, EventArgs e)
        {
            // copy to played PB
            humanPB.Image = pb4.Image;
            humanPlayerValue = 4;
            humanPlayedCards[3] = true;

            // turn over computer player
            computerPlayerMove();

            // calculate score
            updateScores();

            // remove card from pile
            pb4.Hide();
            nextBtn.Enabled = true;
            disableCards();
            cardsTurned++;
            checkFinishedGame();
        }

        private void pb5_Click(object sender, EventArgs e)
        {
            // copy to played PB
            humanPB.Image = pb5.Image;
            humanPlayerValue = 5;
            humanPlayedCards[4] = true;

            // turn over computer player
            computerPlayerMove();

            // calculate score
            updateScores();

            // remove card from pile
            pb5.Hide();
            nextBtn.Enabled = true;
            disableCards();
            cardsTurned++;
            checkFinishedGame();
        }

        private void pb6_Click(object sender, EventArgs e)
        {
            // copy to played PB
            humanPB.Image = pb6.Image;
            humanPlayerValue = 6;
            humanPlayedCards[5] = true;

            // turn over computer player
            computerPlayerMove();

            // calculate score
            updateScores();

            // remove card from pile
            pb6.Hide();
            nextBtn.Enabled = true;
            disableCards();
            cardsTurned++;
            checkFinishedGame();
        }

        private void pb7_Click(object sender, EventArgs e)
        {
            // copy to played PB
            humanPB.Image = pb7.Image;
            humanPlayerValue = 7;
            humanPlayedCards[6] = true;

            // turn over computer player
            computerPlayerMove();

            // calculate score
            updateScores();

            // remove card from pile
            pb7.Hide();
            nextBtn.Enabled = true;
            disableCards();
            cardsTurned++;
            checkFinishedGame();
        }

        private void pb8_Click(object sender, EventArgs e)
        {
            // copy to played PB
            humanPB.Image = pb8.Image;
            humanPlayerValue = 8;
            humanPlayedCards[7] = true;

            // turn over computer player
            computerPlayerMove();

            // calculate score
            updateScores();

            // remove card from pile
            pb8.Hide();
            nextBtn.Enabled = true;
            disableCards();
            cardsTurned++;
            checkFinishedGame();
        }

        private void pb9_Click(object sender, EventArgs e)
        {
            // copy to played PB
            humanPB.Image = pb9.Image;
            humanPlayerValue = 9;
            humanPlayedCards[8] = true;
            
            // turn over computer player
            computerPlayerMove();

            // calculate score
            updateScores();

            // remove card from pile
            pb9.Hide();
            nextBtn.Enabled = true;
            disableCards();
            cardsTurned++;
            checkFinishedGame();
        }

        private void pb10_Click(object sender, EventArgs e)
        {
            // copy to played PB
            humanPB.Image = pb10.Image;
            humanPlayerValue = 10;
            humanPlayedCards[9] = true;

            // turn over computer player
            computerPlayerMove();

            // calculate score
            updateScores();

            // remove card from pile
            pb10.Hide();
            nextBtn.Enabled = true;
            disableCards();
            cardsTurned++;
            checkFinishedGame();
        }

        private void pb11_Click(object sender, EventArgs e)
        {
            // copy to played PB
            humanPB.Image = pb11.Image;
            humanPlayerValue = 11;
            humanPlayedCards[10] = true;

            // turn over computer player
            computerPlayerMove();

            // calculate score
            updateScores();

            // remove card from pile
            pb11.Hide();
            nextBtn.Enabled = true;
            disableCards();
            cardsTurned++;
            checkFinishedGame();
        }

        private void pb12_Click(object sender, EventArgs e)
        {
            // copy to played PB
            humanPB.Image = pb12.Image;
            humanPlayerValue = 12;
            humanPlayedCards[11] = true;

            // turn over computer player
            computerPlayerMove();

            // calculate score
            updateScores();

            // remove card from pile
            pb12.Hide();
            nextBtn.Enabled = true;
            disableCards();
            cardsTurned++;
            checkFinishedGame();
        }

        private void pb13_Click(object sender, EventArgs e)
        {
            // copy to played PB
            humanPB.Image = pb13.Image;
            humanPlayerValue = 13;
            humanPlayedCards[12] = true;

            // turn over computer player
            computerPlayerMove();

            // calculate score
            updateScores();

            // remove card from pile
            pb13.Hide();
            nextBtn.Enabled = true;
            disableCards();
            cardsTurned++;
            checkFinishedGame();
        }

        // ensures that only one card can be played each round
        // (e.g can't play two cards on top of one upturned card

        private void disableCards() 
        {
            pb1.Enabled = false;
            pb2.Enabled = false;
            pb3.Enabled = false;
            pb4.Enabled = false;
            pb5.Enabled = false;
            pb6.Enabled = false;
            pb7.Enabled = false;
            pb8.Enabled = false;
            pb9.Enabled = false;
            pb10.Enabled = false;
            pb11.Enabled = false;
            pb12.Enabled = false;
            pb13.Enabled = false;
            nextBtn.BackColor = Color.Red;

           
        }

        // enables cards again

        private void enableCards()
        {
            pb1.Enabled = true;
            pb2.Enabled = true;
            pb3.Enabled = true;
            pb4.Enabled = true;
            pb5.Enabled = true;
            pb6.Enabled = true;
            pb7.Enabled = true;
            pb8.Enabled = true;
            pb9.Enabled = true;
            pb10.Enabled = true;
            pb11.Enabled = true;
            pb12.Enabled = true;
            pb13.Enabled = true;
            nextBtn.BackColor = Color.Gray;

        }

        private void resetCards()
        {
            // unhide all the players cards
            pb1.Show();
            pb2.Show();
            pb3.Show();
            pb4.Show();
            pb5.Show();
            pb6.Show();
            pb7.Show();
            pb8.Show();
            pb9.Show();
            pb10.Show();
            pb11.Show();
            pb12.Show();
            pb13.Show();


            // and blank out the previous cards

            computerPB.Image = cardBacks.Images[0];
            humanPB.Image = cardBacks.Images[0];


        }

        private void checkFinishedGame()
        {
            string messageTxt;

            // called when all cards have been turned
            if (cardsTurned == 13)
            {
                // display a popup saying who won
                if (humanScore > computerScore)
                {
                    messageTxt = "Well done! you won! \n Press START for a new game";
                    humanWins++;
                    hpRounds.Text = humanWins.ToString();
                }
                else if (humanScore < computerScore)
                {
                    messageTxt = "Hard luck, the computer won\n Press START for a new game";
                    computerWins++;
                    cpRounds.Text = computerWins.ToString();
                }
                else
                    messageTxt = "Game was drawn\n Press START for a new game";

                MessageBox.Show(messageTxt);

                // call start
            }
        }

        private void helpBtn_Click(object sender, EventArgs e)
        {
            string messageTxt = "HOW TO PLAY\n Goofspiel is played using cards from a standard deck of cards, and is typically a two-player game. Each suit is ranked A (low), 2, ..., 10, J, Q, K (high).\n\n 1. One suit is singled out as \"competition\" suit (in this version, we use the HEARTS suit); each of the remaining suits becomes a hand for one player. The hearts are shuffled and placed between the players with one card turned up.\n\n 2. Play proceeds in a series of rounds. The players make \"closed bids\" for the top (face up) heart by selecting a card from their hand (keeping their choice secret from their opponent). Once these cards are selected, they are simultaneously revealed, and the player making the highest bid takes the competition card. In the case of a tie, the competition card is discarded and neither player scores (other variants of this rule are sometimes played, however are not included in this game).\n\n 3. The cards used for bidding are discarded, and play continues with a new upturned heart card.\n\n 4. After 13 rounds, there are no remaining cards and the game ends. Typically, players earn points equal to sum of the ranks of cards won (i.e. Ace is worth one point, 2 is two points, etc, Jack being worth 11, Queen 12, and King worth 13 points). \n\n\n In tnis implementation, the computer (clubs) plays with a fixed strategy coded into the program. The computer plays against a human (spades) who selects cards to play using the mouse\n\n\nINSTRUCTIONS\n\n Press START to start a new game\n\n The computer automatically plays its card, the human player needs to click on a card to select the one they wish to play \n\n Click NEXT CARD to turn the next card in the deck for the next round \n\n Scoring is updated automatically";
            MessageBox.Show(messageTxt);
        }

        private void deckPB_Click(object sender, EventArgs e)
        {

        }

    }
}
