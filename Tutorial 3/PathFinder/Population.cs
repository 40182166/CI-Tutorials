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
    public partial class Population : Form
    {
        public Population(List<string> pop, string bestRoute)
        {
          
            string path;
            InitializeComponent();
            
            listBoxPop.Font = new System.Drawing.Font("Courier New", (float)10.0, System.Drawing.FontStyle.Regular);
            bestRouteText.Text = bestRoute;

            foreach (string text in pop)
            {
                path = text;
                listBoxPop.Items.Add(path);
            }
        }

        private void listBoxPop_SelectedIndexChanged(object sender, EventArgs e)
        {

        }

        private void closeButton_Click(object sender, EventArgs e)
        {
            this.Close();
        }
    }
}
