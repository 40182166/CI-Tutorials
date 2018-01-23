namespace PathFinder
{
    partial class Population
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.listBoxPop = new System.Windows.Forms.ListBox();
            this.closeButton = new System.Windows.Forms.Button();
            this.label1 = new System.Windows.Forms.Label();
            this.bestRouteText = new System.Windows.Forms.TextBox();
            this.SuspendLayout();
            // 
            // listBoxPop
            // 
            this.listBoxPop.FormattingEnabled = true;
            this.listBoxPop.HorizontalScrollbar = true;
            this.listBoxPop.Location = new System.Drawing.Point(21, 73);
            this.listBoxPop.Name = "listBoxPop";
            this.listBoxPop.Size = new System.Drawing.Size(528, 368);
            this.listBoxPop.TabIndex = 0;
            this.listBoxPop.SelectedIndexChanged += new System.EventHandler(this.listBoxPop_SelectedIndexChanged);
            // 
            // closeButton
            // 
            this.closeButton.Location = new System.Drawing.Point(410, 457);
            this.closeButton.Name = "closeButton";
            this.closeButton.Size = new System.Drawing.Size(139, 44);
            this.closeButton.TabIndex = 1;
            this.closeButton.Text = "Close";
            this.closeButton.UseVisualStyleBackColor = true;
            this.closeButton.Click += new System.EventHandler(this.closeButton_Click);
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(22, 24);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(60, 13);
            this.label1.TabIndex = 2;
            this.label1.Text = "Best Route";
            // 
            // bestRouteText
            // 
            this.bestRouteText.Location = new System.Drawing.Point(112, 19);
            this.bestRouteText.Name = "bestRouteText";
            this.bestRouteText.Size = new System.Drawing.Size(436, 20);
            this.bestRouteText.TabIndex = 3;
            // 
            // Population
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(583, 513);
            this.Controls.Add(this.bestRouteText);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.closeButton);
            this.Controls.Add(this.listBoxPop);
            this.Name = "Population";
            this.Text = "Population";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.ListBox listBoxPop;
        private System.Windows.Forms.Button closeButton;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.TextBox bestRouteText;
    }
}