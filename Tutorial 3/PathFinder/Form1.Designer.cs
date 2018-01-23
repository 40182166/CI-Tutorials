namespace PathFinder
{
    partial class Form1
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
            this.goButton = new System.Windows.Forms.Button();
            this.mazePanel = new System.Windows.Forms.Panel();
            this.panel1 = new System.Windows.Forms.Panel();
            this.crossCombo = new System.Windows.Forms.ComboBox();
            this.label10 = new System.Windows.Forms.Label();
            this.label6 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.itersUpDown = new System.Windows.Forms.NumericUpDown();
            this.mutRateText = new System.Windows.Forms.TextBox();
            this.label4 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.maxStepsUpDown = new System.Windows.Forms.NumericUpDown();
            this.label2 = new System.Windows.Forms.Label();
            this.tnSizeUpDown = new System.Windows.Forms.NumericUpDown();
            this.label1 = new System.Windows.Forms.Label();
            this.popSizeUpDown = new System.Windows.Forms.NumericUpDown();
            this.bestLabel = new System.Windows.Forms.Label();
            this.impLabel = new System.Windows.Forms.Label();
            this.bestChromo = new System.Windows.Forms.Label();
            this.resetButton = new System.Windows.Forms.Button();
            this.panel2 = new System.Windows.Forms.Panel();
            this.buttonShowPop = new System.Windows.Forms.Button();
            this.fitnessLbl = new System.Windows.Forms.Label();
            this.iterationsLbl = new System.Windows.Forms.Label();
            this.label9 = new System.Windows.Forms.Label();
            this.label8 = new System.Windows.Forms.Label();
            this.label7 = new System.Windows.Forms.Label();
            this.label11 = new System.Windows.Forms.Label();
            this.label12 = new System.Windows.Forms.Label();
            this.replaceCB = new System.Windows.Forms.ComboBox();
            this.mutateCB = new System.Windows.Forms.ComboBox();
            this.panel1.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.itersUpDown)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.maxStepsUpDown)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.tnSizeUpDown)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.popSizeUpDown)).BeginInit();
            this.panel2.SuspendLayout();
            this.SuspendLayout();
            // 
            // goButton
            // 
            this.goButton.Location = new System.Drawing.Point(378, 337);
            this.goButton.Name = "goButton";
            this.goButton.Size = new System.Drawing.Size(84, 37);
            this.goButton.TabIndex = 0;
            this.goButton.Text = "Go!";
            this.goButton.UseVisualStyleBackColor = true;
            this.goButton.Click += new System.EventHandler(this.goButton_Click);
            // 
            // mazePanel
            // 
            this.mazePanel.Location = new System.Drawing.Point(47, 29);
            this.mazePanel.Name = "mazePanel";
            this.mazePanel.Size = new System.Drawing.Size(414, 295);
            this.mazePanel.TabIndex = 5;
            this.mazePanel.Paint += new System.Windows.Forms.PaintEventHandler(this.mazePanel_Paint);
            // 
            // panel1
            // 
            this.panel1.BackColor = System.Drawing.SystemColors.GradientInactiveCaption;
            this.panel1.Controls.Add(this.mutateCB);
            this.panel1.Controls.Add(this.replaceCB);
            this.panel1.Controls.Add(this.label12);
            this.panel1.Controls.Add(this.label11);
            this.panel1.Controls.Add(this.crossCombo);
            this.panel1.Controls.Add(this.label10);
            this.panel1.Controls.Add(this.label6);
            this.panel1.Controls.Add(this.label5);
            this.panel1.Controls.Add(this.itersUpDown);
            this.panel1.Controls.Add(this.mutRateText);
            this.panel1.Controls.Add(this.label4);
            this.panel1.Controls.Add(this.label3);
            this.panel1.Controls.Add(this.maxStepsUpDown);
            this.panel1.Controls.Add(this.label2);
            this.panel1.Controls.Add(this.tnSizeUpDown);
            this.panel1.Controls.Add(this.label1);
            this.panel1.Controls.Add(this.popSizeUpDown);
            this.panel1.Location = new System.Drawing.Point(19, 380);
            this.panel1.Name = "panel1";
            this.panel1.Size = new System.Drawing.Size(499, 149);
            this.panel1.TabIndex = 6;
            // 
            // crossCombo
            // 
            this.crossCombo.FormattingEnabled = true;
            this.crossCombo.Items.AddRange(new object[] {
            "One point",
            "Two point",
            "Uniform"});
            this.crossCombo.Location = new System.Drawing.Point(383, 29);
            this.crossCombo.Name = "crossCombo";
            this.crossCombo.Size = new System.Drawing.Size(113, 21);
            this.crossCombo.TabIndex = 12;
            this.crossCombo.SelectedIndexChanged += new System.EventHandler(this.crossCombo_SelectedIndexChanged);
            // 
            // label10
            // 
            this.label10.AutoSize = true;
            this.label10.Location = new System.Drawing.Point(323, 31);
            this.label10.Name = "label10";
            this.label10.Size = new System.Drawing.Size(54, 13);
            this.label10.TabIndex = 11;
            this.label10.Text = "Crossover";
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Font = new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label6.Location = new System.Drawing.Point(135, 11);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(174, 13);
            this.label6.TabIndex = 10;
            this.label6.Text = "Genetic Algorithm Parameters";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(17, 109);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(50, 13);
            this.label5.TabIndex = 9;
            this.label5.Text = "Iterations";
            // 
            // itersUpDown
            // 
            this.itersUpDown.Location = new System.Drawing.Point(95, 109);
            this.itersUpDown.Maximum = new decimal(new int[] {
            100000,
            0,
            0,
            0});
            this.itersUpDown.Name = "itersUpDown";
            this.itersUpDown.Size = new System.Drawing.Size(54, 20);
            this.itersUpDown.TabIndex = 8;
            this.itersUpDown.Value = new decimal(new int[] {
            1000,
            0,
            0,
            0});
            // 
            // mutRateText
            // 
            this.mutRateText.Location = new System.Drawing.Point(248, 69);
            this.mutRateText.Name = "mutRateText";
            this.mutRateText.Size = new System.Drawing.Size(48, 20);
            this.mutRateText.TabIndex = 7;
            this.mutRateText.Text = "0.02";
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(165, 76);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(74, 13);
            this.label4.TabIndex = 6;
            this.label4.Text = "Mutation Rate";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(165, 31);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(81, 13);
            this.label3.TabIndex = 5;
            this.label3.Text = "Maximum Steps";
            // 
            // maxStepsUpDown
            // 
            this.maxStepsUpDown.Location = new System.Drawing.Point(248, 27);
            this.maxStepsUpDown.Maximum = new decimal(new int[] {
            50,
            0,
            0,
            0});
            this.maxStepsUpDown.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.maxStepsUpDown.Name = "maxStepsUpDown";
            this.maxStepsUpDown.Size = new System.Drawing.Size(51, 20);
            this.maxStepsUpDown.TabIndex = 4;
            this.maxStepsUpDown.Value = new decimal(new int[] {
            35,
            0,
            0,
            0});
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(11, 71);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(87, 13);
            this.label2.TabIndex = 3;
            this.label2.Text = "Tournament Size";
            // 
            // tnSizeUpDown
            // 
            this.tnSizeUpDown.Location = new System.Drawing.Point(104, 69);
            this.tnSizeUpDown.Maximum = new decimal(new int[] {
            500,
            0,
            0,
            0});
            this.tnSizeUpDown.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.tnSizeUpDown.Name = "tnSizeUpDown";
            this.tnSizeUpDown.Size = new System.Drawing.Size(46, 20);
            this.tnSizeUpDown.TabIndex = 2;
            this.tnSizeUpDown.Value = new decimal(new int[] {
            2,
            0,
            0,
            0});
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(9, 31);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(80, 13);
            this.label1.TabIndex = 1;
            this.label1.Text = "Population Size";
            // 
            // popSizeUpDown
            // 
            this.popSizeUpDown.Location = new System.Drawing.Point(104, 29);
            this.popSizeUpDown.Maximum = new decimal(new int[] {
            500,
            0,
            0,
            0});
            this.popSizeUpDown.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.popSizeUpDown.Name = "popSizeUpDown";
            this.popSizeUpDown.Size = new System.Drawing.Size(46, 20);
            this.popSizeUpDown.TabIndex = 0;
            this.popSizeUpDown.Value = new decimal(new int[] {
            50,
            0,
            0,
            0});
            // 
            // bestLabel
            // 
            this.bestLabel.AutoSize = true;
            this.bestLabel.Location = new System.Drawing.Point(87, 358);
            this.bestLabel.Name = "bestLabel";
            this.bestLabel.Size = new System.Drawing.Size(0, 13);
            this.bestLabel.TabIndex = 7;
            // 
            // impLabel
            // 
            this.impLabel.AutoSize = true;
            this.impLabel.Location = new System.Drawing.Point(286, 354);
            this.impLabel.Name = "impLabel";
            this.impLabel.Size = new System.Drawing.Size(0, 13);
            this.impLabel.TabIndex = 9;
            // 
            // bestChromo
            // 
            this.bestChromo.AutoSize = true;
            this.bestChromo.Location = new System.Drawing.Point(71, 547);
            this.bestChromo.Name = "bestChromo";
            this.bestChromo.Size = new System.Drawing.Size(0, 13);
            this.bestChromo.TabIndex = 10;
            // 
            // resetButton
            // 
            this.resetButton.Location = new System.Drawing.Point(304, 337);
            this.resetButton.Name = "resetButton";
            this.resetButton.Size = new System.Drawing.Size(68, 37);
            this.resetButton.TabIndex = 11;
            this.resetButton.Text = "Clear";
            this.resetButton.UseVisualStyleBackColor = true;
            this.resetButton.Click += new System.EventHandler(this.resetButton_Click);
            // 
            // panel2
            // 
            this.panel2.BackColor = System.Drawing.Color.FromArgb(((int)(((byte)(128)))), ((int)(((byte)(255)))), ((int)(((byte)(255)))));
            this.panel2.Controls.Add(this.buttonShowPop);
            this.panel2.Controls.Add(this.fitnessLbl);
            this.panel2.Controls.Add(this.iterationsLbl);
            this.panel2.Controls.Add(this.label9);
            this.panel2.Controls.Add(this.label8);
            this.panel2.Controls.Add(this.label7);
            this.panel2.Location = new System.Drawing.Point(19, 535);
            this.panel2.Name = "panel2";
            this.panel2.Size = new System.Drawing.Size(499, 60);
            this.panel2.TabIndex = 12;
            // 
            // buttonShowPop
            // 
            this.buttonShowPop.Location = new System.Drawing.Point(396, 6);
            this.buttonShowPop.Name = "buttonShowPop";
            this.buttonShowPop.Size = new System.Drawing.Size(85, 45);
            this.buttonShowPop.TabIndex = 5;
            this.buttonShowPop.Text = "View population";
            this.buttonShowPop.UseVisualStyleBackColor = true;
            this.buttonShowPop.Click += new System.EventHandler(this.buttonShowPop_Click);
            // 
            // fitnessLbl
            // 
            this.fitnessLbl.AutoSize = true;
            this.fitnessLbl.Location = new System.Drawing.Point(111, 23);
            this.fitnessLbl.Name = "fitnessLbl";
            this.fitnessLbl.Size = new System.Drawing.Size(0, 13);
            this.fitnessLbl.TabIndex = 4;
            // 
            // iterationsLbl
            // 
            this.iterationsLbl.AutoSize = true;
            this.iterationsLbl.Location = new System.Drawing.Point(305, 22);
            this.iterationsLbl.Name = "iterationsLbl";
            this.iterationsLbl.Size = new System.Drawing.Size(0, 13);
            this.iterationsLbl.TabIndex = 3;
            // 
            // label9
            // 
            this.label9.AutoSize = true;
            this.label9.Font = new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label9.Location = new System.Drawing.Point(213, 23);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(60, 13);
            this.label9.TabIndex = 2;
            this.label9.Text = "Iterations";
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(101, 25);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(0, 13);
            this.label8.TabIndex = 1;
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Font = new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label7.Location = new System.Drawing.Point(13, 25);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(80, 13);
            this.label7.TabIndex = 0;
            this.label7.Text = "Best Fitness:";
            // 
            // label11
            // 
            this.label11.AutoSize = true;
            this.label11.Location = new System.Drawing.Point(165, 116);
            this.label11.Name = "label11";
            this.label11.Size = new System.Drawing.Size(48, 13);
            this.label11.TabIndex = 13;
            this.label11.Text = "Mutation";
            // 
            // label12
            // 
            this.label12.AutoSize = true;
            this.label12.Location = new System.Drawing.Point(307, 76);
            this.label12.Name = "label12";
            this.label12.Size = new System.Drawing.Size(70, 13);
            this.label12.TabIndex = 14;
            this.label12.Text = "Replacement";
            this.label12.Click += new System.EventHandler(this.label12_Click);
            // 
            // replaceCB
            // 
            this.replaceCB.FormattingEnabled = true;
            this.replaceCB.Items.AddRange(new object[] {
            "replace worst",
            "myReplace"});
            this.replaceCB.Location = new System.Drawing.Point(385, 78);
            this.replaceCB.Name = "replaceCB";
            this.replaceCB.Size = new System.Drawing.Size(110, 21);
            this.replaceCB.TabIndex = 15;
            // 
            // mutateCB
            // 
            this.mutateCB.FormattingEnabled = true;
            this.mutateCB.Items.AddRange(new object[] {
            "Standard",
            "myMutate"});
            this.mutateCB.Location = new System.Drawing.Point(247, 118);
            this.mutateCB.Name = "mutateCB";
            this.mutateCB.Size = new System.Drawing.Size(90, 21);
            this.mutateCB.TabIndex = 16;
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(530, 625);
            this.Controls.Add(this.panel2);
            this.Controls.Add(this.resetButton);
            this.Controls.Add(this.bestChromo);
            this.Controls.Add(this.impLabel);
            this.Controls.Add(this.bestLabel);
            this.Controls.Add(this.panel1);
            this.Controls.Add(this.mazePanel);
            this.Controls.Add(this.goButton);
            this.Name = "Form1";
            this.Text = "PathFinder GA";
            this.panel1.ResumeLayout(false);
            this.panel1.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.itersUpDown)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.maxStepsUpDown)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.tnSizeUpDown)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.popSizeUpDown)).EndInit();
            this.panel2.ResumeLayout(false);
            this.panel2.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button goButton;
        private System.Windows.Forms.Panel mazePanel;
        private System.Windows.Forms.Panel panel1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.NumericUpDown tnSizeUpDown;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.NumericUpDown popSizeUpDown;
        private System.Windows.Forms.TextBox mutRateText;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.NumericUpDown maxStepsUpDown;
        private System.Windows.Forms.NumericUpDown itersUpDown;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.Label bestLabel;
        private System.Windows.Forms.Label impLabel;
        private System.Windows.Forms.Label bestChromo;
        private System.Windows.Forms.Button resetButton;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.Panel panel2;
        private System.Windows.Forms.Label fitnessLbl;
        private System.Windows.Forms.Label iterationsLbl;
        private System.Windows.Forms.Label label9;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.Label label10;
        private System.Windows.Forms.ComboBox crossCombo;
        private System.Windows.Forms.Button buttonShowPop;
        private System.Windows.Forms.Label label11;
        private System.Windows.Forms.Label label12;
        private System.Windows.Forms.ComboBox mutateCB;
        private System.Windows.Forms.ComboBox replaceCB;
    }
}

