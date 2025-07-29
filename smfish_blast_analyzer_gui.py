# smfish_blast_analyzer_gui.py

import pandas as pd
import re
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os


class SmFISHBlastAnalyzerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("smFISH BLAST Results Analyzer")
        self.root.geometry("800x600")
        self.root.resizable(True, True)

        # Variables
        self.blast_file = tk.StringVar()
        self.csv_files = []
        self.probe_names_column = tk.StringVar()
        self.sequence_column = tk.StringVar()
        self.output_dir = tk.StringVar()
        self.combined_df = None
        self.available_columns = []

        self.create_widgets()

    def create_widgets(self):
        # Create notebook for tabs
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Tab 1: File Selection
        file_tab = ttk.Frame(notebook)
        notebook.add(file_tab, text="1. File Selection")

        # Tab 2: Column Configuration
        config_tab = ttk.Frame(notebook)
        notebook.add(config_tab, text="2. Column Configuration")

        # Tab 3: Analysis & Output
        analysis_tab = ttk.Frame(notebook)
        notebook.add(analysis_tab, text="3. Analysis & Output")

        self.create_file_tab(file_tab)
        self.create_config_tab(config_tab)
        self.create_analysis_tab(analysis_tab)

    def create_file_tab(self, parent):
        main_frame = ttk.Frame(parent, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # BLAST file selection
        ttk.Label(
            main_frame, text="1. Select BLAST Results File:", font=("", 10, "bold")
        ).pack(anchor=tk.W, pady=(0, 5))

        blast_frame = ttk.Frame(main_frame)
        blast_frame.pack(fill=tk.X, pady=(0, 15))

        self.blast_entry = ttk.Entry(
            blast_frame, textvariable=self.blast_file, state="readonly"
        )
        self.blast_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 5))

        ttk.Button(blast_frame, text="Browse", command=self.browse_blast_file).pack(
            side=tk.RIGHT
        )

        # CSV files selection
        ttk.Label(
            main_frame,
            text="2. Select CSV Files with Probe Data:",
            font=("", 10, "bold"),
        ).pack(anchor=tk.W, pady=(0, 5))

        csv_frame = ttk.Frame(main_frame)
        csv_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        # CSV files listbox with scrollbar
        list_frame = ttk.Frame(csv_frame)
        list_frame.pack(fill=tk.BOTH, expand=True)

        self.csv_listbox = tk.Listbox(list_frame, height=8)
        scrollbar = ttk.Scrollbar(
            list_frame, orient="vertical", command=self.csv_listbox.yview
        )
        self.csv_listbox.configure(yscrollcommand=scrollbar.set)

        self.csv_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # CSV control buttons
        csv_button_frame = ttk.Frame(csv_frame)
        csv_button_frame.pack(fill=tk.X, pady=(5, 0))

        ttk.Button(
            csv_button_frame, text="Add CSV Files", command=self.add_csv_files
        ).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(
            csv_button_frame, text="Remove Selected", command=self.remove_csv_file
        ).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(
            csv_button_frame, text="Clear All", command=self.clear_csv_files
        ).pack(side=tk.LEFT)

        # Load and combine button
        ttk.Button(
            main_frame,
            text="Load and Combine Files",
            command=self.load_and_combine_files,
            style="Accent.TButton",
        ).pack(pady=20)

    def create_config_tab(self, parent):
        main_frame = ttk.Frame(parent, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Instructions
        instructions = ttk.Label(
            main_frame,
            text="Configure column mappings after loading files:",
            font=("", 10, "bold"),
        )
        instructions.pack(anchor=tk.W, pady=(0, 15))

        # Column selection frame
        config_frame = ttk.Frame(main_frame)
        config_frame.pack(fill=tk.X, pady=(0, 20))
        config_frame.columnconfigure(1, weight=1)

        # Probe names column
        ttk.Label(config_frame, text="Probe Names Column:").grid(
            row=0, column=0, sticky=tk.W, padx=(0, 10), pady=5
        )
        self.probe_combo = ttk.Combobox(
            config_frame, textvariable=self.probe_names_column, state="readonly"
        )
        self.probe_combo.grid(row=0, column=1, sticky=(tk.W, tk.E), pady=5)

        # Sequence column (the one used for BLAST)
        ttk.Label(config_frame, text="Sequence Column (used for BLAST):").grid(
            row=1, column=0, sticky=tk.W, padx=(0, 10), pady=5
        )
        self.sequence_combo = ttk.Combobox(
            config_frame, textvariable=self.sequence_column, state="readonly"
        )
        self.sequence_combo.grid(row=1, column=1, sticky=(tk.W, tk.E), pady=5)

        # Help text
        help_text = tk.Text(main_frame, height=8, wrap=tk.WORD, bg=self.root.cget("bg"))
        help_text.pack(fill=tk.BOTH, expand=True, pady=(10, 0))
        help_text.insert(
            tk.END,
            "Instructions:\n\n"
            "• Probe Names Column: Select the column containing probe identifiers that match the BLAST query names\n\n"
            "• Sequence Column: Select the column containing the actual sequences that were used for BLAST analysis\n\n"
            "Note: The sequence column is used for validation - it should contain the sequences that were actually submitted to BLAST to ensure data integrity.",
        )
        help_text.config(state=tk.DISABLED)

        # Data preview
        ttk.Label(main_frame, text="Data Preview:", font=("", 10, "bold")).pack(
            anchor=tk.W, pady=(15, 5)
        )

        # Preview frame with scrollbars
        preview_frame = ttk.Frame(main_frame)
        preview_frame.pack(fill=tk.BOTH, expand=True)

        # Treeview for data preview
        self.preview_tree = ttk.Treeview(preview_frame)

        # Scrollbars for preview
        v_scrollbar = ttk.Scrollbar(
            preview_frame, orient="vertical", command=self.preview_tree.yview
        )
        h_scrollbar = ttk.Scrollbar(
            preview_frame, orient="horizontal", command=self.preview_tree.xview
        )
        self.preview_tree.configure(
            yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set
        )

        self.preview_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        v_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        h_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)

    def create_analysis_tab(self, parent):
        main_frame = ttk.Frame(parent, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Output directory selection
        ttk.Label(main_frame, text="Output Directory:", font=("", 10, "bold")).pack(
            anchor=tk.W, pady=(0, 5)
        )

        output_frame = ttk.Frame(main_frame)
        output_frame.pack(fill=tk.X, pady=(0, 15))

        self.output_entry = ttk.Entry(
            output_frame, textvariable=self.output_dir, state="readonly"
        )
        self.output_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 5))

        ttk.Button(output_frame, text="Browse", command=self.browse_output_dir).pack(
            side=tk.RIGHT
        )

        # Analysis options
        ttk.Label(main_frame, text="Analysis Options:", font=("", 10, "bold")).pack(
            anchor=tk.W, pady=(0, 5)
        )

        options_frame = ttk.Frame(main_frame)
        options_frame.pack(fill=tk.X, pady=(0, 15))

        self.unique_hits_only = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            options_frame,
            text="Filter for unique hits only (NumberOfHits == 1)",
            variable=self.unique_hits_only,
        ).pack(anchor=tk.W)

        # Run analysis button
        self.analyze_btn = ttk.Button(
            main_frame,
            text="Run Analysis",
            command=self.run_analysis,
            state="disabled",
            style="Accent.TButton",
        )
        self.analyze_btn.pack(pady=20)

        # Progress bar
        self.progress = ttk.Progressbar(main_frame, mode="indeterminate")
        self.progress.pack(fill=tk.X, pady=(0, 10))

        # Status log
        ttk.Label(main_frame, text="Analysis Log:", font=("", 10, "bold")).pack(
            anchor=tk.W, pady=(10, 5)
        )

        log_frame = ttk.Frame(main_frame)
        log_frame.pack(fill=tk.BOTH, expand=True)

        self.log_text = tk.Text(log_frame, height=12, wrap=tk.WORD)
        log_scrollbar = ttk.Scrollbar(
            log_frame, orient="vertical", command=self.log_text.yview
        )
        self.log_text.configure(yscrollcommand=log_scrollbar.set)

        self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        log_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

    def log_message(self, message):
        """Add message to log"""
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()

    def browse_blast_file(self):
        """Browse for BLAST results file"""
        filename = filedialog.askopenfilename(
            title="Select BLAST results text file",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
        )
        if filename:
            self.blast_file.set(filename)
            self.log_message(f"Selected BLAST file: {os.path.basename(filename)}")

    def add_csv_files(self):
        """Add CSV files to the list"""
        filenames = filedialog.askopenfilenames(
            title="Select CSV files",
            filetypes=[
                ("CSV files", "*.csv"),
                ("Excel files", "*.xlsx"),
                ("All files", "*.*"),
            ],
        )
        for filename in filenames:
            if filename not in self.csv_files:
                self.csv_files.append(filename)
                self.csv_listbox.insert(tk.END, os.path.basename(filename))
                self.log_message(f"Added CSV file: {os.path.basename(filename)}")

    def remove_csv_file(self):
        """Remove selected CSV file"""
        selection = self.csv_listbox.curselection()
        if selection:
            index = selection[0]
            filename = self.csv_files.pop(index)
            self.csv_listbox.delete(index)
            self.log_message(f"Removed CSV file: {os.path.basename(filename)}")

    def clear_csv_files(self):
        """Clear all CSV files"""
        self.csv_files.clear()
        self.csv_listbox.delete(0, tk.END)
        self.log_message("Cleared all CSV files")

    def load_and_combine_files(self):
        """Load and combine all CSV files"""
        if not self.csv_files:
            messagebox.showerror("Error", "Please select at least one CSV file")
            return

        try:
            self.log_message("Loading and combining CSV files...")
            combined_dfs = []

            for i, csv_file in enumerate(self.csv_files):
                self.log_message(
                    f"Loading file {i+1}/{len(self.csv_files)}: {os.path.basename(csv_file)}"
                )

                if csv_file.endswith(".xlsx"):
                    df = pd.read_excel(csv_file)
                else:
                    df = pd.read_csv(csv_file)

                df["SourceFile"] = os.path.basename(csv_file)  # Track source file
                combined_dfs.append(df)
                self.log_message(
                    f"  - Loaded {len(df)} rows, {len(df.columns)} columns"
                )

            # Combine all dataframes
            self.combined_df = pd.concat(combined_dfs, ignore_index=True)
            self.available_columns = [
                col for col in self.combined_df.columns if col != "SourceFile"
            ]

            # Update column dropdowns
            self.probe_combo["values"] = self.available_columns
            self.sequence_combo["values"] = self.available_columns

            # Auto-select common column names
            for col in self.available_columns:
                if any(keyword in col.lower() for keyword in ["probe", "name", "id"]):
                    self.probe_names_column.set(col)
                    break

            for col in self.available_columns:
                if any(
                    keyword in col.lower()
                    for keyword in ["seq", "sequence", "dna", "rna"]
                ):
                    self.sequence_column.set(col)
                    break

            # Update preview
            self.update_preview()

            self.log_message(f"Successfully combined {len(self.csv_files)} files:")
            self.log_message(f"  - Total rows: {len(self.combined_df)}")
            self.log_message(f"  - Columns: {', '.join(self.available_columns)}")

            # Set default output directory
            if self.csv_files:
                default_output = os.path.dirname(self.csv_files[0])
                self.output_dir.set(default_output)

            messagebox.showinfo(
                "Success",
                f"Combined {len(self.csv_files)} files successfully!\n"
                f"Total rows: {len(self.combined_df)}",
            )

        except Exception as e:
            self.log_message(f"Error loading files: {str(e)}")
            messagebox.showerror("Error", f"Failed to load files: {str(e)}")

        # ADDED - Check if ready after loading files
        self.check_ready_for_analysis()

    def update_preview(self):
        """Update the data preview"""
        if self.combined_df is None:
            return

        # Clear existing data
        for item in self.preview_tree.get_children():
            self.preview_tree.delete(item)

        # Set up columns (show first 10 columns to avoid overcrowding)
        display_columns = self.available_columns[:10]
        self.preview_tree["columns"] = display_columns
        self.preview_tree["show"] = "headings"

        # Configure column headings and widths
        for col in display_columns:
            self.preview_tree.heading(col, text=col)
            self.preview_tree.column(col, width=100, minwidth=80)

        # Add sample data (first 20 rows)
        for idx, row in self.combined_df.head(20).iterrows():
            values = [str(row[col]) if col in row else "" for col in display_columns]
            self.preview_tree.insert("", "end", values=values)

    def browse_output_dir(self):
        """Browse for output directory"""
        directory = filedialog.askdirectory(title="Select output directory")
        if directory:
            self.output_dir.set(directory)
            self.log_message(f"Output directory set to: {directory}")
            # ADDED - Check if ready after setting output directory
            self.check_ready_for_analysis()

    def check_ready_for_analysis(self):
        """Check if ready for analysis and enable/disable button"""
        ready = (
            self.blast_file.get()
            and self.combined_df is not None
            and self.probe_names_column.get()
            and self.sequence_column.get()
            and self.output_dir.get()
        )

        # ADDED - Actually update the button state
        self.analyze_btn.config(state="normal" if ready else "disabled")
        return ready

    def parse_blast_results(self, blast_text):
        """Parse BLAST result text and extract required information"""
        queries = re.split(r"Query #\d+: ", blast_text)[1:]
        records = []

        for query in queries:
            # Extract probe name
            probe_name_search = re.search(r"^(.+?)\s+Query ID:", query)
            probe_name = (
                probe_name_search.group(1).strip() if probe_name_search else None
            )

            # Count alignment sections
            alignment_headers = re.findall(r"^>([^\n]+)", query, re.MULTILINE)
            num_hits = len(alignment_headers)

            # Extract unique hit name if only one hit
            unique_hit_name = None
            if num_hits == 1:
                unique_hit_name = alignment_headers[0].strip()

            # Extract start and end positions
            start_end_search = re.search(r"Range 1: (\d+) to (\d+)", query)
            start = int(start_end_search.group(1)) if start_end_search else None
            end = int(start_end_search.group(2)) if start_end_search else None

            # Extract percentage identity
            perc_identity_search = re.search(r"Identities:\s*\d+/\d+\((\d+)%\)", query)
            perc_identity = (
                int(perc_identity_search.group(1)) if perc_identity_search else None
            )

            # Extract probe sequence
            query_seqs = re.findall(r"Query\s+\d+\s+([A-Z]+)\s+\d+", query)
            probe_seq = max(query_seqs, key=len) if query_seqs else None

            records.append(
                {
                    "ProbeName": probe_name,
                    "ProbeSequence": probe_seq,
                    "PercentAlignment": perc_identity,
                    "NumberOfHits": num_hits,
                    "UniqueHitName": unique_hit_name,
                    "Start": start,
                    "End": end,
                }
            )

        return pd.DataFrame(records)

    def run_analysis(self):
        """Run the BLAST analysis"""
        if not self.check_ready_for_analysis():
            messagebox.showerror("Error", "Please complete all required fields")
            return

        try:
            self.progress.start()
            self.analyze_btn.config(state="disabled")

            # Parse BLAST results
            self.log_message("Parsing BLAST results...")
            with open(self.blast_file.get(), "r") as f:
                blast_text = f.read()

            blast_df = self.parse_blast_results(blast_text)
            self.log_message(f"Parsed {len(blast_df)} probe results from BLAST")

            # Merge with combined CSV data
            self.log_message("Merging BLAST results with probe data...")
            probe_col = self.probe_names_column.get()
            seq_col = self.sequence_column.get()

            merged_df = self.combined_df.merge(
                blast_df, left_on=probe_col, right_on="ProbeName", how="left"
            )

            self.log_message(f"Merged data: {len(merged_df)} total rows")

            # Apply filtering if requested
            if self.unique_hits_only.get():
                filtered_df = merged_df[merged_df["NumberOfHits"] == 1].copy()
                self.log_message(
                    f"Filtered for unique hits: {len(filtered_df)} rows remaining"
                )
            else:
                filtered_df = merged_df.copy()
                self.log_message("No filtering applied - keeping all results")

            # Remove duplicate ProbeName column
            if "ProbeName" in filtered_df.columns:
                filtered_df.drop(columns=["ProbeName"], inplace=True)

            # Save results
            output_dir = self.output_dir.get()

            # Save BLAST results
            blast_output = os.path.join(output_dir, "blast_results.csv")
            blast_df.to_csv(blast_output, index=False)
            self.log_message(f"BLAST results saved to: {blast_output}")

            # Save merged results
            if self.unique_hits_only.get():
                merged_output = os.path.join(
                    output_dir, "merged_results_unique_hits.csv"
                )
            else:
                merged_output = os.path.join(output_dir, "merged_results_all.csv")
            filtered_df.to_csv(merged_output, index=False)
            self.log_message(f"Merged results saved to: {merged_output}")

            # Generate summary
            self.log_message("\n=== ANALYSIS SUMMARY ===")
            self.log_message(f"Total probes processed: {len(blast_df)}")
            self.log_message(
                f"Probes with unique hits: {len(blast_df[blast_df['NumberOfHits'] == 1])}"
            )
            self.log_message(
                f"Probes with multiple hits: {len(blast_df[blast_df['NumberOfHits'] > 1])}"
            )
            self.log_message(
                f"Probes with no hits: {len(blast_df[blast_df['NumberOfHits'] == 0])}"
            )

            # Hit count distribution
            self.log_message("\n=== HIT COUNT DISTRIBUTION ===")
            hit_counts = blast_df["NumberOfHits"].value_counts().sort_index()
            for hits, count in hit_counts.items():
                self.log_message(f"Probes with {hits} hit(s): {count}")

            self.progress.stop()
            self.analyze_btn.config(state="normal")

            messagebox.showinfo(
                "Analysis Complete",
                f"Analysis completed successfully!\n\n"
                f"Results saved to:\n"
                f"• {os.path.basename(blast_output)}\n"
                f"• {os.path.basename(merged_output)}",
            )

        except Exception as e:
            self.progress.stop()
            self.analyze_btn.config(state="normal")
            error_msg = f"Error during analysis: {str(e)}"
            self.log_message(error_msg)
            messagebox.showerror("Analysis Error", error_msg)


def main():
    root = tk.Tk()
    app = SmFISHBlastAnalyzerGUI(root)

    # FIXED - Correct event bindings for combobox selection
    app.probe_combo.bind(
        "<<ComboboxSelected>>", lambda e: app.check_ready_for_analysis()
    )
    app.sequence_combo.bind(
        "<<ComboboxSelected>>", lambda e: app.check_ready_for_analysis()
    )

    root.mainloop()


if __name__ == "__main__":
    main()
