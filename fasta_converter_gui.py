# fasta_converter_gui.py

import pandas as pd
import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk


class FastaConverterGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("CSV/Excel to FASTA Converter")
        self.root.geometry("600x400")
        self.root.resizable(True, True)

        # Variables
        self.input_file = tk.StringVar()
        self.header_column = tk.StringVar()
        self.sequence_column = tk.StringVar()
        self.output_file = tk.StringVar()
        self.df = None
        self.columns = []

        self.create_widgets()

    def create_widgets(self):
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)

        # File selection section
        ttk.Label(main_frame, text="1. Select Input File:").grid(
            row=0, column=0, sticky=tk.W, pady=(0, 5)
        )

        file_frame = ttk.Frame(main_frame)
        file_frame.grid(
            row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 15)
        )
        file_frame.columnconfigure(0, weight=1)

        self.file_entry = ttk.Entry(
            file_frame, textvariable=self.input_file, state="readonly"
        )
        self.file_entry.grid(row=0, column=0, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(file_frame, text="Browse", command=self.browse_file).grid(
            row=0, column=1
        )

        # Column selection section
        ttk.Label(main_frame, text="2. Select Columns:").grid(
            row=2, column=0, sticky=tk.W, pady=(0, 5)
        )

        # Header column selection
        ttk.Label(main_frame, text="FASTA Header Column:").grid(
            row=3, column=0, sticky=tk.W, padx=(20, 0)
        )
        self.header_combo = ttk.Combobox(
            main_frame, textvariable=self.header_column, state="readonly"
        )
        self.header_combo.grid(row=3, column=1, sticky=(tk.W, tk.E), pady=2)

        # Sequence column selection
        ttk.Label(main_frame, text="Sequence Column:").grid(
            row=4, column=0, sticky=tk.W, padx=(20, 0)
        )
        self.sequence_combo = ttk.Combobox(
            main_frame, textvariable=self.sequence_column, state="readonly"
        )
        self.sequence_combo.grid(row=4, column=1, sticky=(tk.W, tk.E), pady=2)

        # Output file section
        ttk.Label(main_frame, text="3. Output File:").grid(
            row=5, column=0, sticky=tk.W, pady=(15, 5)
        )

        output_frame = ttk.Frame(main_frame)
        output_frame.grid(
            row=6, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 15)
        )
        output_frame.columnconfigure(0, weight=1)

        self.output_entry = ttk.Entry(output_frame, textvariable=self.output_file)
        self.output_entry.grid(row=0, column=0, sticky=(tk.W, tk.E), padx=(0, 5))

        ttk.Button(output_frame, text="Browse", command=self.browse_output).grid(
            row=0, column=1
        )

        # Convert button
        self.convert_btn = ttk.Button(
            main_frame,
            text="Convert to FASTA",
            command=self.convert_to_fasta,
            state="disabled",
        )
        self.convert_btn.grid(row=7, column=0, columnspan=2, pady=20)

        # Progress bar
        self.progress = ttk.Progressbar(main_frame, mode="indeterminate")
        self.progress.grid(
            row=8, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10)
        )

        # Status text
        self.status_text = tk.Text(main_frame, height=8, wrap=tk.WORD)
        self.status_text.grid(
            row=9, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(0, 10)
        )

        # Scrollbar for status text
        scrollbar = ttk.Scrollbar(
            main_frame, orient="vertical", command=self.status_text.yview
        )
        scrollbar.grid(row=9, column=2, sticky=(tk.N, tk.S))
        self.status_text.configure(yscrollcommand=scrollbar.set)

        # Configure row weight for text area
        main_frame.rowconfigure(9, weight=1)

    def log_message(self, message):
        """Add message to status text"""
        self.status_text.insert(tk.END, message + "\n")
        self.status_text.see(tk.END)
        self.root.update_idletasks()

    def browse_file(self):
        """Browse and select input file"""
        filename = filedialog.askopenfilename(
            title="Select CSV or Excel file",
            filetypes=[
                ("CSV files", "*.csv"),
                ("Excel files", "*.xlsx"),
                ("All files", "*.*"),
            ],
        )

        if filename:
            self.input_file.set(filename)
            self.load_file_columns()

            # Auto-generate output filename
            base_name = os.path.splitext(filename)[0]
            output_name = base_name + ".fasta"
            self.output_file.set(output_name)

            self.log_message(f"Selected input file: {filename}")

    def load_file_columns(self):
        """Load and populate column options from the selected file"""
        try:
            filename = self.input_file.get()

            # Read file to get column names
            if filename.endswith(".xlsx"):
                self.df = pd.read_excel(filename)
            elif filename.endswith(".csv"):
                self.df = pd.read_csv(filename)
            else:
                raise ValueError("File must be .xlsx or .csv")

            self.columns = list(self.df.columns)

            # Update comboboxes
            self.header_combo["values"] = self.columns
            self.sequence_combo["values"] = self.columns

            # Auto-select common column names if they exist
            for col in self.columns:
                if any(
                    keyword in col.lower()
                    for keyword in ["name", "id", "probe", "header"]
                ):
                    self.header_column.set(col)
                    break

            for col in self.columns:
                if any(
                    keyword in col.lower()
                    for keyword in ["seq", "sequence", "dna", "rna"]
                ):
                    self.sequence_column.set(col)
                    break

            self.log_message(
                f"Loaded {len(self.df)} rows with columns: {', '.join(self.columns)}"
            )
            self.check_ready_to_convert()

        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {str(e)}")
            self.log_message(f"Error loading file: {str(e)}")

    def browse_output(self):
        """Browse and select output file location"""
        filename = filedialog.asksaveasfilename(
            title="Save FASTA file as",
            defaultextension=".fasta",
            filetypes=[
                ("FASTA files", "*.fasta"),
                ("FASTA files", "*.fa"),
                ("Text files", "*.txt"),
                ("All files", "*.*"),
            ],
        )

        if filename:
            self.output_file.set(filename)
            self.check_ready_to_convert()

    def check_ready_to_convert(self):
        """Check if all required fields are filled and enable convert button"""
        if (
            self.input_file.get()
            and self.header_column.get()
            and self.sequence_column.get()
            and self.output_file.get()
        ):
            self.convert_btn.config(state="normal")
        else:
            self.convert_btn.config(state="disabled")

    def convert_to_fasta(self):
        """Convert the selected data to FASTA format"""
        try:
            self.progress.start()
            self.convert_btn.config(state="disabled")

            # Get selected columns
            header_col = self.header_column.get()
            seq_col = self.sequence_column.get()
            output_path = self.output_file.get()

            # Validate columns exist
            if header_col not in self.df.columns:
                raise ValueError(f"Header column '{header_col}' not found in file")
            if seq_col not in self.df.columns:
                raise ValueError(f"Sequence column '{seq_col}' not found in file")

            self.log_message(f"Converting {len(self.df)} sequences...")
            self.log_message(f"Header column: {header_col}")
            self.log_message(f"Sequence column: {seq_col}")

            # Write FASTA file
            with open(output_path, "w") as fasta_file:
                for index, row in self.df.iterrows():
                    header = str(row[header_col]).strip()
                    sequence = str(row[seq_col]).strip()

                    # Skip empty rows
                    if pd.isna(row[header_col]) or pd.isna(row[seq_col]):
                        continue

                    fasta_file.write(f">{header}\n")
                    fasta_file.write(f"{sequence}\n")

            self.progress.stop()
            self.convert_btn.config(state="normal")

            success_msg = f"Success! FASTA file saved as: {output_path}"
            self.log_message(success_msg)
            messagebox.showinfo("Conversion Complete", success_msg)

        except Exception as e:
            self.progress.stop()
            self.convert_btn.config(state="normal")
            error_msg = f"Error during conversion: {str(e)}"
            self.log_message(error_msg)
            messagebox.showerror("Conversion Error", error_msg)


def main():
    root = tk.Tk()
    app = FastaConverterGUI(root)

    # Bind combobox selection events to check if ready to convert
    app.header_combo.bind(
        "<<ComboboxSelected>>", lambda e: app.check_ready_to_convert()
    )
    app.sequence_combo.bind(
        "<<ComboboxSelected>>", lambda e: app.check_ready_to_convert()
    )

    root.mainloop()


if __name__ == "__main__":
    main()
