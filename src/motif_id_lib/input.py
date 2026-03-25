
import inquirer, re, csv
from termcolor import colored

class Sequence:
    """
    Parses sequence file and stores contents as list[header, sequence]
    """

    def __init__(self,filename) -> None:
        self.filename = filename
        self.sequence = []

    def genBank_parse(self) -> None:
        """
        Parses GenBank file format
        """
        with open(self.filename) as file:
            header = ""
            seq = ""
            for line in file:
                # Combine header lines to match FASTA header
                match1 = re.match(r'^DEFINITION\s(.*)$', line)
                if match1: header = match1.group(1)
                match2 = re.match(r'^VERSION\s*(.*)$', line)
                if match2: header = match2.group(1) + header 

                # Combine Sequence lines
                if "ORIGIN" in line:
                    lines = file.read()
                    seq_matches = re.findall(r'[a-z]+', lines)
                    seq = "".join(seq_matches).upper()
            self.sequence += header, seq

    def fasta_parse(self) -> None:
        """
        Parses FASTA file format
        """
        with open(self.filename) as file:
            header = file.readline().strip()
            seq = ""
            for line in file:
                seq += line.strip()
            self.sequence += header, seq

    def load_sequence(self) -> None:
        """
        Reads sequence file and parses header, sequence. 
        """
        seq_file = self.filename
        if ".gb" in seq_file:
            self.genBank_parse()
            
        elif ".fa" in seq_file:
            self.fasta_parse()
        

class Enzymes:
    """
    Loads restriction enzymes list from database file, filters list based on user input, then stores the final list in "filtered".
    """

    def __init__(self, database) -> None:
        self.db = database
        self.usr_list = []
        self.filtered = {}


    def app_header(self) -> None:
        """
        Prints TUI header in the terminal.
        """
        print(colored("\nREcut: Plasmid Sequence Cutting Tool\n", 'cyan', on_color='on_dark_grey', attrs=[ 'blink']))

    def interface(self, enzymes, interface) -> None:
        """
        Produces checkbox for user to filter enzyme dict. (default set manually)
        """""
        if interface:
            self.app_header()
            usr_enzymes = [
                inquirer.Checkbox('enzymes',
                                    message="Choose your enzymes below (use the arrows and space bar to select):",
                                    choices= enzymes
                                    ),
                ]
            answers = inquirer.prompt(usr_enzymes)
            self.usr_list = answers['enzymes']
        else:
            self.usr_list = enzymes

    def filter_enzymes(self) -> None:
        with open(self.db, "r") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if row["enzyme"] in self.usr_list:
                    enz = row["enzyme"]
                    motif = row["motif"]
                    cut = row["cutInfo"]
                    self.filtered[enz] = [motif, cut]
