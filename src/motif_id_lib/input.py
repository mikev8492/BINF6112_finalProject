
import inquirer, sys, re

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
        # -------------------
        # Error Handling:
        # -------------------
        else:
            cprint("\nInvalid sequence filetype, please use a GenBank or FASTA file.\n", "red", "on_black",file=sys.stderr)
            sys.exit(1)
        

class Enzymes:
    """
    Loads restriction enzymes list from database file, filters list based on user input, then stores the final list in "filtered".
    """

    def __init__(self) -> None:
        self.renzymes = {}
        self.usr_list = []
        self.filtered = {}
        

    def load_REs(self):
        """
        Reads the enzyme database file and stores the list of enzymes as a dict to show the user. 
        """
        with open("src/database/enzymes.csv") as file:
            file.readline()
            for line in file:
                line = line.split(",")
                self.renzymes[line[0]] = [line[1],line[2].strip()]

    def interface(self) -> None:
        """
        Produces checkbox for user to filter enzyme dict. (default set manually)
        """""
        usr_enzymes = [
            inquirer.Checkbox('enzymes',
                                message="Choose your enzymes below (use the arrows and space bar to select):",
                                choices= self.renzymes.keys(),
                                default=['XhoI','BamHI','EcoRI']
                                ),
            ]
        answers = inquirer.prompt(usr_enzymes)
        self.usr_list = answers['enzymes']

    def filter_enzymes(self) -> None:
        for enzyme in self.usr_list:
            val = self.renzymes[enzyme]
            self.filtered[enzyme] = val