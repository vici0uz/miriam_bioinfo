from pathlib import Path
from io import StringIO
import ast
import os
from re import S

from Bio import Entrez
from Bio import SeqIO

COUNTRYS = [
    "Brazil",
    "Venezuela",
    "Mexico",
    "Nicaragua",
    "Peru",
    "Puerto Rico",
    "Colombia",
    "USA",
    "Argentina",
]
def get_settings():
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    print(os.path.join(__location__, "settings.ini"))
    with open(os.path.join(__location__, "settings.ini")) as file:
        file_data = file.read()
        print(file_data)
        try:
            data = ast.literal_eval(file_data)
            return data
        except:
            print("Could not read settings")

def download_data():
    data = get_settings()
    if data:
        Entrez.email = data['email']
        Entrez.api_key = data['api_key']

        for country in COUNTRYS:
            file_path = 'results/%s' % country
            Path(file_path).mkdir(parents=True, exist_ok=True)
            search = Entrez.esearch(db='nucleotide', retmax=99999, term='Dengue virus %s' % (country))
            
            result = Entrez.read(Entrez.esearch(db="nucleotide", retmax=99999, term="Dengue virus %s"%(country,)))["IdList"]
            print(len(result))

            for res in result:
                data_fetched = Entrez.efetch(db="nucleotide", id=res, rettype="fasta", retmode="text").read()
                data_io = StringIO(data_fetched)
                records = SeqIO.parse(data_io, 'fasta') 
                for rec in records:
                    print(rec.id)
                    file_name = '%s/%s.fasta' % (file_path, rec.id)
                    with open(file_name, 'w') as file:
                        file.write(str(data_fetched))
                        file.close()

download_data()
