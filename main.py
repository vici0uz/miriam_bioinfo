import math
import sys
from pathlib import Path
from io import StringIO
import ast
import re
import os
import pandas as pd
import xlsxwriter

from Bio import Entrez, Seq, SeqIO, SeqRecord, GenBank

from pyfaidx import Fasta

NO_DATA = 'NA'

COLS = [
    "accession_n",
    "year",
    "mont",
    "day",
    "country_of_origin",
    "main_administrative_division",
    "city",
    "imported_case",
    "locus",
    "gene",
    "serotype",
    "genotype",
    "source",
    "isolate",
    "observation",
    "ref_id"
]

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

# SEROTYPES = {
#     '1': 'I',
#     '2': 'II',
#     '3': 'III',
#     '4': 'IV'
# }

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

def get_year(data, gb_record):
    match = re.search(r'\b(?:2023|20[0-2][0-9]|19[70][0-9])\b', data)
    if match:
        fecha = match[0]
        return fecha
    print(gb_record)
    return NO_DATA
    

def get_serotype(data):
    match = re.search(r'(?:.*dengue virus )(?:type )?(\d)', data, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    return NO_DATA

def get_month(data):
    match = re.search(r'(?:Jan(?:uary)?|Feb(?:uary)?|Mar(?:ch)?|Apr(?:il)?|May|Jun(?:e)?|Jul(?:y)?|Aug(?:ust)?|Sep(?:tember)?|Oct(?:ober)?|Nov(?:ember)?|Dec(?:ember)?)', data, re.IGNORECASE)
    # match = re.search(r'(?:Jan(?:uary)?|Feb(?:uary)?|Mar(?:ch)?|Apr(?:il)?|May|Jun(?:e)?|Jul(?:y)?|Aug(?:ust)?|Sep(?:tember)?|Oct(?:ober)?|Nov(?:ember)?|Dec(?:ember)?)', data)
    if match:
        return match.group(0)
    return NO_DATA

def get_gene(data):
    match = re.search(r'(polyprotein|partial sequence|envelope|complete)', data, flags=re.IGNORECASE)
    if match:
        print(match)
    return NO_DATA

def get_locus(gb_record):
    res = '%s, %s bp %s linear' % (gb_record.locus, gb_record.size, gb_record.molecule_type)
    return res
def get_features_obj(gb_record):
    if gb_record.features:
        return gb_record.features
    return False
    
def get_isolate(gb_record):
    res = {}
    if gb_record.features:
        features = gb_record.features
        for f in features:
            if f.key == 'source':
                vals = f.qualifiers
                for v in vals:
                    if v.key == '/isolate=':
                        res['isolate'] = str(v.value).replace('"', '')
                    if v.key == '/genotype=':
                        res['genotype'] = v.value
    return res

def process_record(fasta, gb_record, file_name):
    fasta_name = fasta.name
    fasta_data = fasta.long_name

    year = get_year(fasta_data, gb_record)
    serotype = get_serotype(fasta_data)
    month = get_month(fasta_data)
    gene = get_gene(fasta_data)
    locus = get_locus(gb_record)
    isolate = get_isolate(gb_record)['isolate'] 
    day = NO_DATA
    data = {
        'accession': file_name,
        'version': file_name,
        'year': year,
        'serotype': serotype,
        'month': month,
        'locus': locus,
        'isolate': isolate,
        'day': day,
        'main_administrative_division': NO_DATA,
        'city': NO_DATA
    }
    return data

def make_sheet(data, country):
    # workbook = xlsxwriter.Workbook('results/%s/data.xlsx' % country)
    # file_path = '' % country
    # worksheet = workbook.add_worksheet()
    print(data)
    headers = ["Accession n°", "Year", "Month", "Day", "Country of origin", "Main administrative division", "City", "Imported case/sample processing country","Locus", "Gene", "Serotype", "Genotype", "Source", "Isolate", "Observation", "Ref. ID"]

    cols = ["accession", "year", "month", "day", "country", "main_administrative_division", "city", "Imported case/sample processing country","locus", "gene", "serotype", "genotype", "source", "isolate", "observation", "ref_id"]
    df = pd.DataFrame(data, columns=cols)
    with pd.ExcelWriter('results/%s/data.xlsx' %(country,)) as writer:
        df.to_excel(writer, 'Sheet1', index=False)
    print(df)
    # for index, col in enumerate(cols):
    #     worksheet.write_string(0, index, col)
    # for index, row in enumerate(data):
    #     worksheet.write_row((index +1), ['', '', '', '', country, '', '', '', row['locus'], '', row['serotype'], '', '', '','', ''])

    # workbook.close()


def download_data(make_sheets=False):
    data = get_settings()
    if data:
        Entrez.email = data['email']
        Entrez.api_key = data['api_key']

        for country in COUNTRYS:
            file_path = 'results/%s' % country
            list_data = []
            Path(file_path).mkdir(parents=True, exist_ok=True)
            search = Entrez.esearch(db='nucleotide', retmax=99999, term='(Dengue virus) AND %s) AND "Dengue virus"[porgn:__txid12637]' % (country,))
            result = Entrez.read(search)["IdList"]
            data_holder = []
            for index, res in enumerate(result):
                if index == 10:
                    break
                data_fetched = Entrez.efetch(db="nucleotide", id=res, rettype="fasta", retmode="text").read()
                data_summary_fetched = Entrez.efetch(db="nucleotide", id=res, rettype="gb", retmode="text").read()
                

                data_io = StringIO(data_fetched)
                data_gb = StringIO(data_summary_fetched)

                records = SeqIO.parse(data_io, 'fasta')
                gb_record = GenBank.read(data_gb)
                print(gb_record.accession)
                for rec in records:
                    file_name = '%s/%s.fasta' % (file_path, rec.id)
                    with open(file_name, 'w') as file:
                        file.write(str(data_fetched))
                        file.close()
                        fasta_sequence = Fasta(file_name)
                        if make_sheets:
                            print('hace sheet')
                            for record in fasta_sequence:
                                # fasta_name = record.name # 👈
                                data = process_record(record, gb_record, rec.id)
                                data['country'] = country
                                list_data.append(data)
                                
                                # print(data)
                                
            make_sheet(list_data, country)

make_sheets = True  if (sys.argv[1] == 'make_sheets' and len(sys.argv)) > 1 else False

download_data(make_sheets=make_sheets)
