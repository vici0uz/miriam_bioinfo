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
    "Mexico",
    "USA",
    "Venezuela",
    "Brazil",
    "Nicaragua",
    "Peru",
    "Puerto Rico",
    "Colombia",
    "Argentina",
]

GENOTYPE = {
    '1': 'I',
    '2': 'II',
    '3': 'III',
    '4': 'IV'
}

def get_features_obj(genbank_record): # 👌
    print(genbank_record)
    if genbank_record.features:
        return genbank_record.features
    return False

def get_settings(): # 👌
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

def _get_year(data): # 👌
    match = re.search(r'\b(?:2023|20[0-2][0-9]|19[70][0-9])\b', data)
    if match:
        fecha = match[0]
        return fecha
    return False

def _get_month(data):  # 👌
    match = re.search(r'(?:Jan(?:uary)?|Feb(?:uary)?|Mar(?:ch)?|Apr(?:il)?|May|Jun(?:e)?|Jul(?:y)?|Aug(?:ust)?|Sep(?:tember)?|Oct(?:ober)?|Nov(?:ember)?|Dec(?:ember)?)', data, re.IGNORECASE)
    # match = re.search(r'(?:Jan(?:uary)?|Feb(?:uary)?|Mar(?:ch)?|Apr(?:il)?|May|Jun(?:e)?|Jul(?:y)?|Aug(?:ust)?|Sep(?:tember)?|Oct(?:ober)?|Nov(?:ember)?|Dec(?:ember)?)', data)
    if match:
        return match.group(0)
    return False

def _get_day(data): # 👌
    match = re.search(r'(\d{1,2})', data, flags=re.IGNORECASE)
    if match:
        return match.group(0)
    return False

def _get_genotype(genbank_record): # 👌
    data = get_features_values(genbank_record, 'source', 'note')
    res = ''
    if data:
        match = re.search(r'(I|II|III|IV)', data, flags=re.IGNORECASE)
        if match:
            res = str(match[0]).upper()
        else:
            match = re.search(r'(1|2|3|4)', data, flags=re.IGNORECASE)
            if match:
                hold = match[0]
                if(GENOTYPE.get(hold)):
                    res = GENOTYPE[hold]
        return res
    return False

def get_serotype(data):
    match = re.search(r'(?:.*dengue virus )(?:type )?(\d)', data, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    return NO_DATA



def _get_gene(data):
    match = re.search(r'(polyprotein|partial sequence|envelope|complete)', data, flags=re.IGNORECASE)
    if match:
        print(match)
    return False

def get_gene(fasta_record, genbank_record):
    data = fasta_record
    res = _get_gene(data)
    if res:
        return res
    return False

def get_locus(gb_record):
    res = '%s, %s bp %s linear' % (gb_record.locus, gb_record.size, gb_record.molecule_type)
    return res


    
def get_features_values(genbank_record, feature_key, qualifier_key):
    features = get_features_obj(genbank_record)
    if features:
        for f in features:
            print(f)
            if f.key == feature_key:
                vals = f.qualifiers
                for v in vals:
                    if v.key == "/%s=" %(qualifier_key):
                        return v.value
    return False

def get_date(fasta_record, genbank_record): # PREGUNTAR SI LEER FECHA DE FASTA O GENBANK
    year = month = day = NO_DATA
    fecha_genbank = get_features_values(genbank_record, 'source', 'collection_date')
    if fecha_genbank:
        res_year = _get_year(fecha_genbank)
        if res_year:
            year = res_year
        res_month = _get_month(fecha_genbank)
        if res_month:
            month = res_month
        res_day = _get_day(fecha_genbank)
        if res_day:
            day = res_day

    else:
        fecha_fasta = fasta_record.long_name
        res_year = _get_year(fecha_fasta)
        if res_year:
            year = res_year
        res_month = _get_month(fecha_fasta)
        if res_month:
            month = res_month
    res = {
        'year': year, 
        'month': month,
        'day': day
        }
    return res

def _get_isolate(genbank_record):
    data = get_features_values(genbank_record, 'source', 'isolate')
    if data:
        return data
    return False

def get_isolate(genbank_record):
    res = _get_isolate(genbank_record)
    if res:
        res = res.replace('"', '')
        return res
    else:
        return NO_DATA


def get_genotype(genbank_record):
    res = _get_genotype(genbank_record)
    if res:
        return res
    else:
        return NO_DATA


def process_record(record_id, fasta_record, genbank_record):
    fasta_data = fasta_record.long_name
    res_date = get_date(fasta_record, genbank_record)
    serotype = get_serotype(fasta_data)
    genotype = get_genotype(genbank_record)
    locus = get_locus(genbank_record)
    isolate = get_isolate(genbank_record)
    data = {
        'accession': record_id,
        'version': record_id,
        'year': res_date['year'],
        'serotype': serotype,
        'month': res_date['month'],
        'locus': locus,
        'isolate': isolate,
        'day': res_date['day'],
        'main_administrative_division': NO_DATA,
        'city': NO_DATA,
        'genotype': genotype,
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
    # for index, col in enumerate(cols):
    #     worksheet.write_string(0, index, col)
    # for index, row in enumerate(data):
    #     worksheet.write_row((index +1), ['', '', '', '', country, '', '', '', row['locus'], '', row['serotype'], '', '', '','', ''])

    # workbook.close()

def get_pub(res_id):
    data = Entrez.elink(db='pubmed', dbfrom='nucleotide', id=res_id, cmd='acheck').read()
    print(data)

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
                if index == 10: # CORTE DEBUG
                    break
                data_fetched = Entrez.efetch(db="nucleotide", id=res, rettype="fasta", retmode="text").read()
                data_summary_fetched = Entrez.efetch(db="nucleotide", id=res, rettype="gb", retmode="text").read()
                # pubmed = get_pub(res)                
                data_io = StringIO(data_fetched)
                data_gb = StringIO(data_summary_fetched)

                records = SeqIO.parse(data_io, 'fasta')
                genbank_record = GenBank.read(data_gb)
                for rec in records:
                    file_name = '%s/%s.fasta' % (file_path, rec.id)
                    with open(file_name, 'w') as file:
                        file.write(str(data_fetched))
                        file.close()
                        fasta_sequence = Fasta(file_name) # 👌
                        if make_sheets:
                            print('hace sheet')
                            for fasta_record in fasta_sequence:
                                # fasta_name = record.name # 👈
                                record_id = rec.id
                                data = process_record(record_id, fasta_record, genbank_record,)
                                data['country'] = country
                                list_data.append(data)
                                
                                # print(data)
                                
            make_sheet(list_data, country)

make_sheets = True  if (sys.argv[1] == 'make_sheets' and len(sys.argv)) > 1 else False

download_data(make_sheets=make_sheets)
