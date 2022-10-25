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

from pyfaidx import Fasta # ELIMINAR DEPENDENCIA DE pyfaidx

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
    # "Mexico",
    # "USA",
    # "Venezuela",
    # "Brazil",
    # "Nicaragua",
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
    if genbank_record.features:
        return genbank_record.features
    return False

def get_reference_obj(genbank_record):
    if genbank_record.references:
        return genbank_record.references
    return False

def get_settings(): # 👌
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    print(os.path.join(__location__, "settings.ini"))
    with open(os.path.join(__location__, "settings.ini")) as file:
        file_data = file.read()
        try:
            data = ast.literal_eval(file_data)
            return data
        except:
            print("Could not read settings")

def _get_year(data): # 👌
    match = re.search(r'\b(?:2023|20[0-2][0-9]|19[70][0-9])\b', data) # Matchea con años YYYY entre 2023 y 1970
    if match:
        fecha = match[0]
        return fecha
    return False

def _get_month(data):  # 👌
    # matchea con los meses en Ingles
    match = re.search(r'(?:Jan(?:uary)?|Feb(?:uary)?|Mar(?:ch)?|Apr(?:il)?|May|Jun(?:e)?|Jul(?:y)?|Aug(?:ust)?|Sep(?:tember)?|Oct(?:ober)?|Nov(?:ember)?|Dec(?:ember)?)', data, re.IGNORECASE)
    if match:
        return match.group(0)
    return False

def _get_day(data): # 👌
    # matchea con un numero entre de dos cifras (dia)
    match = re.search(r'(\d{1,2})', data, flags=re.IGNORECASE)
    if match:
        return match.group(0)
    return False

def _get_genotype(genbank_record): # 👌
    data = get_features_values(genbank_record, 'source', 'note')
    res = ''
    if data:
        # Match con la numeracion romana - Primera pasada
        match = re.search(r'(I|II|III|IV)', data, flags=re.IGNORECASE)
        if match:
            res = str(match[0]).upper()
        else:
            # Match con numeros arabigos - Segunda Pasada
            match = re.search(r'(1|2|3|4)', data, flags=re.IGNORECASE)
            if match:
                hold = match[0]
                if(GENOTYPE.get(hold)):
                    res = GENOTYPE[hold]
        return res
    return False

def _get_city(genbank_record): # 👌
    data = get_features_values(genbank_record, 'source', 'country')
    if data:
        data = data.replace('"', '')
        if data.find(':'):
            target_str = data.split(':')
            if (len(target_str)>1):
                res = target_str[1]
                return res
    return False

def get_city(genbank_record): # 👌
    res = _get_city(genbank_record)
    if res:
        return res
    else:
        return NO_DATA


def get_source(genbank_record): # 👌
    data = _get_source(genbank_record)
    if data:
        return data
    else:
        return NO_DATA


def get_serotype(data): # 👌
    # Match con "dengue virus type n", solo captura la n
    match = re.search(r'(?:.*dengue virus )(?:type )?(\d)', data, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    return NO_DATA


def _get_gene(data):
    # match con todas las aparaciones de estas cadenas
    match = re.findall(r'(envelope|partial|polyprotein|complete (?!cds))', data, flags=re.IGNORECASE)
    if match:
        res = ''
        for m in match:
            if m == 'complete':
                res += 'Complete genome'
                break
            elif m == 'envelope':
                res += 'Envelope protein'
                for p in match:
                    if p == 'partial':
                        res += ', partial'
                break
            elif m == 'polyprotein':
                res += 'Polyprotein'
                for p in match:
                    if p == 'envelope':
                        res = 'Envelope protein'
                        break
                    elif p == 'partial':
                        res += ', partial'
                break
        return res
    return False


def get_gene(fasta_record):
    data = fasta_record
    res = _get_gene(data)
    if res:
        return res
    else:
        return NO_DATA


def get_locus(genbank_record):
    # Interpolacion de cadenas
    res = '%s, %s bp %s linear' % (genbank_record.locus, genbank_record.size, genbank_record.molecule_type)
    return res


def get_features_values(genbank_record, feature_key, qualifier_key):
    features = get_features_obj(genbank_record)
    if features:
        # features es una lista
        for f in features:
            if f.key == feature_key:
                vals = f.qualifiers
                for v in vals:
                    if v.key == "/%s=" %(qualifier_key):
                        return v.value
    return False

def get_reference_values(genbank_record, reference_key):
    references = get_reference_obj(genbank_record)
    if references:
        for ref in references:
            if getattr(ref, reference_key):
                data = getattr(ref, reference_key)
                return data
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


def _get_source(genbank_record):
    source = get_reference_values(genbank_record, 'journal')
    pubmed = get_pubmed(genbank_record)
    if pubmed:
        return pubmed
    else:
        if source:
            match = re.search(r'(unpublished|direct submission)', source, flags=re.IGNORECASE)
            if match:
                res = match[0]
                return res
            else:
                return 'Direct submission'
    return False     



def process_record(record_id, fasta_record, genbank_record):
    # print(genbank_record)
    fasta_data = fasta_record.long_name
    print(fasta_data)
    res_date = get_date(fasta_record, genbank_record)
    serotype = get_serotype(fasta_data)
    genotype = get_genotype(genbank_record)
    locus = get_locus(genbank_record)
    isolate = get_isolate(genbank_record)
    city = get_city(genbank_record)
    source = get_source(genbank_record)
    gene = get_gene(fasta_data)
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
        'imported_case_sample_processing_country': NO_DATA,
        'observation': NO_DATA,
        'ref_id': NO_DATA,
        'city': city,
        'genotype': genotype,
        'source': source,
        'gene': gene,
    }
    return data

def make_sheet(data, country):
    headers = ["Accession n°", "Year", "Month", "Day", "Country of origin", "Main administrative division", "City", "Imported case/sample processing country","Locus", "Gene", "Serotype", "Genotype", "Source", "Isolate", "Observation", "Ref. ID"]
    cols = ["accession", "year", "month", "day", "country", "main_administrative_division", "city", "imported_case_sample_processing_country","locus", "gene", "serotype", "genotype", "source", "isolate", "observation", "ref_id"]
    df = pd.DataFrame(data, columns=cols)
    with pd.ExcelWriter('results/%s/data.xlsx' %(country,)) as writer:
        df.to_excel(writer, 'Sheet1', index=False, header=True)

        workbook = writer.book
        worksheet = writer.sheets['Sheet1']
        props = {
            'bold': True,
            'text_wrap': True,
            'valign': 'top'
        }
        header_format = workbook.add_format(props)
        for col_num, value in enumerate(df.columns.values):
            worksheet.write(0, col_num, headers[col_num], header_format)
        # writer.close()

def _get_pubmed(genbank_record):
    res = get_reference_values(genbank_record, 'pubmed_id')
    if res:
        return res
    return False


def get_pubmed(genbank_record):
    res = _get_pubmed(genbank_record)
    if res:
        return 'https://pubmed.ncbi.nlm.nih.gov/%s/' % ( res,)
    else:
        return 

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
            for index, res in enumerate(result):
                # if index == 10: # CORTE DEBUG
                #     break
                data_fetched = Entrez.efetch(db="nucleotide", id=res, rettype="fasta", retmode="text").read()
                data_summary_fetched = Entrez.efetch(db="nucleotide", id=res, rettype="gb", retmode="text").read()
                
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
                            for fasta_record in fasta_sequence:
                                record_id = rec.id
                                data = process_record(record_id, fasta_record, genbank_record,)
                                data['country'] = country
                                list_data.append(data)
            make_sheet(list_data, country)
# Lee la entrada del comando
make_sheets = True  if (sys.argv[1] == 'make_sheets' and len(sys.argv)) > 1 else False
# Invoca la función inicial
download_data(make_sheets=make_sheets)
