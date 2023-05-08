import math
import sys
from pathlib import Path
from io import StringIO
import ast
import re
import os
import pandas as pd
import xlsxwriter
import argparse

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

search_term = ''
tax_id = 0
db = ''
organism= ''
listado_paises = []
exclude = ''
sheet = False

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
    match = re.search(r'\b([12]\d{3})\b', data) #match con años entre 1000 y 2999 # TODO: Verificar que no sea un año negativo, matchear solo con 1970- fecha actual
    # match = re.search(r'\b(?:2023|20[0-2][0-9]|19[70][0-9])\b', data) # Matchea con años YYYY entre 2023 y 1970
    print(data)
    print(match)
    if match:
        fecha = match[0]
        return fecha
    return False

def _get_month(data):  # 👌
    # matchea con los meses en Inglés
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
    # match con todas las apariciones de estas cadenas
    match = re.findall(r'(envelope|partial|complete\sgenome|polyprotein)(?!cds)', data, flags=re.IGNORECASE)
    print(match)
    if match:
        res = ''
        for m in match:
            if m == 'complete genome':
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
        
        print(res)
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

def check_country(genbank_record, country):
    data = get_features_values(genbank_record, 'source', 'country')
    location_data = get_features_values(genbank_record, 'source', 'location')
    res = ''
    if data:
        data = data.replace('"', '')
        if data.find(':'):
            target_str = data.split(':')
            if (len(target_str)>=1):
                res = target_str[0]
        else:
            res =  data
    elif location_data:
        data = data.replace('"', '')
        if data.find(':'):
            target_str = data.split(':')
            if (len(target_str)>=1):
                res = target_str[0]
        else:
            res =  data
    if(res.lower() == country.lower()):
        return True
    return False

def download_data(make_sheets=False):
    fasta_xd = ''
    data = get_settings()
    if data:
        Entrez.email = data['email']
        Entrez.api_key = data['api_key']
        
        Path('results').mkdir(parents=True, exist_ok=True)
        with open('results/multi.fasta','w') as main_file:
            main_file.write('')
        for country in listado_paises:
            file_path = 'results/%s' % country
            list_data = []
            Path(file_path).mkdir(parents=True, exist_ok=True)

            term = f'({search_term}) AND {country} AND "{organism}"[porgn:__txid{tax_id}] AND"{min_length}:{max_length}[Sequence Length])'
            print(term)
            if exclude:
                term += f' NOT {exclude}'
            search = Entrez.esearch(db=f'{db}', retmax=99999, term=term)

            # search = Entrez.esearch(db='nucleotide', retmax=99999, term='(Dengue virus) AND %s) AND "Dengue virus"[porgn:__txid12637]' % (country,))
            # search = Entrez.esearch(db='nucleotide', retmax=10, term='(Dengue virus) AND %s) AND "Dengue virus"[porgn:__txid12637]' % (country,))
            # print(Entrez.read(search))
           
            result = Entrez.read(search)["IdList"]
            print(len(result))
            for index, res in enumerate(result):
                if index == 10: # CORTE DEBUG
                    break
                data_fetched = Entrez.efetch(db=f"{db}", id=res, rettype="fasta", retmode="text").read()
                data_summary_fetched = Entrez.efetch(db=f"{db}", id=res, rettype="gb", retmode="text").read()
                
                data_io = StringIO(data_fetched)
                data_gb = StringIO(data_summary_fetched)

                records = SeqIO.parse(data_io, 'fasta')
                
                try:
                    genbank_record = GenBank.read(data_gb)
                    for rec in records:
                        if check_country(genbank_record, country):
                            file_name = '%s/%s.fasta' % (file_path, rec.id)
                            with open(file_name, 'w') as file:
                                fasta_str = str(data_fetched)
                                if fasta_str.find(country) != -1:
                                    pass
                                else:
                                    fasta_lines = fasta_str.split('\n')
                                    words = fasta_lines[0].split()
                                    words.append(f'{country}')
                                    fasta_lines[0] = ' '.join(words)
                                    new_fasta_str = '\n'.join(fasta_lines)
                                    fasta_str = new_fasta_str
                                with open('results/multi.fasta','a') as main_file:
                                    main_file.write(str(fasta_str))
                                    
                                file.write(str(data_fetched))
                                file.close()
                                fasta_sequence = Fasta(file_name) # 👌
                                if sheet:
                                    for fasta_record in fasta_sequence:
                                        record_id = rec.id
                                        data = process_record(record_id, fasta_record, genbank_record,)
                                        data['country'] = country
                                        # x = check_country(genbank_record,)
                                        
                                        list_data.append(data)
                except:
                    print('error')
            if sheet:
                print(list_data)
                make_sheet(list_data, country)

# Lee la entrada del comando

# make_sheets = True  if (sys.argv[1] == 'make_sheets' and len(sys.argv)) > 1 else False
# Invoca la función inicial
# download_data(make_sheets=make_sheets)
parser = argparse.ArgumentParser(description="Args")
parser.add_argument('--sheet', type=bool, help="Export excel datasheet (Only for dengue virus)")
parser.add_argument('--multi', type=bool, help="Export to multifasta")
parser.add_argument('--tax_id', type=int, help="NCBI Taxonomical id")
parser.add_argument('--db', type=str, help="Db search", default='nucleotide')
parser.add_argument('--organism', type=str, help="Organism name")
parser.add_argument('--search_term', type=str, help="Search term in example: Dengue virus")
parser.add_argument('--min_length', type=int, help="Minimum length", default=100)
parser.add_argument('--max_length', type=int, help='Maximum length', default=1000000000)
parser.add_argument('--country_list', type=str, help="A list of country")
parser.add_argument('--exclude', type=str, help="Exclude word")
args = parser.parse_args()

if args.sheet:
    sheet = True
if args.tax_id:
    print(args.tax_id)
    tax_id = args.tax_id
else:
    raise ValueError("Tax id required")

if args.search_term:
    search_term = args.search_term
    print(args.search_term)
else:
    raise ValueError("Search term required")

if args.organism:
    organism = args.organism
else:
    raise ValueError("Organism required")

if args.country_list:
    print(type(args.country_list))
    print(len(args.country_list))
    print(args.country_list)
    listado_paises = args.country_list.split(",")
    print(listado_paises)
    print(len(listado_paises))

if args.exclude:
    exclude = args.exclude
db = args.db
min_length = args.min_length
max_length = args.max_length
download_data(make_sheets=False)