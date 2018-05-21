#! ${params.python3}

from Bio import SeqIO, SeqRecord, Seqs
from collections import defaultdict
import pickle
from ete3 import ncbi_taxonomy
# from ETE3_Utils import defaultdict_string, defaultdict_defaultdict

def get_clade_dict(tax_map):
    ncbi = ncbi_taxonomy.NCBITaxa()
    taxon_clade_mod = defaultdict(str)
    for k, v in tax_map.items():
        lineage = ncbi.get_lineage(v)
        lineage_string = ""
        for l in lineage:
            rank = ncbi.get_rank([l])[l]
            if rank != 'no rank':
                lineage_string = "{} {};".format(lineage_string, ncbi.get_taxid_translator([l])[l])
        taxon_clade_mod[k] = lineage_string
    return taxon_clade_mod


def parse_eggnog_map(mapping_file):
    local_eggnog_map = {}
    with open(mapping_file) as handle:
        for line in handle:
            if line.startswith('-'):
                continue
            else:
                line = line.split()
                local_eggnog_map[line[1]] = line[0]
    return local_eggnog_map


def commonprefix(taxon_paths):
    taxon_path_lists = [t.split('; ') for t in taxon_paths]
    s1 = min(taxon_path_lists)
    s2 = max(taxon_path_lists)
    for i, c in enumerate(s1):
        if c != s2[i]:
            return s1[:i]
    return s1


def remove_duplicate_seqs(fasta, taxon_clade_mod):
    ncbi = ncbi_taxonomy.NCBITaxa()
    seq_dict = defaultdict(list)
    seq_rec_dict = {}
    for rec in SeqIO.parse(fasta,'fasta'):
        seq_dict[str(rec.seq)].append(rec.id)
        seq_rec_dict[rec.id] = rec
    for k,v in seq_dict.items():
        if len(v) > 1:
            tax_strings = []
            seqrec = seq_rec_dict[v[0]]
            for id in v:
                tax_strings.append(taxon_clade_mod[id])
                taxon_clade_mod.pop(id)
                seq_rec_dict.pop(id)
            try:
                tax_path = commonprefix(tax_strings)
                taxon = tax_path[-1].replace(';', '')
                id = ncbi.get_name_translator([taxon])[taxon][0]
                taxon_clade_mod[id] = "; ".join(tax_path)
                seqrec.id = "{}.{}".format(id, taxon)
                seq_rec_dict[seqrec.id] = seqrec
            except:
                print("key error")

def write_class_map(eggnog_id_mapping, basename):
    with open("{}.class".format(basename_files), 'w') as class_map:
        class_map.write("%s\\t%s\\n" % (eggnog_id_mapping, basename))


def write_eggnog_map(eggnog_id_mapping, basename, fasta):
    tax_map = defaultdict(str)
    with open("{}_eggnog.map".format(basename), 'w') as eggnog_map:
        for rec in SeqIO.parse(fasta, 'fasta'):
            s = rec.id.split('.')
            eggnog_map.write("%s\\t%s\\n" % (s[1], eggnog_id_mapping))
            tax_map[rec.id] = int(s[0])
    return tax_map


def write_taxon_map(taxon_clade_mod, basename):
    with open("{}.taxid.map".format(basename), 'w') as tax_map:
        for k,v in taxon_clade_mod.items():
            tax_map.write("{}\t{}".format(k,v))


def write_mapping_files(eggnog_id, base_file, local_eggnog_map, fasta):
    basename_files = base_file.split('.')[1] if len(base_file.split('.')) > 1 else base_file
    write_class_map(local_eggnog_map[eggnog_id], basename_files)
    tax_map = write_eggnog_map(local_eggnog_map[eggnog_id], base_files, fasta)
    taxon_clade_mod = get_clade_dict(tax_map)
    write_taxon_map(taxon_clade_mod, base_files)


def main(eggnog_id, eggnog_map, base_file, fasta):
    local_eggnog_map = parse_eggnog_map(eggnog_map)
    write_mapping_files(eggnog_id, base_file, local_eggnog_map, fasta)


# if __name__ == '__main__':
#     main("${fasta.simpleName}", "${megan_eggnog_map}", "${fasta.baseName}", "${fasta}")
